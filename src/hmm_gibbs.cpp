#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <algorithm>
#include <limits>
#include <armadillo>

#include <fenv.h>
#include <signal.h>

using namespace std;
using namespace arma;

mt19937_64 GENERATOR;

void normalize_simplex(rowvec &x) {
    double sum = 0;
    for (auto &xi : x) {
        assert(xi >= 0);
        sum += xi;
    }

    if (sum < numeric_limits<double>::epsilon()) {
        x.fill(1.0 / x.size());
    } else {
        x /= sum;
    }
}

rowvec rand_dirichlet(const rowvec &alpha) {
    rowvec ret(alpha.size());
    for (unsigned int i = 0; i < ret.size(); ++i) {
        gamma_distribution<> g(alpha[i], 1);
        ret[i] = g(GENERATOR);
    }
    normalize_simplex(ret);
    return ret;
}

template <typename T>
T square(T x) {
    return x * x;
}

static const double INV_SQRT_2PI = 0.3989422804014327;
double normal_log_pdf(double x, double m, double s)
{
    double z = (x - m) / s;
    return log(INV_SQRT_2PI) - log(s) - 0.5 * z * z;
}

class logsumexp {
    private:
        list<double> xs;
        double max = -INFINITY;

    public:
        void add(double x) {
            if (x > max) {
                max = x;
            }
            xs.push_back(x);
        }

        double get_result() const {
            double sum = 0;
            for (const auto &x : xs) {
                sum += exp(x - max);
            }
            return max + log(sum);
        }
};

void normalize_log_simplex(rowvec &x) {
    logsumexp s;
    for (auto &xi : x) {
        s.add(xi);
    }
    x -= s.get_result();
}

int rand_discrete(const rowvec &p) {
    discrete_distribution<> d(p.begin(), p.end());
    return d(GENERATOR);
}

double rand_gamma(double shape, double rate) {
    // The std::gamma_distribution second parameter is scale = 1 / rate.
    gamma_distribution<> d(shape, 1 / rate);
    return d(GENERATOR);
}

double rand_normal(double mean, double sd) {
    normal_distribution<> d(mean, sd);
    return d(GENERATOR);
}

struct model {
    const struct data {
        int nstates;

        vector<rowvec> runs;

        data(const string &datafile) {
            ifstream f(datafile);
            assert(f.good());
            int nruns;
            f >> nstates >> nruns;
            runs.resize(nruns);
            for (auto &observed_states : runs) {
                int chain_length;
                f >> chain_length;
                observed_states.resize(chain_length);
                for (int i = 0; i < chain_length; ++i) {
                    f >> observed_states[i];
                }
            }
        }
    } data;

    int seed;
    int iteration = 0;
    ofstream resultf;

    rowvec initial_prob;
    mat transition_prob;
    rowvec means;
    rowvec sds;

    double ksi;
    double Rsqr;
    double kappa;
    double gamma_alpha = 2;
    double gamma_beta;
    double gamma_beta_g = 0.2;
    double gamma_beta_h;

    mat transition_prob_alpha;
    rowvec state_sum;
    rowvec state_count;
    rowvec sse;

    model(int seed, const string &datafile, const string &resultfile) :
        data(datafile),
        seed(seed),
        resultf(resultfile, ios::binary | ios::out)
    {
        initial_prob.zeros(data.nstates);
        transition_prob.zeros(data.nstates, data.nstates);
        means.zeros(data.nstates);
        sds.zeros(data.nstates);
        transition_prob_alpha.ones(data.nstates, data.nstates);
        state_sum.zeros(data.nstates);
        state_count.zeros(data.nstates);
        sse.zeros(data.nstates);

        double observed_min = 0;
        double observed_max = 0;
        for (const auto &ys : data.runs) {
            for (const auto y : ys) {
                if (y < observed_min) observed_min = y;
                if (y > observed_max) observed_max = y;
            }
        }

        ksi = (observed_min + observed_max) / 2;
        Rsqr = square(observed_max - observed_min);
        kappa = 1 / (100 * Rsqr);
        gamma_beta_h = 10 / Rsqr;

        // Random initial values
        initial_prob = rand_dirichlet(1 + state_count);
        gamma_beta = gamma_beta_g / gamma_beta_h;
        for (int i = 0; i < data.nstates; ++i) {
            means[i] = rand_normal(ksi, 1 / sqrt(kappa));
            sds[i] = 1 / sqrt(rand_gamma(gamma_alpha, gamma_beta));
            transition_prob.row(i) = rand_dirichlet(transition_prob_alpha.row(i));
        }
        means = sort(means);

        resultf << "seed\titeration\tparameter\trow\tcol\tvalue\n";
    }

    void flush() {
        resultf.flush();
    }

    void sweep() {
        ++iteration;
        run_chains();
        generate_params();
    }

    void generate_params() {
        for (int i = 0; i < data.nstates; ++i) {
            double var = square(sds[i]);
            double denom = state_count[i] + var * kappa;
            means[i] = rand_normal(
                    (state_sum[i] + ksi * var * kappa) / denom,
                    sds[i] / sqrt(denom));
        }

        rowvec sds_invsq(data.nstates);
        for (int i = 0; i < data.nstates; ++i) {
            sds_invsq[i] = rand_gamma(
                    gamma_alpha + state_count(i) / 2,
                    gamma_beta + sse(i) / 2);
        }
        sds = 1 / sqrt(sds_invsq);

        gamma_beta = rand_gamma(
                gamma_beta_g + data.nstates * gamma_alpha,
                gamma_beta_h + sum(sds_invsq));

        initial_prob = rand_dirichlet(1 + state_count);
        for (int i = 0; i < data.nstates; ++i) {
            transition_prob.row(i) = rand_dirichlet(transition_prob_alpha.row(i));
        }
    }

    void run_chains() {
        mat log_transition_prob(log(transition_prob));

        state_sum.zeros();
        state_count.zeros();
        sse.zeros();
        transition_prob_alpha.ones();

        for (const auto &observed_states : data.runs) {
            // log_likelihood[k, i] = log p(y_{k+1:n} | X_k = i)
            mat log_likelihood(observed_states.size(), data.nstates);
            log_likelihood.row(observed_states.size() - 1).fill(0);
            for (int k = observed_states.size() - 2; k >= 0; --k) {
                for (int i = 0; i < data.nstates; ++i) {
                    logsumexp s;
                    for (int j = 0; j < data.nstates; ++j) {
                        s.add(log_transition_prob(i, j)
                                + normal_log_pdf(observed_states(k + 1), means(j), sds(j))
                                + log_likelihood(k + 1, j));
                    }
                    log_likelihood(k, i) = s.get_result();
                }
            }

            rowvec log_initial = log(initial_prob);
            for (int i = 0; i < data.nstates; ++i) {
                log_initial(i) +=
                    normal_log_pdf(observed_states(0), means(i), sds(i))
                    + log_likelihood(0, i);
            }
            normalize_log_simplex(log_initial);
            int state = rand_discrete(exp(log_initial));

            state_sum(state) += observed_states(0);
            state_count(state) += 1;
            sse(state) += square(observed_states(0) - means(state));

            int prev_state = state;
            for (int k = 1; k < observed_states.size(); ++k) {
                rowvec log_p = log_transition_prob.row(prev_state);
                for (int i = 0; i < data.nstates; ++i) {
                    log_p(i) +=
                        normal_log_pdf(observed_states(k), means(i), sds(i))
                        + log_likelihood(k, i);
                }
                normalize_log_simplex(log_p);

                prev_state = state;
                state = rand_discrete(exp(log_p));

                state_sum(state) += observed_states(k);
                state_count(state) += 1;
                transition_prob_alpha(prev_state, state)++;

                sse(state) += square(observed_states(k) - means(state));
            }
        }
    }

    void print() {
        print_vector("initial_prob", initial_prob);
        print_vector("means", means);
        print_vector("sds", sds);
        print_value("beta", 0, 0, gamma_beta);
        print_matrix("transition_prob", transition_prob);
    }

    void premature_exit() {
        resultf.flush();
        cerr << "Numerical error at iteration " << iteration << ", exiting." << endl;
    }

    private:
    template<typename T>
    void print_value(const char *name, unsigned row, unsigned col, T value) {
        resultf
            << seed << "\t"
            << iteration << "\t"
            << name << "\t"
            << row << "\t"
            << col << "\t"
            << value << "\n";
    }

    template<typename T>
    void print_vector(const char *name, const Row<T> &x) {
        for (unsigned int i = 0; i < x.size(); ++i) {
            print_value(name, 0, i, x[i]);
        }
    }

    template<typename T>
    void print_matrix(const char *name, const Mat<T> &m) {
        for (unsigned int i = 0; i < m.n_rows; ++i) {
            for (unsigned int j = 0; j < m.n_cols; ++j) {
                print_value(name, i, j, m(i, j));
            }
        }
    }
};

string format_time(double t) {
    double seconds = fmod(t, 60);
    double minutes = fmod((t - seconds) / 60, 60);
    double hours = floor(t / 60 / 60);
    ostringstream s;
    s
        << fixed << setprecision(0)
        << hours << ":"
        << setfill('0') << setw(2)
        << minutes << ":"
        << setfill('0') << setw(2)
        << seconds;
    return s.str();
}

static model *MODEL = NULL;
void sigfpe_callback(int sig) {
    MODEL->premature_exit();
    exit(0);
}

int main(int argc, char **argv) {
    // Enable floating point exceptions
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // Faster iostream
    ios::sync_with_stdio(false);

    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " seed iterations datafile resultfile" << endl;
        return 1;
    }

    const int seed = atoi(argv[1]);
    const int iterations = atoi(argv[2]);
    const string datafile = argv[3];
    const string resultfile = argv[4];

    int burn_in = 0;
    int thin = 1;
    if (argc == 7) {
        burn_in = atoi(argv[5]);
        thin = atoi(argv[6]);
    }

    GENERATOR.seed(seed);

    model m(seed, datafile, resultfile);
    MODEL = &m;
    signal(SIGFPE, sigfpe_callback);

    wall_clock timer;
    timer.tic();
    double last_print = timer.toc();

    {
        if (!burn_in) m.print();
        for (int i = 1; i < iterations; ++i) {
            m.sweep();
            if (i > burn_in && i % thin == 0) {
                m.print();
            }

            if (i % 100 == 0 && (timer.toc() - last_print) > 10) {
                const double elapsed = timer.toc();
                last_print = elapsed;
                const double speed = i / elapsed;
                const double eta = (iterations - i) / speed;
                cerr
                    << setw(static_cast<int>(ceil(log10(iterations))) + 1)
                    << i << "/" << setw(0) << iterations << ", "
                    << setprecision(2) << fixed
                    << setw(6) << i * 100.0 / iterations << "%, "
                    << "elapsed: " << format_time(elapsed) << ", "
                    << setw(6) << speed << " iterations/s, "
                    << "ETA: " << format_time(eta)
                    << endl;
                m.flush();
            }
        }
    }

    cerr << "Done. " << iterations << " iterations. Elapsed: " << format_time(timer.toc()) << endl;

    return 0;
}
