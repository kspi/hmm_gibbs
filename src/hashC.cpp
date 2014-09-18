#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector hashC(IntegerVector cs) {
    unsigned long hash = 5381;
    for (int i = 0; i < cs.size(); ++i) {
        unsigned long c = cs[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    IntegerVector ret(1);
    ret[0] = hash;
    return ret;
}
