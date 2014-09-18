CXX = clang++

#CXXFLAGS = -std=c++11 -g -O0 -Wall
#CXXFLAGS = -std=c++11 -O3 -Wall -funroll-all-loops -ffast-math -fomit-frame-pointer -finline-functions -fstrict-aliasing -static -march=core2
#CXXFLAGS = -std=c++11 -O0 -Wall -static -g -march=core2
CXXFLAGS = -std=c++11 -O3 -Wall -static -march=core2

CHAINS = 4

.PHONY: all
all: run_tests

.PHONY: run_tests
run_tests: hmm_gibbs run_tests.R src/hashC.cpp
	mkdir -p tests
	Rscript run_tests.R 2>&1 | tee tests/run_tests.log

hmm_gibbs: src/hmm_gibbs.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f hmm_gibbs
