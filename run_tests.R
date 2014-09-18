source("lib.R")
library(Rcpp)

sourceCpp("src/hashC.cpp")
hash <- function(s) {
    hashC(utf8ToInt(s))
}

runTest <- function(dir, chain_length, transition_prob, sigma) {
    stopifnot(ncol(transition_prob) == nrow(transition_prob))
    nstates <- nrow(transition_prob)
    means <- seq(0, length=nstates)
    if (file.exists(file.path(dir, "done"))) {
        message("Test ", dir, " already done.\n")
        return()
    }
    message("Running test ", dir, " with ", nstates, " states, ", chain_length, " length chain.")
    cat("Means: "); print(means)
    cat("Sigma: "); print(sigma)
    print(transition_prob)
    system(paste("rm -rf", dir))
    dir.create(dir, recursive=TRUE)
    message("Generating data.")
    generateData(dir, nstates, chain_length, transition_prob, means, sigma)
    message("Running sampler.")
    runAuto(dir)
    message("Processing result.")
    processResult(dir)
    warnings()
    system(paste("touch", file.path(dir, "done")))
    message("Done.\n")
}

runTestMulti <- function(dir, ...) {
    set.seed(hash(dir))
    multiSeeds <- round(runif(4, 0, 10000))
    for (seed in multiSeeds) {
        set.seed(seed)
        runTest(paste0(dir, "_seed", seed), ...)
    }
}

tmStays <- function(n, pStay=0.9) {
    diag(n) * pStay + (1 - diag(n)) * (1 - pStay) / (n - 1)
}

stayTest <- function(states, chainLength, sigma) {
    dir <- sprintf("tests/%dstates_%d_stays_sigma%.2f", states, chainLength, sigma)
    runTestMulti(dir, chainLength, tmStays(states), rep(sigma, states))
}

stayTest(2, 100, 0.1)
stayTest(2, 1000, 0.1)
stayTest(2, 10000, 0.1)
stayTest(2, 100, 0.25)
stayTest(2, 1000, 0.25)
stayTest(2, 10000, 0.25)
stayTest(2, 100, 0.5)
stayTest(2, 1000, 0.5)
stayTest(2, 10000, 0.5)

stayTest(3, 100, 0.1)
stayTest(3, 1000, 0.1)
stayTest(3, 10000, 0.1)
stayTest(3, 100, 0.25)
stayTest(3, 1000, 0.25)
stayTest(3, 10000, 0.25)
stayTest(3, 100, 0.5)
stayTest(3, 1000, 0.5)
stayTest(3, 10000, 0.5)

stayTest(4, 100, 0.1)
stayTest(4, 1000, 0.1)
stayTest(4, 10000, 0.1)
stayTest(4, 100, 0.25)
stayTest(4, 1000, 0.25)
stayTest(4, 10000, 0.25)
stayTest(4, 100, 0.5)
stayTest(4, 1000, 0.5)
stayTest(4, 10000, 0.5)


tmA <- rbind(c(0.6, 0.3, 0.1),
             c(0.1, 0.8, 0.1),
             c(0.1, 0.3, 0.6))

testA <- function(chainLength, sigma) {
    dir <- sprintf("tests/3states_%d_A_sigma%.2f", chainLength, sigma)
    runTestMulti(dir, chainLength, tmA, rep(sigma, 3))
}

testA(100, 0.1)
testA(1000, 0.1)
testA(10000, 0.1)
testA(100, 0.25)
testA(1000, 0.25)
testA(10000, 0.25)
testA(100, 0.5)
testA(1000, 0.5)
testA(10000, 0.5)
