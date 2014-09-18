library(modeest)
library(data.table)
library(foreach)
library(ggplot2)
library(reshape)
library(gridExtra)
library(ggthemes)
library(RHmm)


autocor <- function(x, lag=1) {
    cor(head(x, -lag), tail(x, -lag))
}
autocorrelations <- function(x, n=30) {
    data.table(delta=1:n,
               autocor=sapply(1:n, function(lag) autocor(x, lag)))
}

effectiveSampleSize <- function(x) {
    acs <- autocorrelations(x, n=100)
    length(x) / (1 + 2 * sum(acs$autocor))
}

thinFactor <- function(x) {
    ess <- effectiveSampleSize(x)
    max(1, factor <- length(x) / ess)
}

autoThin <- function(iteration, value) {
    indices <- seq(1, length(value), ceiling(thinFactor(value)))
    data.table(iteration=iteration[indices],
               value=value[indices])
}

modeEst <- function (x) {
    mlv(x, method='shorth')$M
}

modeHDI <- function (x, p=0.95) {
    n <- max(1, ceiling((1 - p) * length(x)))
    y <- sort(x)
    a <- head(y, n)
    b <- tail(y, n)
    i <- order(b - a)[1]
    list(mode=modeEst(x),
         lowerHDI=a[i],
         upperHDI=b[i])
}

readResult <- function(filename) {
    r <- fread(filename, stringsAsFactors=TRUE, showProgress=TRUE)
    r$row <- r$row + 1
    r$col <- r$col + 1
    for (s in unique(r$seed)) {
        ms <- r[parameter == 'means' & seed == s,
                median(value),
                by='col']$V1
        mapping <- rank(ms, ties.method="first")
        r[parameter != "beta" & seed == s, col := mapping[col]]
        r[parameter == 'transition_prob' & seed == s, row := mapping[row]]
    }
    r <- r[order(parameter, row, col, seed, iteration)]
    setkeyv(r, c("parameter", "row", "col"))
    r
}

transposeData <- function(data) {
    data[,
         local({
             l <- as.list(value)
             names(l) <- paste(parameter, row, col, sep='_')
             l
         }),
         by=c('seed', 'iteration')]
}

runGibbs <- function(seeds, iterations, dataFilename, resultFilename, burnIn=NULL, thinFactor=NULL) {
    additionalArgs <- if (!is.null(burnIn)) {
        paste(burnIn, thinFactor)
    } else {
        ""
    }
    system(paste0("parallel -u -j8 ./hmm_gibbs {} ",
                  iterations, " ",
                  dataFilename, " ",
                  resultFilename, ".{} ",
                  additionalArgs, " ",
                  " ::: ",
                  paste(seeds, collapse=" ")))
    system(paste0("./join_result ", resultFilename, " ", resultFilename, ".*"))
}

runAuto <- function(dir, samples=4000, nchains=8) {
    preSeeds <- round(runif(8, 1, 10000))
    seeds <- round(runif(nchains, 1, 10000))

    runGibbs(preSeeds, 2000, file.path(dir, "data.txt"), file.path(dir, "pre_result.txt"))
    message("Estimation run complete.")

    pre <- readResult(file.path(dir, "pre_result.txt"))

    thinningFactor <- ceiling(quantile(pre[,
                                       thinFactor(value),
                                       by=c('parameter', 'row', 'col', 'seed')]$V1,
    0.75))

    message("Thinning factor: ", thinningFactor)

    pairwiseApply <- function(l, f) {
        pairs <- combn(seq_along(l), 2)
        ns <- foreach(i=pairs[1,], j=pairs[2,], .combine=c) %do% {
            paste(i, j, sep=".")
        }
        vs <- foreach(i=pairs[1,], j=pairs[2,]) %do% {
            f(l[[i]], l[[j]])
        }
        names(vs) <- ns
        vs
    }

    multiKS <- function(x, y) {
        stopifnot(ncol(x) == ncol(y))
        foreach(i=1:ncol(x), .combine=min) %do% {
            ks.test(x[, i], y[, i])$p.value
        }
    }

    pre[, bin := cut(iteration, breaks=20, labels=FALSE)]

    estR <- function(x, seed) {
        min(data.table(y=x, seed)[, list(R=abs(mean(x) - mean(y))), by='seed']$R)
    }
    Rs <- pre[, list(iteration=tail(iteration, 1), R=estR(value, seed)), by=c('parameter', 'row', 'col', 'bin')]
    burnIn <- Rs[, list(burnIn=iteration[R > median(R)][1]), by=c('parameter', 'row', 'col')]

    burnInIteration <- 2 * (max(burnIn$burnIn, na.rm=TRUE) + 2)

    message("Burn-in iterations: ", burnInIteration)

    iterations <- burnInIteration + thinningFactor * samples %/% length(seeds)

    message("Target sample size: ", samples)
    message("MCMC iterations for each chain: ", iterations)

    runGibbs(seeds, iterations, file.path(dir, "data.txt"), file.path(dir, "result.txt"), burnInIteration, thinningFactor)
}

runMarkovChain <- function(n, initial, transition) {
    nstates <- length(initial)
    output <- integer(n)
    output[1] <- sample.int(nstates, 1, prob=initial)
    if (n > 1) {
        for (i in 2:n) {
            prevState <- output[i - 1]
            prob <- transition[prevState, ]
            output[i] <- sample.int(nstates, 1, prob=prob)
        }
    }
    output
}

catvector <- function(x, ...) {
    for (i in 1:length(x)) {
        cat(x[i], ...)
        if (i != length(x)) {
            cat(" ", ...)
        }
    }
    cat("\n", ...)
}

catmatrix <- function(m, ...) {
    for (i in 1:nrow(m)) {
        catvector(m[i, ], ...)
    }
}

generateData <- function(dir, nstates, chain_length, transition_prob, means, sigma) {
    stopifnot(all(apply(transition_prob, 1, sum) == 1))

    rdatafile <- file.path(dir, "data.Rdata")
    datafile <- file.path(dir, "data.txt")

    burn_in_latent <- runMarkovChain(10000, rep(1, nstates), transition_prob)
    initial_prob <- as.vector(table(burn_in_latent)) / length(burn_in_latent)

    hidden_states <- list()
    visible_states <- list()

    for (i in seq_along(chain_length)) {
        l <- chain_length[i]
        hidden_states[[i]] <- runMarkovChain(l, initial_prob, transition_prob)
        visible_states[[i]] <- rnorm(l, mean=means[hidden_states[[i]]], sd=sigma[hidden_states[[i]]])
    }

    save("nstates", "chain_length",
         "transition_prob", "initial_prob",
         "means", "sigma",
         "hidden_states", "visible_states",
         file=rdatafile)


    cat(file=datafile)
    catvector(c(nstates, length(chain_length)), file=datafile, append=TRUE)
    for (i in seq_along(chain_length)) {
        catvector(c(chain_length[i], visible_states[[i]]), file=datafile, append=TRUE)
    }
}

rdirichlet <- function(n, alpha) {
    t(replicate(n, {
        x <- rgamma(length(alpha), 1, alpha)
        x / sum(x)
    }))
}

summarizeResult <- function(result) {
    summarize <- function(value) {
        c(list(count=length(value),
               mean=mean(value),
               median=median(value),
               sd=sd(value),
               iqr=IQR(value),
               mad=mad(value),
               autocor=autocor(value),
               ess=effectiveSampleSize(value)),
          modeHDI(value))
    }
    ret <- result[,
                  summarize(value),
                  by=c('parameter', 'row', 'col')]
    setkeyv(ret, c("parameter", "row", "col"))
    ret["initial_prob", mean := mean / sum(mean)]
    ret["initial_prob", median := median / sum(median)]
    ret["initial_prob", mode := mode / sum(mode)]
    ret["transition_prob", mean := mean / sum(mean), by='row']
    ret["transition_prob", median := median / sum(median), by='row']
    ret["transition_prob", mode := mode / sum(mode), by='row']
}

processResult <- function(dir) {
    load(file.path(dir, "data.Rdata"))

    data <- readResult(file.path(dir, "result.txt"))

    parameterInterval <- function (name) {
        data[name,
             modeHDI(value),
             by=c('row', 'col')]
    }

    no_y_axis <- theme(axis.ticks = element_blank(), axis.text.y = element_blank())

    paramRange <- function(name, positive=FALSE) {
        r <- quantile(data[name, value]$value, c(0.1, 0.9))
        w <- r[2] - r[1]
        pr <- c(r[1] - w/2, r[2] + w/2)
        if (positive) pr[1] <- max(0, pr[1])
        pr
    }

    paramBinWidth <- function(name, positive=FALSE) {
        r <- paramRange(name, positive)
        w <- r[2] - r[1]
        w / 100
    }

    paramAutoCor <- function (param, title) {
        d <- data[param,
                  autocorrelations(value),
                  by=c("seed", "row", "col")]
        (ggplot(d, aes(delta, autocor, fill=as.factor(seed)))
         + ggtitle(title)
         + xlab("Shift") + ylab("Correlation")
         + facet_grid(row ~ col)
         + geom_bar(stat="identity", position="dodge")
         + scale_y_continuous(limits=c(-0.3, 1), breaks=seq(-0.4, 1, 0.2))
         + guides(fill=FALSE)
         )
    }

    bwFit <- HMMFit(visible_states, nStates=nstates)
    bwMeanIdx <- order(bwFit$HMM$distribution$mean)

    stateSeq <- 1:nstates

    theme_set(theme_few())

    pdf(file.path(dir, "stats.pdf"), width=19, height=6)
    (ggplot(data.frame(visible=unlist(visible_states), hidden=unlist(hidden_states)),
            aes(visible, fill=as.factor(hidden)))
    + ggtitle("Overall visible state distribution") + xlab("") + ylab("")
    + geom_histogram()
    + guides(fill=FALSE)
    )
    grid.arrange((ggplot(data["initial_prob"])
                  + ggtitle("Initial probability") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + geom_histogram(aes(value), binwidth=0.01)
                  + scale_x_continuous(limits=c(0, 1))
                  + geom_vline(aes(xintercept=x), data.frame(x=initial_prob, col=stateSeq), color='red')
                  + geom_vline(aes(xintercept=x), data.frame(x=bwFit$HMM$initProb[bwMeanIdx], col=stateSeq), color='green')
                  + geom_vline(aes(xintercept=mode), parameterInterval('initial_prob'), color='black')
                  + no_y_axis
                  ),
                 (ggplot(data["initial_prob"])
                  + ggtitle("Initial probability") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + geom_line(aes(iteration, value, color=as.factor(seed)))
                  + scale_y_continuous(limits=c(0, 1))
                  + geom_hline(aes(yintercept=y), data.frame(y=initial_prob, col=stateSeq), color='red')
                  ),
                 paramAutoCor('initial_prob', "Initial probability"),
                 ncol=3)
    grid.arrange((ggplot(data["transition_prob"], aes(x=value))
                  + ggtitle("Transition probability") + xlab("") + ylab("")
                  + facet_grid(row ~ col)
                  + geom_histogram(binwidth=0.01)
                  + scale_x_continuous(limits=c(0, 1))
                  + geom_vline(aes(xintercept=value), melt(transition_prob, varnames=c("row", "col")), color='red')
                  + geom_vline(aes(xintercept=value), melt(bwFit$HMM$transMat[bwMeanIdx, bwMeanIdx], varnames=c("row", "col")), color='green')
                  + geom_vline(aes(xintercept=mode), parameterInterval('transition_prob'), color='black')
                  + no_y_axis
                  ),
                 (ggplot(data["transition_prob"], aes(iteration, value, color=as.factor(seed)))
                  + ggtitle("Transition probability") + xlab("") + ylab("")
                  + facet_grid(row ~ col)
                  + geom_line()
                  + scale_y_continuous(limits=c(0, 1))
                  + geom_hline(aes(yintercept=value), melt(transition_prob, varnames=c("row", "col")), color='red')
                  ),
                 paramAutoCor('transition_prob', "Transition probability"),
                 ncol=3)
    grid.arrange((ggplot(data["means"])
                  + ggtitle("Means") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + geom_histogram(aes(value), binwidth=paramBinWidth("means"))
                  + scale_x_continuous(limits=paramRange("means"))
                  + geom_vline(aes(xintercept=x), data.frame(x=means, col=stateSeq), color='red')
                  + geom_vline(aes(xintercept=x), data.frame(x=bwFit$HMM$distribution$mean[bwMeanIdx], col=stateSeq), color='green')
                  + geom_vline(aes(xintercept=mode), parameterInterval('means'), color='black')
                  + no_y_axis
                  ),
                 (ggplot(data["means"])
                  + ggtitle("Means") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + geom_line(aes(iteration, value, color=as.factor(seed)))
                  + coord_cartesian(ylim=paramRange("means"))
                  + geom_hline(aes(yintercept=x), data.frame(x=means, col=stateSeq), color='red')
                  ),
                 paramAutoCor('means', "Means"),
                 ncol=3)
    grid.arrange((ggplot(data["sds"])
                  + ggtitle("Standard deviation") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + geom_histogram(aes(value), binwidth=paramBinWidth("sds", positive=TRUE))
                  + scale_x_continuous(limits=paramRange("sds", TRUE))
                  + geom_vline(aes(xintercept=x), data.frame(x=sigma, col=stateSeq), color='red')
                  + geom_vline(aes(xintercept=x), data.frame(x=sqrt(bwFit$HMM$distribution$var[bwMeanIdx]), col=stateSeq), color='green')
                  + geom_vline(aes(xintercept=mode), parameterInterval('sds'), color='black')
                  + no_y_axis
                  ),
                 (ggplot(data["sds"])
                  + ggtitle("Standard deviation") + xlab("") + ylab("")
                  + facet_grid(col ~ .)
                  + coord_cartesian(ylim=paramRange("sds", TRUE))
                  + geom_line(aes(iteration, value, color=as.factor(seed)))
                  + geom_hline(aes(yintercept=y), data.frame(y=sigma, col=stateSeq), color='red')
                  ),
                 paramAutoCor('sds', "Standard deviation"),
                 ncol=3)
    dev.off()
}
