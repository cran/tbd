#' @title Estimate the confidence interval of SACE using bootstrap.
#'
#' @description Give quantiles of bootstrap samples SACE.
#'
#' @param object an object of class \code{sace}.
#' @param nboot a positive integer. The number of bootstrap samples desired.
#' @param seed an integer vector with length \code{nboot}. Seed to generate samples.
#' @param alpha confidence level.
#' @param max.step see documentation of \link[tbd]{sace}.
#' @param singular.ok see documentation of \link[tbd]{sace}.
#' @param print.progress logical. Need progress be printed?
#'
#' @export
#' @return a list with 4 elements:
#' \item{nskip}{number of failures during bootstrap.}
#' \item{sace.boot.record}{a vector with length \code{nboot}-\code{skip}. SACE estimates of all bootstrap samples.}
#' \item{boot.sd}{scaler. Standard deviation of SACE estimates of all bootstrap samples.}
#' \item{ci}{a vector with length 2. Estimated confidence interval.}
#' @author Zhixuan Shao <shaozhixuansh@pku.edu.cn>

boot.ci <- function(object, nboot = 1000, seed = 100:(100 + nboot - 1), alpha = 0.05, max.step = 1000, singular.ok = FALSE, print.progress = TRUE) {
    if (!inherits(object, "sace")) stop("Object must be of class 'sace'")
    if (length(seed)!=nboot) stop("Length of seed should equal nboot!")
    n <- object$n
    data <- object$data
    record <- vector()
    nskip <- 0
    for (i in 1:nboot) {
        if(print.progress) cat("This is the ", i, "th step \n")
        set.seed(seed[i])
        bt.index <- sample(1:n, n, replace = TRUE)
        tmp <- try(sace(Z = data$Z, S = data$S, Y = data$Y, X = data$X, A = data$A, subset = bt.index, max.step = max.step, singular.ok = singular.ok, need.variance = FALSE)$sace)
        if ('try-error' %in% class(tmp)) {
            nskip = nskip + 1
            next
        }
        else (record[length(record) + 1] <- tmp)
        }
    return(list(nskip = nskip, sace.boot.record = record, boot.sd = sd(record), ci = quantile(record, c(alpha / 2, 1 - alpha / 2)))) # quantile bootstrap interval
}