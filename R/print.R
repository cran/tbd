#' @title Print results of \code{sace}
#'
#' @description \code{print.sace} prints estimation of the SACE (survivor average causal effect).
#'
#' @param x an object of class \code{sace}.
#' @param ... additional arguments.
#' @method print sace
#' @export
#' @return the input object is returned silently.

print.sace <- function(x, ...) {
    if (!inherits(x, "sace")) stop("Object must be of class 'sace'")

    cat('Call:\n')
    print(x$CALL)
    cat('\n')
    cat("sample size:", x$n, "\n")
    cat("average potential outcomes among control group:", x$mu_0_LL,' ')
    if (x$need.variance) { cat("(s.e. ", sqrt(x$mu_0_LL.var), ")", sep = "") }
    cat('\n')
    cat("average potential outcomes among treatment group:", x$mu_1_LL,' ')
    if (x$need.variance) { cat("(s.e. ", sqrt(x$mu_1_LL.var), ")", sep = "") }
    cat('\n')
    cat("SACE (survivor average causal effect):", x$sace,' ')
    if (x$need.variance) { cat("(s.e. ", sqrt(x$sace.var), ")", sep = "") }
    cat('\n')
}