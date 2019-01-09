#' @title Summarize results of \code{sace}
#'
#' @description \code{summary.sace} summary estimation of the SACE (survivor average causal effect) and all other model parameters.
#'
#' @note 
#' If \code{need.variance} is \code{TRUE}, \code{sace} must have been called with \code{need.variance == TRUE}, so that the information needed was recorded.
#'
#' @param object an object of class \code{sace}.
#' @param ... additional arguments.
#' @method summary sace
#' @export
#' @return the input object is returned silently.


summary.sace <- function(object, ...) {
    if (!inherits(object, "sace")) stop("Object must be of class 'sace'")
    cat('Call:\n')
    print(object$CALL)
    cat('\n')
    cat("sample size:", object$n, "\n")
    cat('\n')
    cat("----- Parameter Estimation -----", "\n")
    cat('\n')
    if (object$need.variance) {
        coef.summary <- function(coef, coef.var) {
            coef.table <- cbind(coef, sqrt(diag(coef.var)))
            coef.table <- cbind(coef.table, coef.table[, 1] / coef.table[, 2])
            coef.table <- cbind(coef.table, 2 * (1 - pnorm(abs(coef.table[, 3]))))
            coef.table[, c(3, 4)] <- round(coef.table[, c(3, 4)], digits = 3)
            colnames(coef.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            return(coef.table)
        }
        
        cat('beta:', '\n')
        print(coef.summary(object$beta, object$beta.var))
        cat('\n')
        cat('gamma:', '\n')
        print(coef.summary(object$gamma, object$gamma.var))
        cat('\n')
        cat('alpha_1:', '\n')
        print(coef.summary(object$alpha_1, object$alpha_1.var))
        cat('\n')
        cat('alpha_2:', '\n')
        print(coef.summary(object$alpha_2, object$alpha_2.var))
    }
    else {
        print(object$beta)
        print(object$gamma)
        print(object$alpha_1)
        print(object$alpha_2)
    }
    cat('\n')
    cat("----- Average Potential Outcomes and SACE among Always-Survivor Group -----", '\n')
    cat('\n')

    if (object$need.variance) {
        cat("average potential outcomes among control group: mu_0_LL", '\n')
        print(coef.summary(object$mu_0_LL, object$mu_0_LL.var))
        cat('\n')
        cat("average potential outcomes among treatment group: mu_1_LL", '\n')
        print(coef.summary(object$mu_1_LL, object$mu_1_LL.var))
        cat('\n')
        cat("survivor average causal effect: SACE", '\n')
        print(coef.summary(object$sace, object$sace.var))
    }
    else {
        cat("average potential outcomes among control group:", object$mu_0_LL, '\n')
        cat("average potential outcomes among treatment group:", object$mu_1_LL, '\n')
        cat("SACE (survivor average causal effect):", object$sace, '\n')
    }
}

#summary(t)