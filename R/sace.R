#' @title Estimation of causal effects with outcomes truncated by death
#'
#' @description \code{sace} estimates survivor average causal effects (SACE) with outcomes truncated by death.
#'
#' @details
#' This function \code{sace}, gives estimation of average causal effects (ACE) with outcomes truncated by death. The identification of SACE relies on the existence of a substitution variable and requires the assumptions of monotonicity, ignorability, exclusion restriction, and relevance. While the naive estimates given by the coefficient of \code{Z} from \code{lm(Y ~ Z + X + A, subset = S == 1)} are restricted among survivors and therefore may be subject to selection bias, this method gives consistent estimates of the SACE (survivor average causal effect), defined as the average causal effect among the subgroup consisting of subjects who would survive under either exposure, i.e. among the always-survivor group (\eqn{G=LL}). See references for details of the assumptions and the model parameterizations.
#'
#' Parameters \code{beta} and \code{gamma} are estimated by MLE, using \link[stats]{optim}.
#'
#' If \code{need.variance == TRUE}, the asymptotic variance estimators of both parameters and estimators will be given. This requires the \pkg{numDeriv} package.
#'
#' @note
#' The length of vectors \code{Z}, \code{Y}, \code{S}, as well as the row number of matrix \code{X} and \code{A} must equal the sample size \code{n}.
#'
#' @param Z a logical vector. Exposure indicator. Convetionally, \code{1} means treatment and \code{0} means control. Must not have missing values.
#' @param S a logical vector. Survival indicator. \code{1} means survival and \code{0} means death. Must not have missing values.
#' @param Y a numeric vector. (Univariate) outcomes. May have \code{NA} where \eqn{S=0} (since \eqn{Y} is not well-defined where \eqn{S=0}).
#' @param X an optional numeric matrix or vector. Baseline covariates. 
#' @param A an optional numeric matrix or vector. Substitution variable(s) which satisfies the assumptions of "exclusion restriction" and "substitution relevance". See references. If \code{A == NULL}, then the naive method, namely OLS, will be used.
#' @param subset an optional vector specifying a subset of obervations to be used.
#' @param optim.method The method to be used for maximum likelihood optimization. See \link[stats]{optim}.
#' @param max.step integer. Maximum iterating steps of maximum likelihood optimization.
#' @param singular.ok logical. Refers to the OLS estimation of the coefficients \code{alpha_1} and \code{alpha_2} using \link[stats]{lm}. If \code{FALSE} (default), a singular fit raises an error.
#' @param need.variance logical. Is variance of parameters and estimators needed? See details.
#' @param hessian logical. If \code{TRUE}, the hessian returned by \link[stats]{optim} will be used to compute the information matrix. If \code{FALSE}, the matrix will be calculated by an explicit formula.
#' @export
#' @return a list with following elements:
#' \item{CALL}{function call.}
#' \item{data}{data used (within the specified subset).}
#' \item{optim.method}{method used for optimization.}
#' \item{need.variance}{is variance of parameters and estimators needed?}
#' \item{n}{sample size.}
#' \item{mu_0_LL}{average potential outcomes among control group, \eqn{E[ Y(0) | G=LL ]}.}
#' \item{mu_1_LL}{average potential outcomes among treatment group, \eqn{E[ Y(1) | G=LL ]}.}
#' \item{sace}{survivor average causal effect, equals \code{mu_1_LL}-\code{mu_0_LL}.}
#' \item{beta}{\eqn{Pr{S(1)=1| X,A}=expit(\beta_0+X' \beta_1+ A \beta_2)}, estimated by MLE.}
#' \item{gamma}{\eqn{Pr{S(0)=1| X,A}/Pr{S(1)=1| X,A}=expit(\gamma_0+X' \gamma_1+ A \gamma_2)}, estimated by MLE.}
#' \item{beta_gamma.convergence}{indicator of convergence of MLE optimization of beta and gamma. 0 means convergence. See \link[stats]{optim}.}
#' \item{alpha_1}{\eqn{E[Y(0)| Z=0, G=LL, X, A ]=\alpha_{10}+X' \alpha_{11}+ A \alpha_{12}}, coefficients of \code{lm(Y ~ 1 + X + A, subset = Z == 0)}.}
#' \item{alpha_2}{\eqn{E[Y(1)| Z=1, G=LL, X, A ]=\alpha_{20}+X' \alpha_{21}+ G \alpha_{22}}, coefficients of \code{lm(Y ~ 1 + X + W.expit, subset = (Z == 1 & S == 1))}.}
#' The following items will be given only if \code{need.variance == TRUE}:
#' \item{beta.var}{estimated asymptotic covariance matrix of beta.}
#' \item{gamma.var}{estimated asymptotic covariance matrix of gamma.}
#' \item{relevance.Pvalue}{P value of the asymptotic chi-squared test on the relevance assumption for the substitution variable. A large P value suggests that the relevance assumption may not hold, namely, the substitution variable(s) may have little impact on the latent survival type.}
#' \item{alpha_1.var}{estimated asymptotic covariance matrix of alpha_1.}
#' \item{alpha_2.var}{estimated asymptotic covariance matrix of alpha_2.}
#' \item{mu_0_LL.var}{estimated asymptotic variance of mu_0_LL.}
#' \item{mu_1_LL.var}{estimated asymptotic variance of mu_1_LL.}
#' \item{sace.var}{estimated asymptotic variance of the SACE.}
#' @author Linbo Wang <linbo.wang@utoronto.ca>
#' @author Zhixuan Shao <shaozhixuansh@pku.edu.cn>
#' @references Linbo Wang, Xiao-Hua Zhou, Thomas S. Richardson; Identification and estimation of causal effects with outcomes truncated by death, Biometrika, Volume 104, Issue 3, 1 September 2017, Pages 597-612, \url{https://doi.org/10.1093/biomet/asx034}
#' @examples
#' attach(simulated_data)
#' X <- cbind(X.X1, X.V2, X.V3)
#' sace.result <- sace(Z, S, Y, X, A)
#' sace


sace <- function(Z, S, Y, X, A, subset, optim.method = "BFGS", max.step = 1000, singular.ok = TRUE, need.variance = TRUE, hessian=TRUE) {
    ### 0. Checks arguments ###################################

    ## 0.1 Checks data type ###################################

    if (!is.vector(Z)) stop("Z should be a vector.")
    Z <- as.logical(Z)
    if (!is.vector(S)) stop("S should be a vector.")
    S <- as.logical(S)
    if (!is.vector(Y)) stop("Y should be a vector.")
    if (missing(X)) X <- NULL
    else {
        X <- as.matrix(X)
        if (!is.matrix(X)) stop("X must be either a vector or a matrix.")
    }
    if (missing(A)) A <- NULL
    else {
        A <- as.matrix(A)
        if (!is.matrix(A)) stop("A must be either a vector or a matrix.")
        }

    ## 0.2 Length check ###################################

    if (missing(subset)) { # check if Z,S,Y,A,X have same length
        n <- length(Z)
        if (length(S) != n) stop("S must have the same length as Z's.")
        if (length(Y) != n) stop("Y must have the same length as Z's.")
        if ((!is.null(X)) & (nrow(X) != n)) stop("X must have the same row number as Z's length.")
        if ((!is.null(A)) & (nrow(A) != n)) stop("A must have the same row number as Z's length.")
    }

    ## 0.3 Subset ##################################

    else { # subset is provided
        Z <- Z[subset]
        S <- S[subset]
        Y <- Y[subset]
        if (!is.null(A)) A <- A[subset,]
        if (!is.null(X)) X <- X[subset,]
        n <- length(Z)
    }

    ## 0.4 Missing values ###############################

    if (sum(is.na(Z))) stop("Z should not have missing values.")
    if (sum(is.na(S))) stop("S should not have missing values.")
    if (sum(is.na(Y[S]))) stop("Y should not have missing values where S==1.")
    if ((!is.null(X)) & (sum(is.na(X)))) stop("X should not have missing values.")
    if (is.null(A)) {
        warning("A is not provided, naive method (OLS) is used.")
        if (is.null(X)) return(lm(Y ~ Z, subset = S))
        else return(lm(Y ~ Z + X, subset = S))
        }
    else { if (sum(is.na(A))) stop("A should not have missing values.") }


    ### 1. Sets up data ###########################
    s1 = S == 1;
    s0_z1 = S == 0 & Z == 1
    s1_z0 = S == 1 & Z == 0;
    s0_z0 = S == 0 & Z == 0
    sz = cbind(s1, s0_z1, s1_z0, s0_z0)

    W <- cbind(rep(1, n), X, A)
    d <- ncol(W)

    ### 2. Defining interior functions ###########################################

    expit <- function(x) exp(x) / (1 + exp(x))

    nLL_beta <- function(beta, gamma, W, sz) {
        ebeta = expit(W %*% beta)
        egamma = expit(W %*% gamma)
        loglike = sum(sz[, 1] * log(ebeta) + sz[, 2] * log(1 - ebeta) + sz[, 3] * log(egamma) + sz[, 4] * log(1 - ebeta * egamma))
        return(loglike)
    }

    beta_gamma.gr <- function(beta, gamma, W, sz) {
        ebeta = as.vector(expit(W %*% beta))
        egamma = as.vector(expit(W %*% gamma))
        loglike.partial.beta <- colSums((sz[, 1] * (1 - ebeta) + sz[, 2] * (-ebeta) + sz[, 4] / (1 - ebeta * egamma) * (-egamma) * (ebeta) * (1 - ebeta)) * W)
        loglike.partial.gamma <- colSums((sz[, 3] * (1 - egamma) + sz[, 4] / (1 - ebeta * egamma) * (-ebeta) * (egamma) * (1 - egamma)) * W)
        return (c(loglike.partial.beta,loglike.partial.gamma))
    }

    beta_gamma_joint.hessian <- function(beta, gamma, W, sz) {
        Wbeta <- as.vector(W %*% beta)
        ebeta <- expit(Wbeta)
        Wgamma <- as.vector(W %*% gamma)
        egamma <- expit(Wgamma)
        d <- length(beta)
        loglike.hessian.beta <- matrix(rep(0, d ^ 2), nrow = d)
        loglike.hessian.gamma <- matrix(rep(0, d ^ 2), nrow = d)
        loglike.hessian.beta_gamma <- matrix(rep(0, d ^ 2), nrow = d)
        for (i in 1:n) {
            loglike.hessian.beta <- loglike.hessian.beta - (sz[i, 1] + sz[i, 2] + sz[i, 4] * egamma[i] / (1 - ebeta[i] * egamma[i]) * (-ebeta[i] + (1 - ebeta[i]) / (1 - ebeta[i] * egamma[i]))) * ebeta[i] * (1 - ebeta[i]) * outer(W[i,], W[i,])
            loglike.hessian.gamma <- loglike.hessian.gamma - (sz[i, 3] + sz[i, 4] * ebeta[i] / (1 - ebeta[i] * egamma[i]) * (-egamma[i] + (1 - egamma[i]) / (1 - ebeta[i] * egamma[i]))) * egamma[i] * (1 - egamma[i]) * outer(W[i,], W[i,])
            loglike.hessian.beta_gamma <- loglike.hessian.beta_gamma - sz[i, 4] * ebeta[i] * (1 - ebeta[i]) * egamma[i] * (1 - egamma[i]) / (1 - ebeta[i] * egamma[i]) ^ 2 * outer(W[i,], W[i,])
        }
        loglike.hessian.beta_gamma_joint <- matrix(rep(0, (2 * d) ^ 2), nrow = 2 * d)
        loglike.hessian.beta_gamma_joint[1:d, 1:d] <- loglike.hessian.beta
        loglike.hessian.beta_gamma_joint[-(1:d), - (1:d)] <- loglike.hessian.gamma
        loglike.hessian.beta_gamma_joint[1:d, - (1:d)] <- t(loglike.hessian.beta_gamma)
        loglike.hessian.beta_gamma_joint[-(1:d), 1:d] <- loglike.hessian.beta_gamma
        return(loglike.hessian.beta_gamma_joint)
    }

    ### 3. Parameter estimation ##############################################

    ####-------- 3.1 Estimating $\beta$ and $\gamma$--------
    #lm.s1_z1 <- glm(S~X+A,subset=Z==1,family=binomial)

    #beta <- rep(0, d) #initial values
    #gamma <- rep(0, d) #initial values
    #Diff <- function(x, y) sum((x - y) ^ 2) / sum(x ^ 2 + thres) # normalized euclidean distance
    #diff <- thres + 1
    #step <- 0
    #opt1 <- opt2 <- NULL
    #while (diff > thres & step < max.step) {
        #opt1 <- optim(beta, nLL_beta, gamma = gamma, W = W, sz = sz,
                      #method = optim.method,
                      #gr = beta.gr,
                      #hessian = TRUE,
                      #control = list(fnscale = -1, maxit = max.step))
        #diff1 <- Diff(opt1$par, beta)
        #beta <- opt1$par
        #cat("loglik1:", nLL_beta(beta, gamma, W, sz), '\n')
        #cat("loglik.partial.beta:", beta.gr(beta, gamma, W, sz), '\n')
        #require(numDeriv)
        #cat("loglik.partial.beta.numer:", grad(nLL_beta, beta, gamma = gamma, W = W, sz = sz),'\n')
        #opt2 <- optim(gamma, nLL_gamma, beta = beta, W = W, sz = sz,
                      #method = optim.method,
                      #gr = gamma.gr,
                      #hessian = TRUE,
                      #control = list(fnscale = -1, maxit = max.step))
        #diff <- max(diff1, Diff(opt2$par, gamma))
        #gamma <- opt2$par
        #cat("loglik2:", nLL_beta(beta, gamma, W, sz), '\n')
        #cat("loglik.partial.gamma:", gamma.gr(gamma, beta, W, sz), '\n')
        #cat("loglik.partial.gamma.numer:", grad(nLL_gamma, gamma,beta= beta, W = W, sz = sz,),'\n')
        #step <- step + 1
    #}
    opt3 <- optim(c(rep(0, d), rep(0, d)),
                  function(beta_gamma, W, sz) nLL_beta(beta_gamma[1:d], beta_gamma[-(1:d)], W, sz),
                  W = W, sz = sz,
                  gr = function(beta_gamma, W, sz) beta_gamma.gr(beta_gamma[1:d], beta_gamma[-(1:d)], W, sz),
                  method = optim.method,
                  hessian = need.variance&hessian,
                  control = list(fnscale = -1, maxit = 2 * max.step))
    if (opt3$convergence != 0) { warning(paste("Optimization of beta and gamma didn't converge in", max.step, "steps !")) }
    beta <- opt3$par[1:d]
    gamma <- opt3$par[-(1:d)]
    #cat("opt3$beta:", beta, '\n')
    #cat("opt3$gamma:", gamma, '\n')
    #cat("loglik_3:", nLL_beta(beta,gamma , W, sz), '\n')
    #cat("loglik.partial.beta:", beta.gr(beta, gamma, W, sz), '\n')
    #cat("loglik.partial.beta.numer:", grad(nLL_beta, beta, gamma = gamma, W = W, sz = sz), '\n')
    #cat("loklike_1:", opt1$value, "\n")
    #cat("loglike_2:", opt2$value, "\n")
    #cat("hessian_beta_formula:", '\n')
    #print(beta.hessian(beta, gamma, W, sz))
    #cat("hessian_beta_numer:", '\n')
    ##print(opt1$hessian)
    #print(hessian(nLL_beta, beta, gamma = gamma, W = W, sz = sz))
    #cat("hessian_gamma_formula:", '\n')
    #print(gamma.hessian(gamma, beta, W, sz))
    #cat("hessian_beta_gamma_formula:", '\n')
    #print(beta_gamma.hessian(beta, gamma, W, sz))
    #cat("hessian_beta_gamma_numer:", '\n')
    #print(opt3$hessian)

    ## 3.2 Estimating alpha_1 (for Z==0) ################

    if (is.null(X)) lm.y.z0 <- lm(Y ~ 1 + A, subset = Z == 0, singular.ok = singular.ok)
    else lm.y.z0 <- lm(Y ~ 1 + X + A, subset = Z == 0, singular.ok = singular.ok)
                     
    alpha_1 <- lm.y.z0$coef

    ## 3.3 Estimating alpha_2 (for Z==1) ##################


    #W.expit <- expit(W %*% gamma) # W.expit is the estimate of G
    #if (is.null(X)) lm.y.z1 <- lm(Y ~ 1 + W.expit, subset = (Z == 1 & S == 1), singular.ok = singular.ok)
    #else lm.y.z1 <- lm(Y ~ 1 + X + W.expit, subset = (Z == 1 & S == 1), singular.ok = singular.ok)
    

    alpha_2.gamma <- function(gamma) {
        W.expit <- expit(W %*% gamma) # W.expit is the estimate of G
        if (is.null(X)) lm.y.z1 <- lm(Y ~ 1 + W.expit, subset = (Z == 1 & S == 1), singular.ok = singular.ok)
        else lm.y.z1 <- lm(Y ~ 1 + X + W.expit, subset = (Z == 1 & S == 1), singular.ok = singular.ok)
        return(lm.y.z1)
    }
    lm.y.z1<-alpha_2.gamma(gamma)

    alpha_2 <- lm.y.z1$coef
    #if (length(alpha_2) < d) alpha_2[d] <- 0


    ### 4. SACE estimation ###################################################

    ## 4.1 E[Y(0)|G=LL] \equiv sace_z0

    #if (is.null(X)) mu_0_LL_W <- predict(lm.y.z0, data.frame(A = A))
    #else mu_0_LL_W <- predict(lm.y.z0, data.frame(X = X, A = A)))
    mu_0_LL.fun <- function(alpha_1.beta.gamma, W) {
        alpha_1 <- alpha_1.beta.gamma[1:d]
        beta <- alpha_1.beta.gamma[(d + 1):(2 * d)]
        gamma <- alpha_1.beta.gamma[(2 * d + 1):(3 * d)]
        mu_0_LL_W <- W %*% alpha_1 # n*1
        LL_W <- expit(W %*% beta) * expit(W %*% gamma) #n*1
        mu_0_LL <- sum(mu_0_LL_W * LL_W) / sum(LL_W) #scaler
        return(mu_0_LL)
    }

    mu_0_LL <- mu_0_LL.fun(c(alpha_1, beta, gamma), W)

    ## 4.2 E[Y(1)|G=LL] \equiv sace_z1

    mu_1_LL.fun <- function(alpha_2.beta.gamma, W) {
        alpha_2 <- alpha_2.beta.gamma[1:length(alpha_2)]
        beta <- alpha_2.beta.gamma[(length(alpha_2) + 1):(length(alpha_2) + d)]
        gamma <- alpha_2.beta.gamma[(length(alpha_2) + d + 1):(length(alpha_2) + 2 * d)]
        Coef <- cbind(rep(1, n), X, rep(1, n))
        mu_1_LL_W <- Coef %*% alpha_2 #n*1
        LL_W <- expit(W %*% beta) * expit(W %*% gamma) #n*1
        mu_1_LL <- sum(mu_1_LL_W * LL_W) / sum(LL_W) #scaler
        return(mu_1_LL)
    }
    
    
    mu_1_LL <- mu_1_LL.fun(c(alpha_2, beta, gamma), W)

    sace <- mu_1_LL - mu_0_LL
    
    names(beta) <- names(alpha_1)
    names(gamma) <- names(alpha_1)

    results <- list(CALL = match.call(),
                    data = list(Z = Z, S = S, Y = Y, X = X, A = A), n = n,
                    optim.method = optim.method,
                    need.variance=need.variance,
                    mu_0_LL = mu_0_LL, mu_1_LL = mu_1_LL,
                    sace = sace,
                    beta = beta, gamma = gamma,
                    beta_gamma.convergence = opt3$convergence,
                    alpha_1 = alpha_1, alpha_2 = alpha_2)

    #########-----------estimate the asymptotic variance-----------####
    if (need.variance) {
        #require(numDeriv)
        failure <- FALSE
        if (hessian == TRUE) {
            beta_gamma.var <- try(solve(-opt3$hessian))
            if ('try-error' %in% class(beta_gamma.var)) {
                warning("The retured Hessian is not reversible. Use alternative method.")
                failure <- TRUE
            }
            else {
                pos_definite <- try(prod(diag(beta_gamma.var) >= 0))
                if ('try-error' %in% class(pos_definite) | is.na(pos_definite) | (!pos_definite)) {
                    warning("The returned Hessian is not positive definite! Use alternative method.")
                    failure <- TRUE
                }
            }
        }
        if (hessian == FALSE | failure == TRUE) {
            #ebeta = as.vector(expit(W %*% beta))
            #egamma = as.vector(expit(W %*% gamma))
            #score_at_estimated <- cbind((sz[, 1] * (1 - ebeta) + sz[, 2] * (-ebeta) + sz[, 4] / (1 - ebeta * egamma) * (-egamma) * (ebeta) * (1 - ebeta)) * W, (sz[, 3] * (1 - egamma) + sz[, 4] / (1 - ebeta * egamma) * (-ebeta) * (egamma) * (1 - egamma)) * W)
            #sum_tmp <- 0
            #for (i in ncol(score_at_estimated)) {
                #sum_tmp <- sum_tmp + outer(score_at_estimated[i,], score_at_estimated[i,])
            #}
            #print(sum_tmp)
            #beta_gamma.var <- try(solve(sum_tmp))
            beta_gamma.var <- try((-1) * solve(beta_gamma_joint.hessian(beta, gamma, W, sz)))
            if ('try-error' %in% class(beta_gamma.var)) {
                warning("Failed to estimate variance: Information matrix not reversible!")
                results <- c(results, list(sace.var = Inf))
                class(results) <- c("sace", "list")
                return(results)
            }
        }
        beta.var <- beta_gamma.var[1:d, 1:d]
        gamma.var <- beta_gamma.var[-(1:d), - (1:d)]
        #if ((P_value <- 2 * (1 - max(pnorm(abs(beta[d] / sqrt(beta.var[d, d]))), pnorm(abs(gamma[d] / sqrt(gamma.var[d, d])))))) > 0.10) {
        P_value <- pchisq(c(tail(beta, ncol(A)), tail(gamma, ncol(A))) %*% solve(beta_gamma.var[c(tail(1:d, ncol(A)), tail((d + 1):(2 * d), ncol(A))), c(tail(1:d, ncol(A)), tail((d + 1):(2 * d), ncol(A)))]) %*% c(tail(beta, ncol(A)), tail(gamma, ncol(A))), 2 * ncol(A), lower.tail = FALSE)
        if (P_value > 0.10) {
        #if ((P_value <- 2 * (1 - pnorm(abs(gamma[d] / sqrt(gamma.var[d, d]))))) > 0.20) {
           
            warning(paste("Substitution variable (A) had insignificant effect on survival (S)!\n P Value =", P_value))
        }

        alpha_1.var <- vcov(lm.y.z0)

        alpha_2.var <- vcov(lm.y.z1)
        # additonal variance of alpha_2 due to the randomness of W.expit
        alpha_2.partial.gamma <- matrix(nrow = length(alpha_2), ncol = length(gamma))
        for (i in 1:length(alpha_2)) {
            alpha_2.partial.gamma[i,] <- grad(function(x)(alpha_2.gamma(x)$coef)[i], gamma)
        }
        alpha_2.var <- alpha_2.var + alpha_2.partial.gamma %*% gamma.var %*% t(alpha_2.partial.gamma)

        mu_0_LL.grad <- try(grad(mu_0_LL.fun, c(alpha_1, beta, gamma), W = W))
        if ('try-error' %in% class(mu_0_LL.grad)) {
            warning("Failed to estimate variance!")
            results <- c(results, list(sace.var = Inf))
            class(results) <- c("sace", "list")
            return(results)
        }

        mu_1_LL.grad <- try(grad(mu_1_LL.fun, c(alpha_2, beta, gamma), W = W))
        if ('try-error' %in% class(mu_1_LL.grad)) {
            warning("Failed to estimate variance!")
            results <- c(results, list(sace.var = Inf))
            class(results) <- c("sace", "list")
            return(results)
        }

        alpha_1_alpha_2_beta_gamma.var <- matrix(rep(0, (3 * d + length(alpha_2)) ^ 2), nrow = 3 * d + length(alpha_2))
        alpha_1_alpha_2_beta_gamma.var[1:d, 1:d] <- alpha_1.var
        alpha_1_alpha_2_beta_gamma.var[(d + 1):(d + length(alpha_2)), (d + 1):(d + length(alpha_2))] <- alpha_2.var
        alpha_1_alpha_2_beta_gamma.var[(d + length(alpha_2) + 1):(3 * d + length(alpha_2)), (d + length(alpha_2) + 1):(3 * d + length(alpha_2))] <- beta_gamma.var

        mu_0_LL.var <- mu_0_LL.grad %*% alpha_1_alpha_2_beta_gamma.var[-((d + 1):(d + length(alpha_2))), - ((d + 1):(d + length(alpha_2)))] %*% mu_0_LL.grad

        mu_1_LL.var <- mu_1_LL.grad %*% alpha_1_alpha_2_beta_gamma.var[-(1:d), - (1:d)] %*% mu_1_LL.grad

        sace.grad <- rep(0, 3 * d + length(alpha_2))
        sace.grad[1:d] <- (-1) * mu_0_LL.grad[1:d]
        sace.grad[(d + 1):(d + length(alpha_2))] <- mu_1_LL.grad[1:length(alpha_2)]
        sace.grad[(d + length(alpha_2) + 1):(3 * d + length(alpha_2))] <- mu_1_LL.grad[-(1:length(alpha_2))] - mu_0_LL.grad[-(1:d)]

        sace.var <- sace.grad %*% alpha_1_alpha_2_beta_gamma.var %*% sace.grad
        
        results <- c(results, list(beta.var = beta.var, gamma.var = gamma.var, relevance.Pvalue = P_value, alpha_1.var = alpha_1.var, alpha_2.var = alpha_2.var, mu_0_LL.var = mu_0_LL.var, mu_1_LL.var = mu_1_LL.var, sace.var = sace.var))
        
    }
    class(results) <- c("sace", "list")
    return(results)
}