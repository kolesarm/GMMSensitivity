tol <- .Machine$double.eps^0.75

## Find interval containing zero of a function, then find the zero
##
## Given function \code{f} find \code{x0} such that \code{f(x0)==0}
## @param f function whose root we're looking for
## @param ival upper endpoint of initial interval in which to search
## @param negative logical: should the lower endpoint be \code{1/ival} (if the
##     root is guaranteed to be positive), or \code{-ival}?
FindZero <- function(f, ival=1.1, negative=TRUE) {
    minval <- function(ival) if (negative==TRUE) -ival else min(1/ival, 1e-3)

    while(sign(f(ival))==sign(f(minval(ival))))
            ival <- 2*ival
    stats::uniroot(f, c(minval(ival), ival), tol=tol)$root
}

#' J-test of overidentifying restrictions under local misspecification
#'
#' Computes J-test of overidentifying restrictions with critical value adjusted
#' to allow for local misspecification, when the parameter \eqn{c} takes the
#' form \eqn{c=B\gamma}{c=B*gamma} with the \eqn{\ell_p}{lp} norm of
#' \eqn{\gamma}{gamma} bounded by \eqn{M}.
#'
#' The test assumes initial estimator in \code{eo} is optimal under correct
#' specification, computed using \code{eo$Sig} as the weight matrix. The test is
#' based on a J statistic using critical values that account for local
#' misspecification; see appendix B in Armstrong and Kolesár (2020) for details.
#' @param eo List containing initial estimates with the following components:
#'
#'     \describe{
#'
#'     \item{Sig}{Estimate of variance of the moment condition, matrix with
#'        dimension \eqn{d_g} by \eqn{d_g}, where \eqn{d_g} is the number of
#'        moments}
#'
#'     \item{G}{Estimate of derivative of the moment condition, matrix with
#'     dimension \eqn{d_g} by \eqn{d_\theta}{d_theta}, where
#'     \eqn{d_\theta}{d_theta} is the dimension of \eqn{\theta}{theta}}
#'
#'     \item{n}{sample size}
#'
#'     \item{g_init}{Moment condition evaluated at initial estimate}
#'
#'     }
#' @inheritParams OptEstimator
#' @examples
#' ## Replicates first line of Table 1 in Armstrong and Kolesár (2020)
#' ## 1. Compute matrix B when instrument D/F # cars is invalid
#' I <- vector(mode="logical", length=nrow(blp$G))
#' I[6] <- TRUE
#' B <- (blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ))[, I, drop=FALSE]
#' ## 2. Make sure Sig corresponds to inverse of weight matrix
#' eo <- list(G=blp$G, Sig=solve(blp$W), n=blp$n, g_init=blp$g_init)
#' Jtest(eo, B, M=1, p=2, alpha=0.05)
#' Jtest(eo, B, M=1, p=Inf, alpha=0.05)
#' @return List with three elements: \describe{
#'
#' \item{J}{Value of J statistic}
#'
#' \item{p0}{P-value of usual J test}
#'
#' \item{pC}{P-value for J-test that  allows for local misspecification}
#'
#' \item{Mmin}{Minimum value of \eqn{M} for which the J-test does not reject} }
#' @references{
#'
#' \cite{Armstrong, T. B., and M. Kolesár (2020): Sensitivity Analysis Using
#' Approximate Moment Condition Models,
#' \url{https://arxiv.org/abs/1808.07387}}
#'
#' }
#' @export
Jtest <- function(eo, B, M=1, p=2, alpha=0.05) {
    J <- eo$n*drop(crossprod(eo$g_init, solve(eo$Sig, eo$g_init)))
    ## Sigma^{-1/2}
    e <- eigen(eo$Sig)
    Sig12 <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)
    SG <- Sig12 %*% eo$G
    R <- diag(nrow(eo$G))-SG %*% solve(crossprod(SG), t(SG))
    A <- R %*% Sig12 %*% B

    kbar <- switch(as.character(p),
                   "2"=max(eigen(crossprod(A))$values),
                   "1"=max(colSums(A^2)),
                   "Inf"=-stats::optim(par=rep(1, ncol(A)), function(x)
                    -drop(crossprod(A%*%x))/max(x^2), method="BFGS")$value)

    ## p-value under correct specification
    p0 <- 1-stats::pchisq(q=J, df=nrow(eo$G)-ncol(eo$G), ncp=0)

    ## Value of M for which p-value is 0.05
    M0 <- 0
    if (p0<=alpha) {
        M0 <- FindZero(function(M)
            0.95-stats::pchisq(q=J, df=nrow(eo$G)-ncol(eo$G),
                               ncp=M^2*kbar),
            ival=2, negative=FALSE)
    }

    list(J=J, p0=p0,
         pC=1-stats::pchisq(q=J, df=nrow(eo$G)-ncol(eo$G), ncp=M^2*kbar),
         Mmin=M0)
}
