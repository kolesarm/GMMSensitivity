#' Construct optimal sensitivity and compute worst-case bias and sd for a given
#' lambda
#' @param eo object containig parameter estimates based on theta_init
#' @param B matrix
#' @param MM matrix \eqn{M'M}
#' @param lam weight lambda on bias
#' @param K diameter
#' @param alpha for CI coverage
#' @keywords internal
l2sens <- function(eo, B, MM=diag(ncol(B)), lam, K, alpha=0.05) {
    ## Optimal inverse weight
    BMMB <- B %*% solve(MM, t(B))
    Wi <- eo$Sig + (lam*K^2)*BMMB

    ## Worst-case bias and variance, and sensitivity
    GWG <- crossprod(eo$G, solve(Wi, eo$G))
    k <- drop(-eo$H %*% solve(GWG, t(solve(Wi, eo$G))))
    bias <- K * sqrt(drop(crossprod(k, BMMB %*% k)) / eo$n)
    se <- sqrt(drop(crossprod(k, eo$Sig %*% k)) / eo$n)

    ## Critical value
    cvb <- function(b) {
        if (b>10)
            b + qnorm(1-alpha)
        else
            sqrt(qchisq(1-alpha, df = 1, ncp = b^2))
    }

    list(bias=bias, se=se, hl=cvb(bias/se) * se, k=k)
}


#' One-step estimator based on optimal sensitivity under ell_2 constraints
#'
#' Computes the optimal sensitivity and the optimal estimator when the set
#' \eqn{C} takes the form \eqn{B*gamma} where the ell_2 (Euclidean) norm of
#' M*gamma is bounded by K. Here B is a full rank matrix, and M is a conformable
#' matrix such that M'M is invertible.
#' @param opt.criterion Optimality criterion for choosing optimal bias-variance
#'     tradeoff. The options are:
#'
#'    \describe{
#'
#'    \item{\code{"MSE"}}{maximum MSE}
#'
#'    \item{\code{"FLCI"}}{Length of (fixed-length) two-sided
#'        confidence intervals.}
#'
#'    \item{\code{"Valid"}}{Optimal estimator under valid moments. This returns
#'    the original estimator, with confidence intervals adjusted for possible
#'    misspecification}
#'
#'     }
#'
#' @param alpha determines confidence level, \code{1-alpha}, for
#'     constructing/optimizing confidence intervals.
#' @param eo List containing initial estimates with the following components:
#'
#'     \describe{
#'
#'     \item{Sig}{Estimate of variance of the moment condition, matrix with dimension
#'        d_g by d_g, where d_g is the number of moments}
#'
#'     \item{G}{Estimate of defivative of the moment condition, matrix with
#'     dimension d_g by d_theta, where d_theta is the dimension of \eqn{theta}}
#'
#'     \item{H}{Estimate of defivative of \eqn{h(theta_{0})}. A vector of length d_theta}
#'     \item{n}{sample size}
#'     \item{h_init}{Estimate of \eqn{h(theta_0)}}
#'
#'     \item{g_init}{Moment condition evaluated at initial estimate}
#'
#'     }
#' @param K diameter of set \eqn{C}
#' @param B matrix \eqn{B} with dimensions d_g by d_gamma that determines the
#'     set \eqn{C}, where d_gamma is the number of invalid moments
#' @param MM invertible matrix \eqn{M'M} with dimensions d_gamma by d_gamma that
#'     determines the set \eqn{C}.
#' @export
l2opt <- function(eo, B, MM=diag(ncol(B)), K, alpha=0.05,
                  opt.criterion="FLCI") {

    if (opt.criterion=="MSE")
        lam <- 1
    else if (opt.criterion=="FLCI") {
        crit <- function(lam) l2sens(eo, B, MM, lam, K, alpha)$hl
        ## minimum should be near one
        lam <- optimize(crit, c(0, 5), tol=1e-10)$minimum
    } else {
        lam <- 0
    }
    r <- l2sens(eo, B, MM, lam, K, alpha)

    k_opt <- r$k
    h_opt <- eo$h_init + sum(k_opt * eo$g_init)

    structure(list(k=k_opt, h=h_opt, hl=r$hl, bias=r$bias, se=r$se,
                   opt.criterion=opt.criterion, alpha=alpha, lam=lam),
              class="GMMEstimate")
}

#' @export
print.GMMEstimate <- function(x, digits = getOption("digits"), ...) {

    fmt <- function(x) format(x, digits=digits, width=digits+1)

    r <- as.data.frame(x[c("h", "bias", "se", "hl")])
    r <- cbind(r, l=r$h-r$hl, u=r$h+r$hl)
    r <- fmt(r)
    r$ci <- paste0("(", r$l, ", ",  r$u, ")")
    r$hl <- r$l <- r$u <- NULL
    names(r) <- c("Estimate", "Max. bias", "SE", "CI")
    print(knitr::kable(r))
    cat("\n")
    invisible(x)
}
