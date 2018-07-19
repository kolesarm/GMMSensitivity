#' One-step estimator based on optimal sensitivity under ell_p constraints
#'
#' Computes the optimal sensitivity and the optimal estimator when the set
#' \eqn{C} takes the form \eqn{B*gamma} where the ell_p norm of M*gamma is
#' bounded by K.
#'
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
#'     \item{k_init}{Initial sensitivity}
#'     \item{g_init}{Moment condition evaluated at initial estimate}
#'
#'     }
#' @param B matrix \eqn{B} with full rank and dimension d_g by d_gamma that
#'     determines the set \eqn{C}, where d_gamma is the number of invalid
#'     moments, and d_g is the number of moments
#' @param M invertible matrix \eqn{M} with dimensions d_gamma by d_gamma that
#'     determines the set \eqn{C}.
#' @param K diameter of set \eqn{C}
#' @param p Parameter determining which ell_p norm to use, one of \code{1},
#'     \code{2}, or \code{Inf}.
#' @param spath Optionally, the solution path, output of \code{lph} to speed up
#'     computation. For \code{p==1} and \code{p==Inf} only.
#' @param alpha determines confidence level, \code{1-alpha}, for
#'     constructing/optimizing confidence intervals.
#' @param opt.criterion Optimality criterion for choosing optimal bias-variance
#'     tradeoff. The options are:
#'
#'    \describe{
#'
#'    \item{\code{"MSE"}}{Mean squared errror}
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
#' @export
OptEstimator <- function(eo, B, M=diag(ncol(B)), K, p=2, spath=NULL, alpha=0.05,
                         opt.criterion="FLCI") {
    if (opt.criterion=="Valid") {
        r <- BuildEstimator(eo$k_init, eo, B, M, K, p, alpha)
        r$opt.criterion <- opt.criterion
        return(r)
    }
    if (p==2)
        return(l2opt(eo, B, M, K, alpha, opt.criterion))

    if (is.null(spath))
        spath <- lph(eo, B, M, p)
    spath <- spath[, -1, drop=FALSE]    # drop lambda/barB
    ep <- BuildEstimator(spath, eo, B, M, K, p, alpha)

    ## Index of criterion to optimize
    idx <- if (opt.criterion=="MSE") {
               which.max(names(ep)=="rmse")
           } else if (opt.criterion=="FLCI") {
               which.max(names(ep)=="hl")
           }
    i <- which.min(ep[[idx]])
    ## Assume minimum is either in [i, i+1] or [i-1, i], which is true if
    ## criterion is convex. Sometimes solution can be stuck, so not enough to
    ## take i+1 (but since which.min takes first minimum, we can always take
    ## i-1)
    ip <- i+1
    while(ip<=nrow(spath) && isTRUE(all.equal(spath[i, ], spath[ip, ])))
        ip <- ip+1

    if (ip<=nrow(spath)) {
        f1 <- function(w)
            BuildEstimator((1-w)*spath[i, ]+w*spath[i+1, ], eo, B, M,
                           K, alpha)[[idx]]
        opt1 <- stats::optimize(f1, interval=c(0, 1))
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            BuildEstimator((1-w)*spath[i-1, ]+w*spath[i, ], eo, B, M,
                           K, alpha)[[idx]]
        opt0 <- stats::optimize(f0, interval=c(0, 1))
    } else {
        opt0 <- list(minimum=1, objective=min(ep[[idx]]))
    }

    kopt <- if (opt1$objective < opt0$objective) {
                (1-opt1$minimum)*spath[i, ] +
                    opt1$minimum*spath[min(i+1, nrow(spath)), ]
            } else {
                (1-opt0$minimum)*spath[max(i-1, 1), ]+opt0$minimum*spath[i, ]
            }

    r <- BuildEstimator(kopt, eo, B, M, K, alpha)
    r$opt.criterion <- opt.criterion
    r
}


#' Build estimator given sensitivity
#' @param k Sensitivity, or matrix of sensitivities
#' @param p 1, 2, or Inf
#' @keywords internal
BuildEstimator <- function(k, eo, B, M, K, p=Inf, alpha=0.05) {
    if (!is.matrix(k)) k <- matrix(k, nrow=1)
    hhat <- drop(eo$h_init + k %*% eo$g_init)
    bf <- function(k, p)
        if (p==Inf) {
            sum(abs(solve(t(M), crossprod(B, k))))
        } else if (p==1) {
            max(abs(solve(t(M), crossprod(B, k))))
        } else {
            sqrt(sum(solve(t(M), crossprod(B, k))^2))
        }
    se <- sqrt(apply(k, 1, function(k)
        drop(crossprod(k, eo$Sig %*% k))) / eo$n)
    bias <- if (ncol(B)==0) 0*se else K * apply(k, 1, bf, p) / sqrt(eo$n)
    hl <- cvb(bias/se, alpha) * se

    structure(list(k=k, h=hhat, bias=bias, se=se, hl=hl, alpha=alpha,
                   rmse=bias^2+se^2),
              class="GMMEstimate")
}


#' Critical value
#' @keywords internal
cvb <- function(b, alpha=0.05)
    sqrt(stats::qchisq(1-alpha, df = 1, ncp = b^2))


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
