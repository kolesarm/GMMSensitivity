#' One-step estimator based on optimal sensitivity under \eqn{\ell_p}{lp}
#' constraints
#'
#' Computes the optimal sensitivity and the optimal estimator when the set
#' \eqn{\mathcal{C}} takes the form \eqn{c=B\gamma}{c=B*gamma} with the
#' \eqn{\ell_p}{lp} norm of \eqn{\gamma}{gamma} bounded by \eqn{M}.
#'
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
#'     \item{H}{Estimate of derivative of \eqn{h(\theta)}{h(theta)}. A vector of
#'     length \eqn{d_\theta}{d_theta}}
#'
#'     \item{n}{sample size}
#'
#'     \item{h_init}{Initial estimate of \eqn{h(\theta)}{h(theta)}}
#'
#'     \item{k_init}{Initial sensitivity}
#'
#'     \item{g_init}{Moment condition evaluated at initial estimate}
#'
#'     }
#' @param B matrix \eqn{B} with full rank and dimension \eqn{d_g} by
#'     \eqn{d_\gamma}{d_gamma} that determines the set \eqn{\mathcal{C}}, where
#'     \eqn{d_\gamma}{d_gamma} is the number of invalid moments, and \eqn{d_g}
#'     is the number of moments
#' @param M Bound on the norm of \eqn{\gamma}{gamma}
#' @param p Parameter determining which \eqn{\ell_p}{lp} norm to use, must equal
#'     \code{1}, \code{2}, or \code{Inf}.
#' @param spath Optionally, the solution path, output of \code{lph} to speed up
#'     computation. For \code{p==1} and \code{p==Inf} only.
#' @param alpha determines confidence level, \code{1-alpha}, for
#'     constructing/optimizing confidence intervals.
#' @param opt.criterion Optimality criterion for choosing optimal bias-variance
#'     tradeoff. The options are:
#'
#'    \describe{
#'
#'    \item{\code{"MSE"}}{Minimize worst-case mean squared error of the
#'    estimator.}
#'
#'    \item{\code{"FLCI"}}{Length of (fixed-length) two-sided confidence
#'        intervals.}
#'
#'    \item{\code{"Valid"}}{Optimal estimator under valid moments. This returns
#'    the original estimator, with confidence intervals adjusted for possible
#'    misspecification}
#'
#'     }
#' @examples
#' ## Replicates estimates in first line of Figure 1 in Armstrong and Kolesár
#' ## (2020)
#' ## 1. Compute matrix B when all instruments are invalid
#' I <- vector(mode="logical", length=nrow(blp$G))
#' I[c(6:13, 20:31)] <- TRUE
#' B <- blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ)[, I, drop=FALSE]
#' ## 2. Collect initial estimates
#' blp$k_init <- -drop(blp$H %*% solve(crossprod(blp$G, blp$W %*% blp$G),
#'                                     crossprod(blp$G, blp$W)))
#' eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
#'            k_init=blp$k_init, h_init= blp$h_init)
#' OptEstimator(eo, B, M=sqrt(sum(I)), p=2, alpha=0.05, opt.criterion="Valid")
#' OptEstimator(eo, B, M=sqrt(sum(I)), p=2, alpha=0.05, opt.criterion="FLCI")
#' @return Object of class \code{"GMMEstimate"}, which is a list with at least
#'     the following components:
#'
#' \describe{
#'
#' \item{h}{Point estimate}
#'
#' \item{bias}{Worst-case bias of estimator}
#'
#' \item{se}{Standard error of estimator}
#'
#' \item{hl}{Half-length of confidence interval, so that the confidence interval
#' takes the form \eqn{h +- hl}}
#'
#' }
#' @references{
#'
#' \cite{Armstrong, T. B., and M. Kolesár (2020): Sensitivity Analysis Using
#' Approximate Moment Condition Models,
#' \url{https://arxiv.org/abs/1808.07387}}
#'
#' }
#' @export
OptEstimator <- function(eo, B, M, p=2, spath=NULL, alpha=0.05,
                         opt.criterion="FLCI") {
    be <- function(k) BuildEstimator(k, eo, B, M, p, alpha)
    if (opt.criterion=="Valid") {
        r <- be(eo$k_init)
        r$opt.criterion <- opt.criterion
        return(r)
    }
    if (p==2)
        return(l2opt(eo, B, M, alpha, opt.criterion))

    if (is.null(spath))
        spath <- lph(eo, B, p)
    ep <- be(spath)

    ## Index of criterion to optimize
    idx <- which.max(names(ep)==switch(opt.criterion, MSE = "rmse",
                                       FLCI = "hl"))
    i <- which.min(ep[[idx]])
    ## Assume minimum is either in [i, i+1] or [i-1, i], which is true if
    ## criterion is convex. Sometimes solution can be stuck, so not enough to
    ## take i+1 (but since which.min takes first minimum, we can always take
    ## i-1)
    ip <- i+1
    while(ip<=nrow(spath) && isTRUE(all.equal(spath[i, ], spath[ip, ])))
        ip <- ip+1

    if (ip <= nrow(spath)) {
        f1 <- function(w)
            be((1-w)*spath[i, ]+w*spath[i+1, ])[[idx]]
        opt1 <- stats::optimize(f1, interval=c(0, 1), tol=tol)
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            be((1-w)*spath[i-1, ]+w*spath[i, ])[[idx]]
        opt0 <- stats::optimize(f0, interval=c(0, 1), tol=tol)
    } else {
        opt0 <- list(minimum=1, objective=min(ep[[idx]]))
    }

    kopt <- if (opt1$objective < opt0$objective) {
                (1-opt1$minimum)*spath[i, ] +
                    opt1$minimum*spath[min(i+1, nrow(spath)), ]
            } else {
                (1-opt0$minimum)*spath[max(i-1, 1), ]+opt0$minimum*spath[i, ]
            }
    r <- be(kopt)
    r$opt.criterion <- opt.criterion

    r
}


## Build estimator given sensitivity
## @param k Sensitivity, or matrix of sensitivities
## @param p 1, 2, or Inf
BuildEstimator <- function(k, eo, B, M, p=Inf, alpha=0.05) {
    if (!is.matrix(k)) k <- matrix(k, nrow=1)
    hhat <- drop(eo$h_init + k %*% eo$g_init)
    se <- sqrt(apply(k, 1, function(k)
        drop(crossprod(k, eo$Sig %*% k))) / eo$n)

    bf <- function(k, p)
        if (p==Inf) {
            sum(abs(crossprod(B, k)))
        } else if (p==1) {
            max(abs(crossprod(B, k)))
        } else {
            sqrt(sum(crossprod(B, k)^2))
        }
    bias <- if (ncol(B)==0) 0*se else M * apply(k, 1, bf, p) / sqrt(eo$n)
    hl <- cvb(bias/se, alpha) * se

    structure(list(k=k, h=hhat, bias=bias, se=se, hl=hl, alpha=alpha,
                   rmse=bias^2+se^2),
              class="GMMEstimate")
}


## Critical value
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
