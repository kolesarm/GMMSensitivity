% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opt.R
\name{OptEstimator}
\alias{OptEstimator}
\title{One-step estimator based on optimal sensitivity under \eqn{\ell_p}{lp}
constraints}
\usage{
OptEstimator(
  eo,
  B,
  M,
  p = 2,
  spath = NULL,
  alpha = 0.05,
  opt.criterion = "FLCI"
)
}
\arguments{
\item{eo}{List containing initial estimates with the following components:

    \describe{

    \item{Sig}{Estimate of variance of the moment condition, matrix with
       dimension \eqn{d_g} by \eqn{d_g}, where \eqn{d_g} is the number of
       moments}

    \item{G}{Estimate of derivative of the moment condition, matrix with
    dimension \eqn{d_g} by \eqn{d_\theta}{d_theta}, where
    \eqn{d_\theta}{d_theta} is the dimension of \eqn{\theta}{theta}}

    \item{H}{Estimate of derivative of \eqn{h(\theta)}{h(theta)}. A vector of
    length \eqn{d_\theta}{d_theta}}

    \item{n}{sample size}

    \item{h_init}{Initial estimate of \eqn{h(\theta)}{h(theta)}}

    \item{k_init}{Initial sensitivity}

    \item{g_init}{Moment condition evaluated at initial estimate}

    }}

\item{B}{matrix \eqn{B} with full rank and dimension \eqn{d_g} by
\eqn{d_\gamma}{d_gamma} that determines the set \eqn{\mathcal{C}}, where
\eqn{d_\gamma}{d_gamma} is the number of invalid moments, and \eqn{d_g}
is the number of moments}

\item{M}{Bound on the norm of \eqn{\gamma}{gamma}}

\item{p}{Parameter determining which \eqn{\ell_p}{lp} norm to use, must equal
\code{1}, \code{2}, or \code{Inf}.}

\item{spath}{Optionally, the solution path, output of \code{lph} to speed up
computation. For \code{p==1} and \code{p==Inf} only.}

\item{alpha}{determines confidence level, \code{1-alpha}, for
constructing/optimizing confidence intervals.}

\item{opt.criterion}{Optimality criterion for choosing optimal bias-variance
    tradeoff. The options are:

   \describe{

   \item{\code{"MSE"}}{Minimize worst-case mean squared error of the
   estimator.}

   \item{\code{"FLCI"}}{Length of (fixed-length) two-sided confidence
       intervals.}

   \item{\code{"Valid"}}{Optimal estimator under valid moments. This returns
   the original estimator, with confidence intervals adjusted for possible
   misspecification}

    }}
}
\value{
Object of class \code{"GMMEstimate"}, which is a list with at least
    the following components:

\describe{

\item{h}{Point estimate}

\item{bias}{Worst-case bias of estimator}

\item{se}{Standard error of estimator}

\item{hl}{Half-length of confidence interval, so that the confidence interval
takes the form \eqn{h +- hl}}

}
}
\description{
Computes the optimal sensitivity and the optimal estimator when the set
\eqn{\mathcal{C}} takes the form \eqn{c=B\gamma}{c=B*gamma} with the
\eqn{\ell_p}{lp} norm of \eqn{\gamma}{gamma} bounded by \eqn{M}.
}
\examples{
## Replicates estimates in first line of Figure 1 in Armstrong and Kolesár
## (2020)
## 1. Compute matrix B when all instruments are invalid
I <- vector(mode="logical", length=nrow(blp$G))
I[c(6:13, 20:31)] <- TRUE
B <- blp$ZZ \%*\% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ)[, I, drop=FALSE]
## 2. Collect initial estimates
blp$k_init <- -drop(blp$H \%*\% solve(crossprod(blp$G, blp$W \%*\% blp$G),
                                    crossprod(blp$G, blp$W)))
eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
           k_init=blp$k_init, h_init= blp$h_init)
OptEstimator(eo, B, M=sqrt(sum(I)), p=2, alpha=0.05, opt.criterion="Valid")
OptEstimator(eo, B, M=sqrt(sum(I)), p=2, alpha=0.05, opt.criterion="FLCI")
}
\references{
{

\cite{Armstrong, T. B., and M. Kolesár (2020): Sensitivity Analysis Using
Approximate Moment Condition Models,
\url{https://arxiv.org/abs/1808.07387}}

}
}
