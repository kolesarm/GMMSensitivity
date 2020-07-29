#' Efficiency bounds under \eqn{\ell_p}{ell_p} constraints
#'
#' Computes the asymptotic efficiency of two-sided fixed-length confidence
#' intervals at \eqn{c=0}, as well as the efficiency of one-sided confidence
#' intervals that optimize a given \code{beta} quantile of excess length, when
#' the set \eqn{\mathcal{C}} is characterized by \eqn{\ell_p} constraints.
#'
#' The set \eqn{\mathcal{C}} takes the form \eqn{B\gamma}{B*gamma} where the
#' \eqn{\ell_p}{lp} norm of \eqn{\gamma}{gamma} is bounded by \eqn{M}.
#' @inheritParams OptEstimator
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
#'     \item{g_init}{Moment condition evaluated at initial estimate}
#'
#'     }
#' @param beta Quantile of excess length that a one-sided confidence interval is
#'     optimizing.
#' @return A list with two elements, \code{"onesided"} for efficiency of
#'     one-sided CIs and \code{"twosided"} for efficiency of two-sided CIs
#' @examples
#' ## Replicates first line of Table 2 in Armstrong and Kolesár (2020)
#' ## First compute matrix B
#' I <- vector(mode="logical", length=nrow(blp$G))
#' I[6] <- TRUE
#' B <- (blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ))[, I, drop=FALSE]
#' eo <- list(H=blp$H, G=blp$G, Sig=solve(blp$W), n=blp$n, g_init=blp$g_init)
#' EffBounds(eo, B, M=1, p=Inf, beta=0.5, alpha=0.05)
#' @references{
#'
#' \cite{Armstrong, T. B., and M. Kolesár (2020): Sensitivity Analysis Using
#' Approximate Moment Condition Models,
#' \url{https://arxiv.org/abs/1808.07387}}
#'
#' }
#' @export
EffBounds <- function(eo, B, M, p=2, beta=0.5, alpha=0.05) {
    ## One-sided
    zal <- stats::qnorm(1-alpha)
    del0 <- zal + stats::qnorm(beta)
    spath <- if (p!=2) lph(eo, B, p) else NULL
    mo <- function(del) modulus(del, eo, B, M, p, spath)

    e1 <- mo(del0)
    eff1 <- mo(2*del0)$omega/(e1$omega+del0*e1$domega)

    integrand <- function(z)
        vapply(z, function(z) mo(2*(zal-z))$omega * stats::dnorm(z),
               numeric(1))
    lo <- -zal                          # lower endpoint
    while(integrand(lo)>1e-10) lo <- lo-2

    den <- 2*OptEstimator(eo, B, M=M, p, alpha=alpha,
                          opt.criterion="FLCI")$hl*sqrt(eo$n)
    eff2 <- stats::integrate(integrand, lo,
                             zal, abs.tol=1e-6)$value / den

    list(onesisded=eff1, twosided=eff2)
}


## Modulus
modulus <- function(delta, eo, B, M, p=2, spath=NULL) {
    ## drop lambda/barB and #dropped moments
    if (p != 2) {
        if (is.null(spath))
            spath <- lph(eo, B, p)
    }

    Vk <- function(k) drop(crossprod(k, eo$Sig %*% k))
    SG <- solve(eo$Sig, eo$G)       # Sigma^{-1} * G
    kv <- -drop(SG %*% solve(crossprod(eo$G, SG), eo$H))

    ## c_delta. If Bk=0, give c_delta that leads to largest possible delta.
    cd <- function(k) {
        Bk <- drop(crossprod(B, k))
        R <- crossprod(B, solve(eo$Sig, B)-SG %*%
                            solve(crossprod(eo$G, SG), crossprod(SG, B)))
        gam <- rep(NA, length=length(Bk))
        if (rcond(R)>100*.Machine$double.eps) {
            gam <- solve(R, drop(crossprod(B, k-kv)))
        } else {
            if (p==2) {
                gam <- -Bk
            } else if (p==Inf) {
                A <- abs(Bk) >= 1e-8
                gam[A] <- -M*sign(Bk[A])
                if (sum(!A)>0) {
                    den <- crossprod(gam[A], R[A, A, drop=FALSE] %*% gam[A]) -
                        crossprod(R[!A, A, drop=FALSE] %*% gam[A],
                                  solve(R[!A, !A, drop=FALSE],
                                        R[!A, A, drop=FALSE] %*% gam[A]))
                    Bkv <- drop(crossprod(B, kv))
                    num <- Vk(k)-Vk(kv) -
                        crossprod(Bkv[!A], solve(R[!A, !A, drop=FALSE],
                                                 Bkv[!A]))
                    lam <- sqrt(drop(num)/drop(den))
                    gam[!A] <- -solve(R[!A, !A, drop=FALSE],
                                      R[!A, A, drop=FALSE] %*% gam[A] +
                                      Bkv[!A]/lam)
                }
            } else if (p==1) {
                A <- (abs(Bk) < max(abs(Bk)) -1e-8)
                gam[A] <- 0
                gam[!A] <- solve(R[!A, !A], drop(crossprod(B, k-kv))[!A])
            }
        }

        den <- switch(as.character(p),
                      "Inf" = max(abs(gam)),
                      "2"=sqrt(sum(gam^2)),
                      "1"=sum(abs(gam)))
        M * drop(B %*% gam) / den
    }

    del <- function(k) {
        Sc <- solve(eo$Sig, cd(k))
        GSc <- crossprod(eo$G, Sc)
        top <- crossprod(cd(k), Sc) -
            crossprod(GSc, solve(crossprod(eo$G, SG), GSc))
        if(!isTRUE(all.equal(Vk(kv), Vk(k)))) {
            2*sqrt(drop(top)/(1-Vk(kv)/Vk(k)))
        } else {
            .Machine$double.xmax        # i.e. Inf, but uniroot dislikes this
        }
    }
    om <- function(k)
        Vk(kv)/sqrt(Vk(k)) * del(k) - 2*sum(kv*cd(k))

    ## W(kappa) * G for ell_2 constraints
    W2G <- function(kap)
        solve(eo$Sig, eo$G) - kap*solve(eo$Sig, B) %*%
            solve(kap*crossprod(B, solve(eo$Sig, B))+(1-kap)*diag(ncol(B)),
                  crossprod(B, solve(eo$Sig, eo$G)))
    GW2G <- function(kap) crossprod(eo$G, W2G(kap))

    ## k_kappa, where kappa is either parametrizing the weight matrix, or else
    ## linearly interpolates the solution path
    kp <- function(kap) {
        if (p==2)
            return(-drop(W2G(kap) %*% solve(GW2G(kap), eo$H)))
        ix <- 1+kap*(nrow(spath)-1)
        (ix-floor(ix))*spath[ceiling(ix), ]+(1-ix+floor(ix))*spath[floor(ix), ]
    }

    kapmax <- 1
    ## If GW2G(1) is not invertible, set kappamax=1-1e-6, the endpoint is never
    ## reached, of for p\neq 2, the second-to-last step in the solution path
    if (rcond(GW2G(1)) < 100*.Machine$double.eps) {
        if (p==2) {
            kapmax <- 1-1e-6
        } else {
            kapmax <- (nrow(spath)-2)/(nrow(spath)-1)
        }
    }

    if (delta < del(kp(kapmax))) {
        if (kapmax==1) {
            do <- sqrt(Vk(kp(kapmax)))
            return(list(omega=do*delta, domega=do, kappa=kapmax))
        } else {
            ## brute force
            mb <- mod_cvx(delta, eo, B, M, p)
            return(list(omega=mb$omega, domega=mb$domega, kappa=-1))
        }
    } else {
        kap <- stats::uniroot(function(kap) del(kp(kap))-delta,
                              c(0, kapmax), tol=tol)$root

        return(list(omega=om(kp(kap)), domega=sqrt(Vk(kp(kap))), kappa=kap))
    }
}


## Modulus using CVX
mod_cvx <- function(delta, eo, B, M, p) {
    th <- CVXR::Variable(length(eo$H))
    ga <- CVXR::Variable(ncol(B))
    ob <- CVXR::Maximize(2*sum(eo$H*th))
    con <- list(CVXR::p_norm(t(solve(chol(eo$Sig))) %*%
                             (B%*%ga - eo$G%*%th))  <= delta/2,
                CVXR::p_norm(ga, p=p)<=M)
    pr <- CVXR::solve(CVXR::Problem(ob, con))
    th0 <- drop(pr$getValue(th))
    ga0 <- drop(pr$getValue(ga))

    ## k_{delta} up to scale: k_{delta}/lambda_1
    k0 <- solve(eo$Sig, drop(B %*% ga0 - eo$G %*% th0))
    ## Get scale
    lam1 <- mean(eo$H[eo$H!=0]/drop(-crossprod(eo$G, k0))[eo$H!=0])
    ## Alternative: domega=sqrt(drop(crossprod(k0, eo$Sig %*% k0)))*lam1

    list(omega=pr$value, domega=delta*lam1/2)
}
