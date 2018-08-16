#' Modulus
modulus <- function(eo, B, K, p=2, spath=NULL, delta) {
    ## drop lambda/barB and #dropped moments
    if (p != 2) {
        if (is.null(spath))
            spath <- lph(eo, B, p)
        spath <- spath[, -c(1, ncol(spath)), drop=FALSE]
    }

    Vk <- function(k) drop(crossprod(k, eo$Sig %*% k))
    SG <- solve(eo$Sig, eo$G)       # Sigma^{-1} * G
    kv <- -drop(SG %*% solve(crossprod(eo$G, SG), eo$H))

    ## c_delta
    cd <- function(k) {
        Bk <- drop(crossprod(B, k))
        if (p==Inf) {
            -K*drop(B %*% sign(Bk))
        } else if (p==1) {
            j <- which.max(abs(Bk))
            -sign(Bk[j]) *K*drop(B %*% (1:ncol(B)==j))
        } else {
            -K*drop(B %*% Bk)/ sqrt(sum(Bk^2))
        }
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
        (1-ix+floor(ix))*spath[floor(ix), ]+(ix-floor(ix))*spath[ceiling(ix), ]
    }

    kapmax <- 1
    ## For p==2, if GW2G(1) is not invertible, set kappamax=1-1e-6
    if (p==2 & min(eigen(GW2G(1))$values) <= 1e4*.Machine$double.eps)
        kapmax <- 1-1e-6
    ## delta needs to be decreasing in kappa
    eps <- 1e-6
    while (del(kp(kapmax)) > del(kp(kapmax-eps))) {
        kapmax <- kapmax-eps
        eps <- 2*eps
    }
    if (delta < del(kp(kapmax))) {
        do <- sqrt(Vk(kp(kapmax)))
        return(list(omega=do*delta, domega=do, kappa=kapmax))
    } else {
        kap <- stats::uniroot(function(kap) del(kp(kap))-delta,
                              c(0, kapmax), tol=tol)$root
        return(list(omega=om(kp(kap)), domega=sqrt(Vk(kp(kap))), kappa=kap))
    }
}


#' Efficiency bounds under ell_p constraints
#'
#' Computes the asymptotic efficiency of two-sided fixed-length confidence
#' intervals, as well as the efficiency of one-sided confidence intervals that
#' optimize a given \code{beta} quantile of excess length.
#'
#' The set \eqn{C} takes the form \eqn{B*gamma} where the ell_p norm of gamma is
#' bounded by K.
#' @inheritParams OptEstimator
#' @param beta Quantile of excess length of one-sided confidence interval to
#'     optimize
#' @export
EffBounds <- function(eo, B, K, p=2, beta=0.5, alpha=0.05) {
    ## One-sided
    zal <- stats::qnorm(1-alpha)
    del0 <- zal + stats::qnorm(1-beta)
    spath <- if (p!=2) lph(eo, B, p) else NULL
    e1 <- modulus(eo, B, K, p, spath, del0)
    eff1 <- modulus(eo, B, K, p, spath, 2*del0)$omega/(e1$omega+del0*e1$domega)

   integrand <- function(z)
       sapply(z, function(z) modulus(eo, B, K, p, spath,
                                     2*(zal-z))$omega * stats::dnorm(z))
    lo <- -zal                          # lower endpoint
    while(integrand(lo)>1e-10) lo <- 2*lo

    den <- 2*OptEstimator(eo, B, K=K, p, alpha=alpha,
                          opt.criterion="FLCI")$hl*sqrt(eo$n)
    eff2 <- stats::integrate(integrand, lo,
                             zal, abs.tol=1e-6)$value / den

    list(onesisded=eff1, twosided=eff2)
}
