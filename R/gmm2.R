#' Construct optimal sensitivity and compute worst-case bias and sd for a given
#' kappa
#' @keywords internal
l2sens <- function(eo, B, kap, K, alpha=0.05) {
    ## Optimal inverse weight
    BB <- if (ncol(B)==0) 0 else tcrossprod(B)
    Wi <- (1-kap)*eo$Sig + kap*K^2*BB
    k <- drop(-eo$H %*% solve(crossprod(eo$G, solve(Wi, eo$G)),
                              t(solve(Wi, eo$G))))
    BuildEstimator(k, eo, B, K, p=2, alpha)
}


#' One-step estimator based on optimal sensitivity under ell_2 constraints
#' @keywords internal
l2opt <- function(eo, B, K, alpha=0.05,
                  opt.criterion="FLCI") {

    if (opt.criterion=="MSE") {
        kap <- 1/2
        r <- l2sens(eo, B, kap, K, alpha)
    } else if (opt.criterion=="FLCI") {
        crit <- function(kap) l2sens(eo, B, kap, K, alpha)$hl
        ## minimum should be near 1/2
        kap <- stats::optimize(crit, c(0, 1), tol=1e-12)$minimum
        r <- l2sens(eo, B, kap, K, alpha)
    }
    r$opt.criterion <- opt.criterion
    r$kap <- kap
    r
}


#' Modulus under l2 constraints
#' @keywords internal
l2mod <- function(eo, B, K, delta) {
    ## Calculate deltabar
    ## B' Sigma^-1 B
    BSB <- crossprod(B, solve(eo$Sig, B))
    ## WG
    WG <- function(kap) solve(eo$Sig, eo$G) - kap*solve(eo$Sig, B) %*%
        solve( kap*BSB+(1-kap)*diag(ncol(B)), crossprod(B, solve(eo$Sig, eo$G)))
    GWG <- function(kap) crossprod(eo$G, WG(kap))
    k <- function(kap) -drop(WG(kap) %*% solve(GWG(kap), eo$H))
    Vk <- function(kap) drop(crossprod(k(kap), eo$Sig %*% k(kap)))

    ## Check it's possible to have such matrix:
    if (min(eigen(GWG(1))$values) <= 1e4*.Machine$double.eps) {
        deltabar <- 0
    } else {
        den <- solve(BSB, crossprod(B, solve(eo$Sig, eo$G))) %*%
            solve(GWG(1), eo$H)
        deltabar <- sqrt(4*K^2*Vk(1)/sum(den^2))
    }

    del <- function(kap)
        2*(1-kap)*K / kap * sqrt(Vk(kap)/sum(drop(crossprod(B, k(kap)))^2))

    om <- function(kap)
        del(kap)*drop(crossprod(eo$H, solve(GWG(kap), eo$H)))/sqrt(Vk(kap))

    if (delta<max(deltabar, del(1-1e-6))) {
        return(list(omega=sqrt(Vk(1))*delta, domega=sqrt(Vk(1)), kappa=1))
    } else {
        kap <- stats::uniroot(function(kap) del(kap)-delta,
                              c(0, 1-1e-6), tol = .Machine$double.eps^0.5)$root
        return(list(omega=om(kap), domega=sqrt(Vk(kap)), kappa=kap))
    }
}


#' Efficiency bounds under ell_2 constraints
#'
#' Computes the asymptotic efficiency of two-sided fixed-length confidence
#' intervals, as well as the efficiency of one-sided confidence intervals that
#' optimize a given \code{beta} quantile of excess length.
#'
#' The set \eqn{C} takes the form \eqn{B*gamma} where the ell_2 norm of gamma is
#' bounded by K.
#' @export
#' @inheritParams OptEstimator
#' @param beta Quantile of excess length of one-sided confidence interval to
#'     optimize
l2eff <- function(eo, B, K, beta=0.5, alpha=0.05) {
    ## One-sided
    del <- stats::qnorm(1-alpha) + stats::qnorm(1-beta)
    e1 <- l2mod(eo, B, K, del)
    effo <- l2mod(eo, B, K, 2*del)$omega/(e1$omega+del*e1$domega)

    ## Two-sided
    set.seed(42)
    Z <- stats::rnorm(1000)
    dels <- 2*(stats::qnorm(1-alpha)-Z[Z<stats::qnorm(1-alpha)]) # deltas
    oms <- sapply(dels, function(del) l2mod(eo, B, K, del)$omega)
    efft <- (1-alpha)*mean(oms)/ (2*OptEstimator(eo, B, K=K, 2,
                                            alpha=alpha,
                                            opt.criterion="FLCI")$hl*sqrt(eo$n))

    list(onesisded=effo, twosided=efft)
}
