#' Construct optimal sensitivity and compute worst-case bias and sd for a given
#' kappa
#' @keywords internal
l2sens <- function(eo, B, M=diag(ncol(B)), kap, K, alpha=0.05) {
    ## Optimal inverse weight
    BMMB <- if (ncol(B)==0) 0 else B %*% solve(crossprod(M), t(B))
    Wi <- (1-kap)*eo$Sig + kap*K^2*BMMB
    k <- drop(-eo$H %*% solve(crossprod(eo$G, solve(Wi, eo$G)),
                              t(solve(Wi, eo$G))))
    BuildEstimator(k, eo, B, M, K, p=2, alpha)
}


#' One-step estimator based on optimal sensitivity under ell_2 constraints
#' @keywords internal
l2opt <- function(eo, B, M=diag(ncol(B)), K, alpha=0.05,
                  opt.criterion="FLCI") {

    if (opt.criterion=="MSE") {
        kap <- 1/2
        r <- l2sens(eo, B, M, kap, K, alpha)
    } else if (opt.criterion=="FLCI") {
        crit <- function(kap) l2sens(eo, B, M, kap, K, alpha)$hl
        ## minimum should be near 1/2
        kap <- stats::optimize(crit, c(0, 1), tol=1e-12)$minimum
        r <- l2sens(eo, B, M, kap, K, alpha)
    }
    r$opt.criterion <- opt.criterion
    r$kap <- kap
    r
}


#' Modulus
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
        return(list(domega=sqrt(Vk(1)), omega=sqrt(Vk(1))*delta, kappa=1))
    } else {
        kap <- stats::uniroot(function(kap) del(kap)-delta,
                              c(0, 1-1e-6))$root
        return(list(domega=sqrt(Vk(kap)), omega=om(kap), kappa=kap))
    }
}


#' old
l2modo <- function(eo, B, kap, K) {

    BB <- tcrossprod(B)
    SG <- solve(eo$Sig, eo$G)
    GSG <- crossprod(eo$G, SG)
    ## W*G
    if (kap > 0.1) {
        WG <- SG - solve(eo$Sig, B) %*%
            solve((1-kap)/(kap*K)* diag(ncol(B)) +
                  crossprod(B, solve(eo$Sig, B)), crossprod(B, SG))
    } else {
        WG <- solve((1-kap)*eo$Sig + kap*K^2*BB, eo$G)
    }
    k <- drop(-eo$H %*% solve(crossprod(eo$G, WG), t(WG)))
    kv <- drop(-eo$H %*% solve(GSG, t(SG)))
    kSk <-  drop(crossprod(k,  eo$Sig %*% k))
    kSkv <- drop(crossprod(kv, eo$Sig %*% kv))
    ## Worst-case bias occurs at:
    cc <- -K * drop(B%*% crossprod(B, k) /
                    sqrt(drop(crossprod(crossprod(B, k)))))
    Sc <- solve(eo$Sig, cc)
    del <- drop(crossprod(cc, Sc)) -
        drop(crossprod(crossprod(eo$G, Sc), solve(GSG, crossprod(eo$G, Sc))))

    del <- 2 * sqrt(del / (1-kSkv/kSk))
    om <- kSkv * del / sqrt(kSk) - 2 * sum(cc*kv)
    list(delta=del, omega=om, domega=sqrt(kSk))
}

#' efficiency
#' @param tol Search in interval [tol, 1-tol]
l2eff <- function(eo, B, K, tu=1e-6, beta=0.5, alpha=0.05) {
    del <- stats::qnorm(1-alpha) + stats::qnorm(1-beta)
    kap1 <- stats::uniroot(function(kap)
        l2mod(eo, B, kap, K=K)$delta-del, c(0, 1-tu))$root
    kap2 <- stats::uniroot(function(kap)
        l2mod(eo, B, kap, K=K)$delta-2*del, c(0, 1-tu))$root
    e1 <- l2mod(eo, B, kap1, K=K)
    effo <- l2mod(eo, B, kap2, K=K)$omega / (e1$omega+del*e1$domega)

    ## Two-sided
    set.seed(42)
    Z <- rnorm(1000)
    Z <- Z[Z<qnorm(0.95)]
    oms <- vector(length=length(Z))
    for (j in seq_along(oms)) {
        del <- 2*(qnorm(1-alpha)-Z[j])
        kap <- stats::uniroot(function(kap)
        l2mod(eo, B, kap, K=K)$delta-del, interval=c(0, 1-tu))$root
        oms[j] <- l2mod(eo, B, kap, K=K)$omega
    }
    efft <- 0.95*mean(oms)/ (2*OptEstimator(eo, B, M=diag(ncol(B)), K=K, 2,
                                            alpha=alpha,
                                            opt.criterion="FLCI")$hl*sqrt(eo$n))

    list(onesisded=effo, twosided=efft)
}
