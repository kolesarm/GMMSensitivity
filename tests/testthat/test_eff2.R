context("L2 Efficiency")

test_that("Check l_2 efficiency calculations using BLP data", {

    ## Alternative efficiency calculation

    #' Modulus under l2 constraints
    #' @keywords internal
    l2mod <- function(eo, B, K, delta) {
        ## Calculate deltabar
        ## B' Sigma^-1 B
        BSB <- crossprod(B, solve(eo$Sig, B))
        ## WG
        WG <- function(kap) solve(eo$Sig, eo$G) - kap*solve(eo$Sig, B) %*%
                                solve(kap*BSB+(1-kap)*diag(ncol(B)),
                                      crossprod(B, solve(eo$Sig, eo$G)))
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

        del2 <- function(kap)
            2*(1-kap)*K / kap * sqrt(Vk(kap)/sum(drop(crossprod(B, k(kap)))^2))

        om2 <- function(kap)
            del2(kap)*drop(crossprod(eo$H, solve(GWG(kap), eo$H)))/sqrt(Vk(kap))

        if (delta<max(deltabar, del2(1-1e-6))) {
            return(list(omega=sqrt(Vk(1))*delta, domega=sqrt(Vk(1)), kappa=1))
        } else {
            kap <- stats::uniroot(function(kap) del2(kap)-delta,
                                  c(0, 1-1e-6), tol=tol)$root
            return(list(omega=om2(kap), domega=sqrt(Vk(kap)), kappa=kap))
        }
    }


    #' Efficiency bounds under ell_2 constraints
    l2eff <- function(eo, B, K, beta=0.5, alpha=0.05) {
        ## One-sided
        zal <- stats::qnorm(1-alpha)
        del <- zal + stats::qnorm(beta)
        e1 <- l2mod(eo, B, K, del)
        effo <- l2mod(eo, B, K, 2*del)$omega/(e1$omega+del*e1$domega)

        integrand <- function(z)
            sapply(z, function(z) l2mod(eo, B, K, 2*(zal-z))$omega *
                                                            stats::dnorm(z))
        lo <- -zal                          # lower endpoint
        while(integrand(lo)>1e-10) lo <- 2*lo

        den <- 2*OptEstimator(eo, B, K=K, 2, alpha=alpha,
                              opt.criterion="FLCI")$hl*sqrt(eo$n)
        eff2 <- stats::integrate(integrand, lo,
                                 zal, abs.tol=1e-6)$value / den


        list(onesisded=effo, twosided=eff2)
    }


    ## List of different specifications
    excluded <- c(6:13, 20:31)
    ivlist <- excluded
    names(ivlist) <- blp$names$iv[excluded]
    ivlist <- as.list(ivlist)
    ivlist$demand_firm <- 6:9
    ivlist$demand_rival <- 10:13
    ivlist$supply_firm <- 20:25
    ivlist$supply_rival <- 26:30
    ivlist$demand <- c(ivlist$demand_firm, ivlist$demand_rival)
    ivlist$supply <- c(ivlist$supply_firm, ivlist$supply_rival, 31)
    ivlist$all <- excluded

    eo <- list(H=blp$H,
               G=blp$G,
               Sig=blp$Sig,             #solve(blp$W)
               n=blp$n,
               g_init=blp$g_init,
               h_init=blp$h_init)

    for (j in 11:length(ivlist)) {
        I <- vector(mode="logical", length=31)
        I[ivlist[[j]]] <- TRUE
        B <- (abs(blp$perturb) * blp$OmZZ)[, I, drop=FALSE]
        K <- sqrt(sum(I))
        eb <- unlist(EffBounds(eo, B, K, p=2, beta=0.5, alpha=0.05))
        e2 <- unlist(l2eff(eo, B, K, beta=0.5, alpha=0.05))
        expect_lt(max(abs(e2-eb)), 1e-5)
    }
})
