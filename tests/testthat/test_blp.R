context("Check BLP analysis")

test_that("Check we match AGS results and modulus matches brute force", {
    excluded <- c(6:13, 20:31)
    ivlist <- c(list(None=vector(mode="numeric", length=0)), as.list(excluded),
            list("All D/F"=6:9, "All D/R"=10:13, "All S/F"=20:25,
                 "All S/R"=26:30, "All excluded demand"=c(6:13),
                 "All excluded supply"=c(20:31), "All excluded"=c(6:13, 20:31)))

    blp$k_init <- -drop(blp$H %*% solve(crossprod(blp$G, blp$W %*% blp$G),
                                        crossprod(blp$G, blp$W)))
    eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
               k_init=blp$k_init, h_init= blp$h_init)

    ## 0.1 Replicate AGS Table II (for blp.tex writeup)
    if(FALSE) {
        t2i <- t2o <- data.frame()
        B1 <- blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb))
        for (j in c(1, 10, 16, 2, 6)) {
            I <- vector(mode="logical", length=nrow(eo$G))
            I[ivlist[[j]]] <- TRUE
            oe <- function(oc) {
                r <- data.frame(OptEstimator(eo, B1[, I, drop=FALSE], 1, 2,
                                             alpha=0.05, opt.criterion=oc)[-1])
                ## Convert to %
                r[c(1:4, 6)] <- 100*r[c(1:4, 6)]
                r
            }
            t2i <- rbind(t2i, oe("Valid"))
            t2o <- rbind(t2o, oe("FLCI"))
        }
        t2i$name <- t2o$name <- names(ivlist[c(1, 10, 16, 2, 6)])
        t2i$group <- "Initial"
        t2o$group <- "Optimal"
        t2o$kap <- NULL
        t2 <- rbind(t2i, t2o)
        print(t2)
    }
    ## Let's back out the biases instead from Figure IV

    ## 1. M=1, graphs, adaptation, Figure IV in AGS ################
    romI <- rom1 <- rom2 <- rofI <- rof1 <- rof2 <- rinI <- rin1 <- rin2 <-
        effh <- effc <- data.frame()

    B0 <- blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ)
    ## romI/rof1/rin: results for optimal CI/optimal MSE estimator/initial
    ## estimator under l1, l2, l_Inf norm; effh/effc: efficiency computed using
    ## cvx/homotopy
    for (j in seq_along(ivlist)) {
        I <- vector(mode="logical", length=nrow(eo$G))
        I[ivlist[[j]]] <- TRUE
        Mp <- function(p)
            if (p==Inf) 1 else if (p==1) sum(I) else if (p==2) sqrt(sum(I))
        oe <- function(p, oc)
            as.data.frame(OptEstimator(eo, B0[, I, drop=FALSE], M=Mp(p), p,
                                       alpha=0.05, opt.criterion=oc)[-1])
        romI <- rbind(romI, oe(Inf, "MSE"))
        rofI <- rbind(rofI, oe(Inf, "FLCI"))
        rinI <- rbind(rinI, oe(Inf, "Valid"))
        rom1 <- rbind(rom1, oe(1, "MSE"))
        rof1 <- rbind(rof1, oe(1, "FLCI"))
        rin1 <- rbind(rin1, oe(1, "Valid"))
        rom2 <- rbind(rom2, oe(2, "MSE"))
        rof2 <- rbind(rof2, oe(2, "FLCI"))
        rin2 <- rbind(rin2, oe(2, "Valid"))

        ## Efficiency of adaptation
        if (j>15) {
            eb <- function(p, cvx) EffBounds(eo, B0[, I, drop=FALSE], Mp(p), p,
                                             beta=0.8, alpha=0.05, cvx=cvx)
            effc <- rbind(effc, data.frame(eb(1, TRUE), eb(2, TRUE),
                                           eb(Inf, TRUE)))
            effh <- rbind(effh, data.frame(eb(1, FALSE), eb(2, FALSE),
                                           eb(Inf, FALSE)))
        }
    }
    sub <- c(10, 16, 2, 6)
    c2 <- rin2[sub, "bias"]*blp$sdZ[do.call(c, ivlist[sub])]
    c1 <- rin1[sub, "bias"]*blp$sdZ[do.call(c, ivlist[sub])]
    cI <- rinI[sub, "bias"]*blp$sdZ[do.call(c, ivlist[sub])]
    expect_equal(c2, c1)
    expect_equal(cI, c1)
    ## Check against Table II in AGS
    expect_equal(round(c2, 4), c(0.1731, 0.2095, 0.1277, 0.2515))
    ## Check under l_infty, biases add up
    expect_equal(sum(rinI[2:5, "bias"]), rinI[22, "bias"])
    expect_equal(sum(rinI[2:9, "bias"]), rinI[26, "bias"])
    ## Check under l_1, l_2, sum of individual biases should be smaller
    expect_lt(sum(rin1[2:9, "bias"]), rin1[26, "bias"])
    expect_lt(sum(rin2[2:9, "bias"]), rin2[26, "bias"])
    expect_lt(sum(rin1[2:5, "bias"]), rin1[22, "bias"])
    expect_lt(sum(rin2[2:5, "bias"]), rin2[22, "bias"])

    ## Adaptation Efficiency
    expect_lt(max(abs(effh-effc)), 4e-5)

})
