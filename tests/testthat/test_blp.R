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
    B0 <- blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ)

    rpoc <- function(p, oc) {
        res <- data.frame()
        for (j in seq_along(ivlist)) {
            I <- vector(mode="logical", length=nrow(eo$G))
            I[ivlist[[j]]] <- TRUE
            Mp <- switch(as.character(p), "Inf"=1, "1"=sum(I), "2"=sqrt(sum(I)))
            oe <- function(p, oc)
                as.data.frame(OptEstimator(eo, B0[, I, drop=FALSE], M=Mp, p,
                                       alpha=0.05, opt.criterion=oc)[-1])
            res <- rbind(res, oe(p, oc))
        }
        res
    }
    rofI <- rpoc(p=Inf, oc="FLCI")
    rof2 <- rpoc(p=2, oc="FLCI")
    rof1 <- rpoc(p=1, oc="FLCI")
    rinI <- rpoc(p=Inf, oc="Valid")
    rin2 <- rpoc(p=2, oc="Valid")
    rin1 <- rpoc(p=1, oc="Valid")
    romI <- rpoc(p=Inf, oc="MSE")
    rom2 <- rpoc(p=2, oc="MSE")
    rom1 <- rpoc(p=1, oc="MSE")

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
})
