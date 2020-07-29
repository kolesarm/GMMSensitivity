context("AGS (2017) replication")

test_that("Replicate Table II and Figure IV in AGS", {

    excinst <- c(6:13, 20:31)
    ## Load data
    Lambda <- -solve(crossprod(blp$G, blp$W %*% blp$G), crossprod(blp$G, blp$W))

    ## Figure IV
    figIV <- (blp$H %*% Lambda %*% blp$ZZ * abs(blp$perturb) / blp$sdZ)[excinst]
    ## Table II
    tblII <- (blp$H %*% Lambda %*% blp$ZZ * blp$perturb)[excinst]
    names(figIV) <- names(tblII) <- blp$names$iv[excinst]

    ## Figure IV and Table II
    subs <- c("supply_firm_const", "supply_rival_const",
              "demand_firm_const", "demand_rival_const")
    expect_equal(unname(round(tblII[subs], 4)),
                 c(-0.1731, 0.2095, -0.1277, 0.2515))
    expect_lt(max(abs(figIV)), 0.15)
    expect_equal(unname(which.min(abs(figIV[1:8]))), 7)
    expect_equal(names(which.min(abs(figIV))), "supply_mpd")

    ## Match value of J statistic from AGS code
    I <- vector(mode="logical", length=31)
    I[6] <- TRUE
    B <- (blp$ZZ %*%
          diag(sqrt(blp$n)*abs(blp$perturb) / blp$sdZ))[, I, drop=FALSE]
    blp$k_init <- -drop(blp$H %*% solve(crossprod(blp$G, blp$W %*% blp$G),
                                        crossprod(blp$G, blp$W)))
    eo <- list(H=blp$H, G=blp$G, Sig=solve(blp$W), n=blp$n, g_init=blp$g_init,
               k_init=blp$k_init, h_init= blp$h_init)
    r <- Jtest(eo, B, M=1, p=2, alpha=0.05)
    expect_equal(round(r$J/blp$n, 5), 0.42715)
    ## Match Table 1 in paper
    expect_equal(round(r$Mmin, 4), 10.2055)
})
