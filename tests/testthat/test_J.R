context("J test")

test_that("Sanity check for J test", {

    ## Construct estimate of initial sensitivity
    blp$k_init <- -drop(blp$H %*%
                        solve(crossprod(blp$G, blp$W %*%
                                               blp$G), crossprod(blp$G, blp$W)))
    eo <- list(H = blp$H, G = blp$G, Sig = blp$Sig, n = blp$n,
               g_init = blp$g_init, k_init = blp$k_init, h_init = blp$h_init)
    ## Matrix B, scaled as described in the paper
    B0 <- blp$ZZ %*% diag(sqrt(blp$n) * abs(blp$perturb)/blp$sdZ)
    ## Value of M

    I <- vector(mode = "logical", length = nrow(eo$G))
    I[6] <- TRUE
    ## With one instrument p, shouldn't matter
    expect_equal(Jtest(eo, B0[, I], M = 1, p = 2, alpha = 0.05)$Mmin,
                 Jtest(eo, B0[, I], M = 1, p = 1, alpha = 0.05)$Mmin)
    eo$g_init <- eo$g_init/5
    ## Rows corresponding to invalid instruments
    I <- vector(mode = "logical", length = nrow(eo$G))
    I[c(6:13, 20:31)] <- TRUE
    ## p0 shouldn't depend on p
    expect_equal(Jtest(eo, B0[, I], M = sqrt(sum(I)), p = 2, alpha = 0.05)$p0,
                 Jtest(eo, B0[, I], M = sqrt(sum(I)), p = 1, alpha = 0.05)$p0)
})
