context("Lp Efficiency")

test_that("Check l_1, l_2, l_inf efficiency agree with one invalid moment", {
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

    eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
               h_init=blp$h_init)

    for (j in 1:10) {
        I <- vector(mode="logical", length=31)
        I[ivlist[[j]]] <- TRUE
        B <- (abs(blp$perturb) * blp$ZZ)[, I, drop=FALSE]
        e1 <- unlist(EffBounds(eo, B, 1, p=1, beta=0.8, alpha=0.05))
        e2 <- unlist(EffBounds(eo, B, 1, p=2, beta=0.8, alpha=0.05))
        eI <- unlist(EffBounds(eo, B, 1, p=Inf, beta=0.8, alpha=0.05))
        expect_lt(max(abs(e1-e2)), 1e-6)
        expect_lt(max(abs(eI-e2)), 1e-6)
    }
})


test_that("Check modulus solution against brute force", {
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

    eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
               h_init=blp$h_init)

    dels <- c(seq(0, 2, by=0.3), 4, 6)
    for(p in c(1, 2, Inf)) {
        for (j in c(11, 21:27)) {
            I <- vector(mode="logical", length=31)
            I[ivlist[[j]]] <- TRUE
            B <- (abs(blp$perturb) * blp$ZZ)[, I, drop=FALSE]
            M <- if (p==2) sqrt(sum(I)) else if (p==Inf) 1 else sum(I)
            mb <- vapply(dels, function(d) mod_cvx(d, eo, B, M, p)$omega,
                         numeric(1))
            mo <- vapply(dels, function(d) modulus(d, eo, B, M, p)$omega,
                         numeric(1))
            expect_lt(max(abs(mb-mo)), 5e-6)
        }
    }
})
