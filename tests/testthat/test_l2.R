context("L2 constraints")

test_that("Replicate initial analysis", {

    ## Soonwoo's analysis with Omega=W^{-1}
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
    ## Note soonwoo forgot to include 31 in supply and in all

    ## Replication
    eo <- list(H=blp$H,
               G=blp$G,
               Sig=solve(blp$W),
               n=blp$n,
               g_init=blp$g_init,
               h_init=blp$h_init)

    res <- data.frame()
    for (j in 1:length(ivlist)) {
        I <- vector(mode="logical", length=31)
        I[ivlist[[j]]] <- TRUE
        B <- (abs(blp$perturb) * blp$OmZZ)[, I, drop=FALSE]
        K <- sqrt(sum(I))

        opt <- l2opt(eo, B, M=diag(ncol(B)), K, alpha=0.05,
                     opt.criterion="FLCI")
        res <- rbind(res, as.data.frame(opt[-1]))
    }


    ## Should do at least as well as soonwoo
    diff <- res$hl-c(0.0355009384620, 0.0350651386714, 0.0349848074023,
                 0.0363107703323, 0.0369459777305, 0.0351213613424,
                 0.0351985143498, 0.0405290362408, 0.0362155968419,
                 0.0361851150386, 0.0349969735236, 0.0352932950704,
                 0.0351043515741, 0.0379122528181, 0.0364961662495,
                 0.0356789198822, 0.0349998936919, 0.0360211909142,
                 0.0350762430871, 0.0348995214706, 0.0367805352046,
                 0.0465497658943, 0.0403230027212, 0.0392036204477,
                 0.0495483930567, 0.0427226902476, 0.0817657198387)
    expect_lt(max(diff), 0)
    ## Should be approx equal

    diffh <- round(res$h, 3) -
        c(0.336, 0.329, 0.328, 0.351, 0.345, 0.332, 0.314, 0.321, 0.359, 0.355,
          0.33, 0.337, 0.332, 0.405, 0.365, 0.356, 0.33, 0.353, 0.332, 0.327,
          0.358, 0.224, 0.46, 0.431, 0.161, 0.501, 0.374)
    expect_lt(max(abs(diffh)), 0.002)
    ## Lambdas
    diffl <- round(res$kap/(1-res$kap), 3) -
        c(0.993, 0.998, 0.999, 0.993, 0.976, 0.997, 0.996, 0.977, 0.979, 0.979,
          0.999, 0.994, 0.997, 0.999, 0.979, 0.991, 0.999, 0.984, 0.997, 1.015,
          0.997, 0.918, 0.973, 0.98, 0.91, 0.965, 0.347)
    expect_lt(max(abs(diffl)), 0.02)

})
