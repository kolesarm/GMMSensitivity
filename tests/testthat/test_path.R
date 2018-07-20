context("L1 and LInf constraints")

test_that("Check l_infty and l_1 solution paths using BLP data", {

    #' CVX minimization problem for l_infty
    linfbrute <- function(eo, B, lams) {
        ks <- matrix(NA, nrow=length(lams), ncol=nrow(eo$G))
        k <- CVXR::Variable(nrow(eo$G))
        for (j in seq_len(length(lams))) {
            ob <- CVXR::Minimize(CVXR::p_norm(chol(eo$Sig)%*%k)^2/2 +
                                 lams[j]*CVXR::p_norm(t(B) %*% k, p=1))
            pr <- CVXR::Problem(ob, list(-eo$H==t(eo$G)%*%k))
            ks[j, ] <- solve(pr)$getValue(k)
        }
        cbind(lams, ks)
    }

    #' CVX minimization problem for l_1
    #' @param Bbound bound on bias
    l1brute <- function(eo, B, Bbound) {
        ks <- matrix(NA, nrow=length(Bbound), ncol=nrow(eo$G))
        k <- CVXR::Variable(nrow(eo$G))
        for (j in seq_len(length(Bbound))) {
            ob <- CVXR::Minimize(CVXR::p_norm(chol(eo$Sig)%*%k)^2/2)
            pr <- CVXR::Problem(ob, list(-eo$H==t(eo$G)%*%k,
                                         CVXR::p_norm(t(B) %*% k, p=Inf) <=
                                         Bbound[j]))
            ks[j, ] <- solve(pr)$getValue(k)
        }
        cbind(Bbound, ks)
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

    res <- data.frame()
    for (j in 1:length(ivlist)) {
        I <- vector(mode="logical", length=31)
        I[ivlist[[j]]] <- TRUE
        B <- (abs(blp$perturb) * blp$OmZZ)[, I, drop=FALSE]
        K <- sqrt(sum(I))
        pathIo <- lph(eo, B, p=Inf)
        pathIb <- linfbrute(eo, B, pathIo[, 1])
        path1o <- lph(eo, B, p=1)
        path1b <- l1brute(eo, B, path1o[, 1])

        expect_lt(max(abs(pathIo[, -ncol(pathIo)]-pathIb)), 0.003)
        expect_lt(max(abs(path1o[, -ncol(path1o)]-path1b)), 0.002)
    }

})

test_that("Check optimal path under no misspecification", {
    I <- vector(mode="logical", length=31)
    B <- (abs(blp$perturb) * blp$OmZZ)[, I, drop=FALSE]
    W <- solve(blp$Sig)
    k_opt <- -blp$H %*% solve(crossprod(blp$G, W %*% blp$G),
                              crossprod(blp$G, W))
    r2 <- OptEstimator(eo, B, K=1, p=2, alpha=0.05, opt.criterion="FLCI")
    r1 <- OptEstimator(eo, B, K=2, p=1, alpha=0.05, opt.criterion="FLCI")
    rI <- OptEstimator(eo, B, K=1, p=Inf, alpha=0.05, opt.criterion="MSE")
    expect_equal(k_opt, r1$k)
    expect_equal(k_opt, r2$k)
    expect_equal(k_opt, rI$k)
})

test_that("Drop invalid instrument under large misspecification", {
    B <- matrix(0, nrow=31)

    B[6] <- 1
    kI <- OptEstimator(eo, B, K=500, p=Inf, alpha=0.05, opt.criterion="FLCI")$k
    k2 <- OptEstimator(eo, B, K=500, p=2, alpha=0.05, opt.criterion="MSE")$k
    k1 <- OptEstimator(eo, B, K=500, p=1, alpha=0.05, opt.criterion="FLCI")$k

    W <- solve(blp$Sig[-6, -6])
    k_opt <- drop(-blp$H %*% solve(crossprod(blp$G[-6, ], W %*% blp$G[-6, ]),
                     crossprod(blp$G[-6, ], W)))
    expect_lt(max(abs(k_opt- kI[-6])), 1e-4)
    expect_lt(max(abs(k_opt- k1[-6])), 1e-4)
    expect_lt(max(abs(k_opt- k2[-6])), 1e-4)
})
