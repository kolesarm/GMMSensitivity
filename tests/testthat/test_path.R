context("Lp constraints")

test_that("Check l_infty solution path", {

    #' CVX minimization problem for l_infty
    linfbrute <- function(eo, B, lams) {
        ks <- matrix(NA, nrow=length(lams), ncol=nrow(eo$G))
        k <- CVXR::Variable(nrow(eo$G))
        for (j in seq_len(length(lams))) {
            ob <- CVXR::Minimize(CVXR::p_norm(chol(eo$Sig)%*%k)^2/2 +
                                 lams[j]*CVXR::p_norm( t(B) %*% k, p=1))
            pr <- CVXR::Problem(ob, list(-eo$H==t(eo$G)%*%k))
            ks[j, ] <- solve(pr)$getValue(k)
        }
        cbind(lams, ks)
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
        path <- lph(eo, B, p=Inf)
        pathB <- linfbrute(eo, B, path[, 1])
        expect_lt(max(abs(path-pathB)), 0.003)
    }

})



## df <- data.frame(lambda=rep(lams, dg), k=as.vector(res[, 2:(dg+1)]),
##                  what=rep(as.factor(1:dg), each=length(lams)))
## pl1 <- qplot(x=lambda, y=k, color=what, geom="line", data=df)
## pdf("testlassoinf.pdf", width=5, height=4)
## print(directlabels::direct.label(pl1+theme_bw(), "first.qp"))
## dev.off()
