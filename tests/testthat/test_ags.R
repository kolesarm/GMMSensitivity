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

})

## test_that("Check we match AGS results", {
##     ## TODO: Objective function: 0.42715, Table II and Figure IV
## })



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
