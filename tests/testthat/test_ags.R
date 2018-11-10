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
