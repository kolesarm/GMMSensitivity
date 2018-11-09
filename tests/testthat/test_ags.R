context("AGS (2017) replication")

test_that("Replicate Table 2 and Figure 2 in AGS", {

excinst <- c(6:13, 20:31)
## Load data
Lambda <- -solve(crossprod(blp$G, blp$W %*% blp$G), crossprod(blp$G, blp$W))

## Figure 2
figIV <- (blp$H %*% Lambda %*% blp$OmZZ * abs(blp$perturb) / blp$sdZ)[excinst]
## Table 2
tblIV <- (blp$H %*% Lambda %*% blp$OmZZ * blp$perturb)[excinst]
names(figIV) <- names(tblIV) <- blp$names$iv[excinst]

## Figuere 2 and Table 2
figIV
subs <- c("supply_firm_const", "supply_rival_const",
          "demand_firm_const", "demand_rival_const")

expect_equal(unname(round(tblIV[subs], 4)),
             c(-0.1731, 0.2095, -0.1277, 0.2515))
expect_lt(max(abs(figIV)), 0.15)

})
