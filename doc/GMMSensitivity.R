## ---- include=FALSE, cache=FALSE-----------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## ---- echo=FALSE---------------------------------------------------------
library("GMMSensitivity")

## ------------------------------------------------------------------------
## Construct estimate of initial sensitivity
blp$k_init <- -drop(blp$H %*% solve(crossprod(blp$G, blp$W %*% blp$G),
                                        crossprod(blp$G, blp$W)))
## list collecting initial estimates of H, G, Sigma, n, g(thetahat), initial
## sensitivity k, and initial estimate of average markup h(thetahat)
eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
           k_init=blp$k_init, h_init= blp$h_init)

## Rows corresponding to invalid instruments
I <- vector(mode="logical", length=nrow(eo$G))
I[c(6:13, 20:31)] <- TRUE

## Matrix B, scaled as described in the paper
B0 <- blp$ZZ %*% diag(sqrt(blp$n)*abs(blp$perturb)/blp$sdZ)

## Value of M
M0 <- sqrt(sum(I))

## Select columns of B0 corresponding to invalid instruments
OptEstimator(eo, B0[, I], M=M0, p=2, alpha=0.05, opt.criterion="FLCI")

## ------------------------------------------------------------------------
EffBounds(eo, B0[, I], M=M0, p=2)$twosided

## ------------------------------------------------------------------------
OptEstimator(eo, B0[, I], M=M0, p=2, alpha=0.05, opt.criterion="Valid")

## ------------------------------------------------------------------------
Jtest(eo, B0[, I], M=M0, p=2, alpha=0.05)

## ------------------------------------------------------------------------
I <- vector(mode="logical", length=nrow(eo$G))
I[6] <- TRUE
OptEstimator(eo, B0[, I, drop=FALSE], M=1, p=2, alpha=0.05, opt.criterion="FLCI")
EffBounds(eo, B0[, I, drop=FALSE], M=1, p=2)$twosided
OptEstimator(eo, B0[, I, drop=FALSE], M=1, p=2, alpha=0.05, opt.criterion="Valid")
Jtest(eo, B0[, I, drop=FALSE], M=1, p=2, alpha=0.05)

