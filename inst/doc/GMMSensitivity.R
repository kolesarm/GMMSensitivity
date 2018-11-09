## ---- include=FALSE, cache=FALSE-----------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## ---- echo=TRUE----------------------------------------------------------
library("GMMSensitivity")
help("blp", type="text")

## ------------------------------------------------------------------------
## Construct estimate of initial sensitivity
blp$k_init <- -drop(blp$H %*% solve(crossprod(blp$G, blp$W %*% blp$G),
                                        crossprod(blp$G, blp$W)))

eo <- list(H=blp$H, G=blp$G, Sig=blp$Sig, n=blp$n, g_init=blp$g_init,
           k_init=blp$k_init, h_init= blp$h_init)

## Rows corresponding to invalid instruments
I <- vector(mode="logical", length=nrow(eo$G))
I[c(6:13, 20:31)] <- TRUE

## Matrix B, scaled as described in the paper
B <- (abs(blp$perturb) * blp$OmZZ)[, I, drop=FALSE]

## Value of K
K0 <- sqrt(sum(I))

OptEstimator(eo, B, K=K0, p=2, alpha=0.05, opt.criterion="FLCI")

## ------------------------------------------------------------------------
EffBounds(eo, B, K=K0, p=2)$twosided

## ------------------------------------------------------------------------
OptEstimator(eo, B, K=K0, p=2, alpha=0.05, opt.criterion="Valid")

## ------------------------------------------------------------------------
Jtest(eo, B, K=K0, p=2, alpha=0.05)

