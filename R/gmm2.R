#' Construct optimal sensitivity and compute worst-case bias and sd for a given
#' kappa
#' @keywords internal
l2sens <- function(eo, B, kap, K, alpha=0.05) {
    ## Optimal inverse weight
    BB <- if (ncol(B)==0) 0 else tcrossprod(B)
    Wi <- (1-kap)*eo$Sig + kap*K^2*BB
    k <- drop(-eo$H %*% solve(crossprod(eo$G, solve(Wi, eo$G)),
                              t(solve(Wi, eo$G))))
    BuildEstimator(k, eo, B, K, p=2, alpha)
}


#' One-step estimator based on optimal sensitivity under ell_2 constraints
#' @keywords internal
l2opt <- function(eo, B, K, alpha=0.05,
                  opt.criterion="FLCI") {

    if (opt.criterion=="MSE") {
        kap <- 1/2
        r <- l2sens(eo, B, kap, K, alpha)
    } else if (opt.criterion=="FLCI") {
        crit <- function(kap) l2sens(eo, B, kap, K, alpha)$hl
        ## minimum should be near 1/2
        kap <- stats::optimize(crit, c(0, 1), tol=tol)$minimum
        r <- l2sens(eo, B, kap, K, alpha)
    }
    r$opt.criterion <- opt.criterion
    r$kap <- kap
    r
}
