#' Orthogonalized homotopy solution for l_1
#' @param I vector of indicators which instruments are invalid
#' @keywords internal
l1h0 <- function(G, Sig, h, I) {
    dt <- ncol(G)
    dg <- nrow(G)
    GSG <- function(A) crossprod(G[A, ], solve(Sig[A, A], G[A, ]))

    ## Initialize
    mu <- solve(GSG(rep(TRUE, dg)), h)
    k <- drop(-solve(Sig, G %*% mu))
    B <- max(abs(k[I]))
    A <- rep(FALSE, dg)
    A[I] <- abs(k[I])==max(abs(k[I]))
    res <- matrix(c(B, k, A), nrow=1)
    colnames(res) <- c("B", 1:dg, paste0("A", 1:dg))
    joined <- which.max(A)

    while (sum(A) < dg-dt+1 & B>0) {
        ## directions
        SS <- solve(Sig[!A, !A], Sig[!A, A, drop=FALSE])
        mu.d <- drop(solve(GSG(!A), t(G[A,, drop=FALSE]) -
                               crossprod(G[!A, ], SS)) %*% sign(k[A]))
        k.d <- sign(k)
        k.d[!A] <- -solve(Sig[!A, !A], G[!A, ] %*% mu.d) - SS %*% sign(k[A])

        a.d <- drop(Sig %*% k.d+G %*% mu.d)
        a <- drop(Sig %*% k+G %*% mu)

        d1 <- (a/a.d)
        d1[!A | d1<0] <- Inf
        if(joined!=0)
            d1[joined] <- B

        d2 <- rep(Inf, dg)
        d2[k/B>k.d] <- ((B-k)/(1-k.d))[k/B>k.d]
        d2[k/B<=k.d] <- ((B+k)/(1+k.d))[k/B<=k.d]
        d2[A | !I] <- Inf

        d <- min(d2, d1)
        if(d<0)
            error("Taking a negative step")

        B <- B-d
        k <- k-d*k.d
        mu <- mu-d*mu.d

        if (min(d2)<min(d1)) {
            joined <- which.min(d2)
            A[joined] <- TRUE
        } else {
            A[which.min(d1)] <- FALSE
            joined <- 0
        }
        res <- rbind(res, c(B, k, A))
    }
    res
}
