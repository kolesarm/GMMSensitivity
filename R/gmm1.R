#' Orthogonalized homotopy solution for l_1
#' @param I vector of indicators which instruments are invalid
#' @keywords internal
l1h0 <- function(G, Sig, h, I) {
    dt <- ncol(G)
    dg <- nrow(G)
    res <- matrix(0, ncol=2*dg+1, nrow=1)
    colnames(res) <- c("B", 1:dg, paste0("A", 1:dg))

    GSG <- function(A) crossprod(G[A, ], solve(Sig[A, A], G[A, ]))

    ## Initialize
    mu <- solve(GSG(rep(TRUE, dg)), h)
    k <- drop(-solve(Sig, G %*% mu))
    B <- max(abs(k[I]))
    A <- rep(FALSE, dg)
    A[I] <- abs(k[I])==max(abs(k[I]))
    res[1, ] <- c(B, k, A)
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
        B <- B-d
        k <- k-d*k.d
        mu <- mu-d*mu.d

        if (min(d2)<min(d1)){
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

## l1brute <- function(G, Sig, h, B, I) {
##     ks <- matrix(NA, nrow=length(B), ncol=nrow(G))
##     k <- CVXR::Variable(nrow(G))
##     for (j in seq_len(length(B))) {
##         ob <- CVXR::Minimize(CVXR::p_norm(chol(Sig)%*%k)^2/2)
##         pr <- CVXR::Problem(ob, list(-h==t(G)%*%k,
##                                      CVXR::p_norm(k[I], p=Inf) <= B[j]))
##         ks[j, ] <- solve(pr)$getValue(k)
##     }
##     cbind(B, ks)
## }


## ## L_1 penalized solution path for optimal sensitivities
## set.seed(42)
## dg <- 8
## nv <- 1
## I <- c(rep(TRUE, dg-nv), rep(FALSE, nv))
## dt <- 3
## G <- matrix(rnorm(dg*dt), nrow=dg)
## Sig <- cov(matrix(rnorm(2*dg*dg), ncol=dg))
## h <- c(1, rep(0, dt-1))


## res <- l1h(G, Sig, h, I)[, 1:(dg+1)]
## Bs <- res[, 1]
## resb <- l1brute(G, Sig, h, Bs, I)
## diff <- max(abs(res[, 1:(dg+1)]-resb))
## if(diff>1e-3) message("Big difference\n") else print(diff)

## df <- data.frame(B=rep(Bs, dg), k=as.vector(res[, 2:(dg+1)]),
##                  what=rep(as.factor(1:dg), each=length(Bs)))
## pl1 <- qplot(x=B, y=k, color=what, geom="line", data=df)
## pdf("testlasso1.pdf", width=5, height=4)
## print(directlabels::direct.label(pl1+theme_bw(), "last.qp"))
## dev.off()
