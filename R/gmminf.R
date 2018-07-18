#' Compute solution path for l_infty or l_1 constraints
#' @inheritParams OptEstimator
#' @export
lph <- function(eo, B, M=diag(ncol(B)), p=Inf) {
    if (ncol(B)==0)
        return(-eo$H %*% solve(crossprod(eo$G, solve(eo$Sig, eo$G)),
                                    t(solve(eo$Sig, eo$G))))

    ## Get orthogonalized homotopy, first B_{\perp}
    Bp <- if (nrow(B)>ncol(B)) {
              qr.Q(qr(B), complete=TRUE)[, (ncol(B)+1):nrow(B)]
          } else {
              matrix(ncol=0, nrow=nrow(B))
          }
    I <- rep(c(FALSE, TRUE), c(ncol(Bp), ncol(B)))
    T <- rbind(t(Bp), M %*% solve(crossprod(B), t(B)))
    Sigt <- T %*% eo$Sig %*% t(T)
    ## Get path of tilde{k}
    if (p==Inf)
        kts <- linfh0(T %*% eo$G, Sigt, eo$H, I)[, 2:(nrow(B)+1)]
    else
        kts <- l1h0(T %*% eo$G, Sigt, eo$H, I)[, 2:(nrow(B)+1)]
    ## Return sensitivities at each step
    kts %*% T
}

linfstep <- function(s, G, Sig, H, I) {
    dg <- nrow(G)

    d1 <- -s$k/s$k.d
    d1[!s$A | d1<0 | !I] <- Inf
    if(s$joined!=0)
        d1[s$joined] <- Inf

    d2 <- rep(Inf, dg)
    a.d <- drop(Sig %*% s$k.d+G %*% s$mu.d)
    a <- drop(Sig %*% s$k+G %*% s$mu)

    if (sum(a.d[!s$A]>1)>0)
        d2[a.d>1] <- ((s$lam-a)/(a.d-1))[a.d>1]
    if (sum(a.d[!s$A]< -1)>0)
        d2[a.d< -1] <- (-(s$lam+a)/(a.d+1))[a.d< -1]
    d2[s$A] <- Inf

    d <- min(d2, d1)
    s$lam <- s$lam+d
    s$k <- s$k + d*s$k.d
    s$mu <- s$mu + d*s$mu.d

    if (min(d2)<min(d1)){
        s$joined <- which(d2<=min(d2))
        s$A[s$joined] <- TRUE
    } else {
        s$A[which(d1<=min(d1))] <- FALSE
        s$joined <- 0
    }

    GSG <- function(A) crossprod(G[A, ], solve(Sig[A, A], G[A, ]))
    sig <- function(k, mu)
        sign(-drop(Sig %*% k +G %*% mu))*I

    ## New directions
    s$s.d <- sig(s$k, s$mu)[s$A]
    s$mu.d <- -solve(GSG(s$A), drop(crossprod(G[s$A, ], solve(Sig[s$A, s$A], s$s.d))))
    s$k.d[s$A] <- solve(Sig[s$A, s$A], -drop(G[s$A, ] %*% s$mu.d +s$s.d))
    s$k.d[!s$A] <- 0
    s
}

#' Orthogonalized homotopy solution for l_infty
#' @param I vector of indicators which instruments are invalid
#' @keywords internal
linfh0 <- function(G, Sig, H, I) {
    dg <- nrow(G)
    res <- matrix(0, ncol=2*dg+1, nrow=1)
    colnames(res) <- c("lam", 1:dg, paste0("A", 1:dg))

    GSG <- function(A) crossprod(G[A, ], solve(Sig[A, A], G[A, ]))

    ## Initialize
    s <- list()
    s$lam <- 0
    s$A <- rep(TRUE, dg)
    s$mu <- solve(GSG(s$A), H)
    s$k <- drop(-solve(Sig, G %*% s$mu))
    s$joined <- 0
    res[1, -1] <- c(s$k, s$A)
    ## directions
    s$s.d <- sign(s$k) * I
    s$mu.d <- -solve(GSG(s$A), drop(crossprod(G, solve(Sig, s$s.d))))
    s$k.d <- solve(Sig, -drop(G %*% s$mu.d +s$s.d))
    j <- 1
    rr <- list(s)
    while (sum(s$A) > max(sum(!I), ncol(G))) {
        s <- linfstep(rr[[j]], G, Sig, H, I)
        rr[[j+1]] <- s
        j <- j+1
        res <- rbind(res, c(s$lam, s$k, s$A))
    }
    res
}


## linfbrute <- function(G, Sig, H, lams, I) {
##     ks <- matrix(NA, nrow=length(lams), ncol=nrow(G))
##     k <- CVXR::Variable(nrow(G))
##     for (j in seq_len(length(lams))) {
##         ob <- CVXR::Minimize(CVXR::p_norm(chol(Sig)%*%k)^2/2 +
##                              lams[j]*CVXR::p_norm(k[I], p=1))
##         pr <- CVXR::Problem(ob, list(-H==t(G)%*%k))
##         ks[j, ] <- solve(pr)$getValue(k)
##     }
##     cbind(lams, ks)
## }

## ## L_Infty penalized solution path for optimal sensitivities
## set.seed(42)
## dg <- 8
## nv <- 2
## I <- c(rep(TRUE, dg-nv), rep(FALSE, nv))
## dt <- 3
## G <- matrix(rnorm(dg*dt), nrow=dg)
## Sig <- cov(matrix(rnorm(2*dg*dg), ncol=dg))
## h <- c(1, rep(0, dt-1))


## res <- linfh(G, Sig, h, I)[, 1:(dg+1)]
## lams <- res[, 1]
## resb <- linfbrute(G, Sig, h, lams, I)
## diff <- max(abs(res[, 1:(dg+1)]-resb))
## if(diff>1e-3) message("Big difference\n") else print(diff)

## df <- data.frame(lambda=rep(lams, dg), k=as.vector(res[, 2:(dg+1)]),
##                  what=rep(as.factor(1:dg), each=length(lams)))
## pl1 <- qplot(x=lambda, y=k, color=what, geom="line", data=df)
## pdf("testlassoinf.pdf", width=5, height=4)
## print(directlabels::direct.label(pl1+theme_bw(), "first.qp"))
## dev.off()
