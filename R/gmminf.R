## ## L_Infty penalized solution path for optimal sensitivities
## set.seed(42)
## dg <- 8
## nv <- 2
## I <- c(rep(TRUE, dg-nv), rep(FALSE, nv))
## dt <- 3
## G <- matrix(rnorm(dg*dt), nrow=dg)
## Sig <- cov(matrix(rnorm(2*dg*dg), ncol=dg))
## h <- c(1, rep(0, dt-1))

## Homotopy solution
linfh <- function(G, Sig, h, I) {
    dt <- ncol(G)
    dg <- nrow(G)
    res <- matrix(0, ncol=2*dg+1, nrow=1)
    colnames(res) <- c("lam", 1:dg, paste0("A", 1:dg))

    GSG <- function(A) crossprod(G[A, ], solve(Sig[A, A], G[A, ]))
    s <- function(k, mu)
        sign(-drop(Sig %*% k +G %*% mu))*I

    ## Initialize
    lam <- 0
    A <- rep(TRUE, dg)
    mu <- solve(GSG(A), h)
    k <- drop(-solve(Sig, G %*% mu))
    res[1, -1] <- c(k, A)
    joined <- 0

    ## directions
    s.d <- sign(k) * I
    mu.d <- -solve(GSG(A), drop(crossprod(G, solve(Sig, s.d))))
    k.d <- solve(Sig, -drop(G %*% mu.d +s.d))

    while (sum(A) > max(sum(!I), dt)) {
        d1 <- -k/k.d
        d1[!A | d1<0 | !I] <- Inf
        if(joined!=0)
            d1[joined] <- Inf

        d2 <- rep(Inf, dg)
        a.d <- drop(Sig %*% k.d+G %*% mu.d)
        a <- drop(Sig %*% k+G %*% mu)
        if (sum(a.d[!A]>1)>0) {
            d2[a.d>1] <- ((lam-a)/(a.d-1))[a.d>1]
        } else if (sum(a.d[!A]< -1)>0) {
            d2[a.d< -1] <- (-(lam+a)/(a.d+1))[a.d< -1]
        }
        d2[A] <- Inf

        d <- min(d2, d1)
        lam <- lam+d
        k <- k+d*k.d
        mu <- mu+d*mu.d

        if (min(d2)<min(d1)){
            joined <- which.min(d2)
            A[joined] <- TRUE
        } else {
            A[which.min(d1)] <- FALSE
            joined <- 0
        }
        res <- rbind(res, c(lam, k, A))
        ## New directions
        s.d <- s(k, mu)[A]
        mu.d <- -solve(GSG(A), drop(crossprod(G[A, ], solve(Sig[A, A], s.d))))
        k.d[A] <- solve(Sig[A, A], -drop(G[A, ] %*% mu.d +s.d))
        k.d[!A] <- 0
    }
    res
}

linfbrute <- function(G, Sig, h, lams, I) {
    ks <- matrix(NA, nrow=length(lams), ncol=nrow(G))
    k <- CVXR::Variable(nrow(G))
    for (j in seq_len(length(lams))) {
        ob <- CVXR::Minimize(CVXR::p_norm(chol(Sig)%*%k)^2/2 +
                             lams[j]*CVXR::p_norm(k[I], p=1))
        pr <- CVXR::Problem(ob, list(-h==t(G)%*%k))
        ks[j, ] <- solve(pr)$getValue(k)
    }
    cbind(lams, ks)
}

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
