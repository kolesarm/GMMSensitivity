#' Homotopy solution
linfh <- function(eo, B, M=diag(ncol(B))) {
    ## Get orthogonalized homotopy

    ## B_{\perp}
    Bp <- if (nrow(B)>ncol(B)) {
              qr.Q(qr(B), complete=TRUE)[, (ncol(B)+1):nrow(B)]
          } else {
              matrix(ncol=0, nrow=nrow(B))
          }
    I <- rep(c(FALSE, TRUE), c(ncol(Bp), ncol(B)))

    T <- rbind(t(Bp), M %*% solve(crossprod(B), t(B)))

    Sigt <- T %*% eo$Sig %*% t(T)
    r <- linfh0(eo$G, Sigt, eo$H, I)
    kts <- r[, 2:(nrow(B)+1)]
    ## Return sensitivities at each step
    kts %*% T
}

#' Orthogonalized homotopy solution
#' @param I vector of indicators which instruments are invalid
linfh0 <- function(G, Sig, H, I) {
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
    mu <- solve(GSG(A), H)
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
            joined <- which(d2<=min(d2))
            A[joined] <- TRUE
        } else {
            A[which(d1<=min(d1))] <- FALSE
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

linfbrute <- function(G, Sig, H, lams, I) {
    ks <- matrix(NA, nrow=length(lams), ncol=nrow(G))
    k <- CVXR::Variable(nrow(G))
    for (j in seq_len(length(lams))) {
        ob <- CVXR::Minimize(CVXR::p_norm(chol(Sig)%*%k)^2/2 +
                             lams[j]*CVXR::p_norm(k[I], p=1))
        pr <- CVXR::Problem(ob, list(-H==t(G)%*%k))
        ks[j, ] <- solve(pr)$getValue(k)
    }
    cbind(lams, ks)
}

#' Build Estimators along a solution path
#' @param res output of \code{ATTh}
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @param C smoothness constant
#' @param sigma2 vector of variances
#' @param alpha CI coverage
#' @param beta quantile of excess length
linfPath <- function(res, d, y, C=1, sigma2,
                            alpha=0.05, beta=0.8) {
    n <- length(y)
    n0 <- n-sum(d)
    if (length(sigma2)==1) sigma2 <- sigma2*rep(1, n)

    if (length(dim(res)) > 1L) {
        m <- res[, 2:(n0+1)]
        r <- res[, (n0+2):(n+1)]
        mu <- res[, "mu"]
        delta <- res[, "delta"]
        att <- mean(y[d==1]) - drop(m %*% y[d==0]) / rowSums(m)
        maxbias <- rowMeans(r) - apply(m, 1, function(x) sum(x^2)) / rowSums(m)
        sd <- sqrt(mean(sigma2[d==1])/sum(d) + apply(m, 1, function(x)
            sum(x^2*sigma2[d==0])) / rowSums(m)^2)
        omega <- 2*(mu+rowMeans(r))
    } else {
        ## It's a vector
        m <- res[2:(n0+1)]
        r <- res[(n0+2):(n+1)]
        mu <- res["mu"]
        delta <- res["delta"]
        att <- mean(y[d==1]) - sum(m * y[d==0]) / sum(m)
        maxbias <- mean(r) -  sum(m^2) / sum(m)
        sd <- sqrt(mean(sigma2[d==1])/sum(d) + sum(m^2*sigma2[d==0]) / sum(m)^2)
        omega <- 2*(mu+mean(r))
    }
    maxbias <- C*maxbias
    hl <- CVb(maxbias/sd, alpha)$cv * sd
    lower <- att - maxbias - stats::qnorm(1-alpha)*sd
    upper <- att + maxbias + stats::qnorm(1-alpha)*sd
    ## worst-case quantile of excess length
    maxel <- 2*maxbias + sd * (stats::qnorm(1-alpha)+stats::qnorm(beta))

    data.frame(att=att, maxbias=maxbias, sd=sd, lower=lower, upper=upper, hl=hl,
               rmse=sqrt(sd^2+maxbias^2), maxel=maxel, omega=unname(omega),
               delta=unname(delta))
}

#' build optimal estimator given a solution path
#' @param res output of \code{ATTh}
#' @param ep output of \code{ATTEstimatePath}
#' @param y outcome vector
#' @param d vector of treatment indicators
#' @param C smoothness constant
#' @param sigma2 vector of variances
#' @param sigma2final vector of variances for final variance.
#' @param alpha CI coverage
#' @param beta quantile of excess length
#' @param opt.criterion One of \code{"RMSE"}, \code{"OCI"},  \code{"FLCI"}
#' @param UpdateC Update C that's assumed in \code{ep}?
#' @export
ATTOptEstimate <- function(res, ep, d, y, C=1, sigma2,
                           sigma2final=sigma2, alpha=0.05, beta=0.8,
                           opt.criterion="RMSE", UpdateC=TRUE) {
    ## Update estimate path with new value of C
    ep$maxbias <- C*ep$maxbias
    ep[, c("rmse", "lower", "upper", "hl", "maxel")] <-
        cbind(sqrt(ep$sd^2+ep$maxbias^2),
              ep$att - ep$maxbias - stats::qnorm(1-alpha)*ep$sd,
              ep$att + ep$maxbias + stats::qnorm(1-alpha)*ep$sd,
              CVb(ep$maxbias/ep$sd, alpha)$cv * ep$sd,
              2*ep$maxbias + ep$sd * (stats::qnorm(1-alpha)+stats::qnorm(beta)))

    ## Index of criterion to optimize
    idx <- if (opt.criterion=="RMSE") {
               which.max(names(ep)=="rmse")
           } else if (opt.criterion=="OCI") {
               which.max(names(ep)=="maxel")
           } else if (opt.criterion=="FLCI") {
               which.max(names(ep)=="hl")
           }
    i <- which.min(ep[[idx]])
    ## Assume minimum is either in [i, i+1] or [i-1, i], which is true if
    ## criterion is convex. Sometimes solution can be stuck, so not enough to
    ## take i+1 (but since which.min takes first minimum, we can always take i-1)
    ip <- i+1
    while(ip<=nrow(res) & ep$delta[i]==ep$delta[ip])
        ip <- ip+1

    if (ip<=nrow(res)) {
        f1 <- function(w)
            ATTEstimatePath((1-w)*res[i, ]+w*res[i+1, ], d, y, C,
                            sigma2, alpha, beta)[[idx]]
        opt1 <- stats::optimize(f1, interval=c(0, 1))
    } else {
        opt1 <- list(minimum=0, objective=min(ep[[idx]]))
    }

    if (i>1) {
        f0 <- function(w)
            ATTEstimatePath((1-w)*res[i-1, ]+w*res[i, ], d, y, C,
                            sigma2, alpha, beta)[[idx]]
        opt0 <- stats::optimize(f0, interval=c(0, 1))
    } else {
        opt0 <- list(minimum=1, objective=min(ep[[idx]]))
    }

    if (opt1$objective < opt0$objective) {
        resopt <- (1-opt1$minimum)*res[i, ] +
            opt1$minimum*res[min(i+1, nrow(res)), ]
    } else {
        resopt <- (1-opt0$minimum)*res[max(i-1, 1), ]+opt0$minimum*res[i, ]
    }

    r1 <- ATTEstimatePath(resopt, d, y, C, sigma2, alpha, beta)
    r2 <- ATTEstimatePath(resopt, d, y, C, sigma2final, alpha, beta)
    if (r1$delta==max(ep$delta))
        warning("Optimum found at end of path")
    cbind(r1, data.frame(rsd=r2$sd, rlower=r2$lower, rupper=r2$upper,
         rhl=r2$hl, rrmse=r2$rmse, rmaxel=r2$maxel, C=C))
}




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
