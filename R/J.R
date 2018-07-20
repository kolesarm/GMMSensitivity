tol <- .Machine$double.eps^0.75

#' Find interval containing zero of a function, then find the zero
#'
#' Same function as in NPRHonest
#'
#' Given function \code{f} find \code{x0} such that \code{f(x0)==0}
#' @param f function whose root we're looking for
#' @param ival upper endpoint of initial interval in which to search
#' @param negative logical: should the lower endpoint be \code{1/ival} (if the
#'     root is guaranteed to be positive), or \code{-ival}?
#' @keywords internal
FindZero <- function(f, ival=1.1, negative=TRUE) {
    minval <- function(ival) if (negative==TRUE) -ival else min(1/ival, 1e-3)

    while(sign(f(ival))==sign(f(minval(ival))))
            ival <- 2*ival
    stats::uniroot(f, c(minval(ival), ival), tol=tol)$root
}

#' J-test of overidentifying restrictions under local misspecification
#'
#' Computes J-test of overidentifying restrictions with critical value adjusted
#' to allow for local misspecification, when the set \eqn{C} takes the form
#' \eqn{B*gamma} where the ell_p norm of M*gamma is bounded by K. Assumes
#' initial estimator in \code{eo} is optimal under correct specification.
#' @inheritParams
Jtest <- function(eo, B, M=diag(ncol(B)), K, p=2, alpha=0.05) {
    J <- eo$n*drop(crossprod(eo$g_init, solve(eo$Sig, eo$g_init)))
    ## Sigma^{-1/2}
    e <- eigen(eo$Sig)
    Sig12 <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)
    SG <- Sig12 %*% eo$G
    R <- diag(nrow(eo$G))-SG %*% solve(crossprod(SG), t(SG))
    A <- (R %*% Sig12 %*% B %*% solve(M))

    kbar <- if(p==2)
                max(eigen(crossprod(A))$values)
            else if (p==1)
                max(colSums(A^2))
            else
                -optim(par=rep(1, ncol(A)), function(x)
                    -drop(crossprod(A%*%x))/max(x^2), method="BFGS")$value

    ## p-value under correct specification
    p0 <- 1-pchisq(q=J, df=nrow(eo$G)-ncol(eo$G), ncp=0)

    ## Value of K for which p-value is 0.05
    K0 <- if (p0<=alpha) {
              FindZero(function(K) 0.95-pchisq(q=J, df=nrow(eo$G)-ncol(eo$G),
                                               ncp=K^2*kbar),
                       ival=2, negative=FALSE)
          } else {
              0
          }

    ## p-value
    pC <- 1-pchisq(q=J, df=nrow(eo$G)-ncol(eo$G), ncp=K^2*kbar)


    list(p0=p0, pC=pC, Kmin=K0)
}
