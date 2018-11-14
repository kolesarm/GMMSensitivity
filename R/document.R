#' Estimates from Berry, Levinsohn, and Pakes (1995)
#'
#' This dataset contains estimates of the model in Berry, Levinsohn, and Pakes
#' (1995), as implemented by Andrews, Gentzkow, and Shapiro (2017), and it is
#' used to illustrate the confidence intervals implemented in this package.
#' @format A list, consisting 11 objects:
#'
#' \describe{
#'
#' \item{G}{Matrix with 31 rows and 17 columns, estimate of derivative of the
#'          moment condition evaluated at initial esitmate of theta from Berry,
#'          Levinsohn, and Pakes (1995), \eqn{\hat{\theta}}{thetahat}.}
#'
#' \item{H}{Vector of length 17, estimate of derivative of average markup
#'   \eqn{h(\theta)}{h(theta)} with respect to model parameters
#'   \eqn{\theta}{theta}, evaluated at \eqn{\hat{\theta}}{thetahat}.}
#'
#' \item{W}{Weight matrix used to obtain the estimate
#'          \eqn{\hat{\theta}}{thetahat}.}
#'
#' \item{g_init}{Average moment condition, evaluated at
#'               \eqn{\hat{\theta}}{thetahat}.}
#'
#' \item{h_init}{Estimate of the average markup, \eqn{h(\hat{\theta})}{h(thetahat)}.}
#'
#' \item{names}{Two lists, one for names of the moment conditions, and one for
#'              elements of \eqn{\theta}{theta}.}
#'
#' \item{ZZ}{Gram matrix \eqn{Z'Z} of the instruments, used to specify the
#' misspecification set \eqn{C}.}
#'
#' \item{Sig}{Estimate of variance of moment condition.}
#'
#' \item{sdZ}{Vector of standard deviations of the instruments.}
#'
#' \item{perturb}{scaling parameters to give intepretable meaning to violations
#' of supply-side conditions. See vignette \code{vignette("GMMSensitivity",
#' package="GMMSensitivity")} for details.}
#'
#' \item{n}{Sample size, number of model/years.}
#' }
#'
#' See Armstrong and Kolesár (2018) for a detailed description of these objects.
#'
#'
#' @source Replication files for Andrews, Gentzkow, and Shapiro (2017), available at
#' \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/LLARSN/2KFPRA&version=1.1}
#' @references Andrews, I., M. Gentzkow, and J. M. Shapiro (2017): Measuring the Sensitivity of Parameter Estimates to Sample Statistics, Quarterly Journal of Economics, 132, 1553–1592.
#'
#' Armstrong, T. B., and M. Kolesár (2018): Sensitivity Analysis Using Approximate Moment Condition Models, Unpublished manuscript
#'
#' Berry, S. T., J. Levinsohn, and A. Pakes (1995): Automobile Prices in Market Equilibrium, Econometrica, 63, 841–890.
"blp"
