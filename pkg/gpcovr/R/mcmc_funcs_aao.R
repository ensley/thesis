#' Prior covariance matrix for the spline coefficients
#'
#' This is a penalty matrix of second-order differences, as described in Lang &
#' Brezger (2001).
#'
#' @param n The number of spline coefficients
#'
#' @return An n x n penalty matrix
buildmat <- function(n) {
  if(n < 3) stop('n must be greater than 3')
  headpart <- if(n == 3) c(1, -2, 1, -2) else c(1, -2, 1, 0)
  cvec <- 1:n
  for(i in seq_len(n-3)) {
    cvec <- c(cvec, rep(n-1, n), n)
  }
  cvec <- c(cvec, (n-1):1, rep(0, 4))
  matvec <- c(headpart, diff(diff(diff(diff(cvec)))))
  matrix(matvec, nrow = n, ncol = n)
}



#' (Log) density of the joint prior distribution on the spline coefficients
#'
#' This is a \code{k}-dimensional multivariate normal distribution, where
#' \code{k} is the dimension of \code{beta}. The covariance is currently set to
#' be the same as the proposal covariance (see \code{\link{aao_dq_beta}}).
#'
#' @param beta The spline coefficients
#' @param tau A constant that scales the covariance matrix (also will be
#'   estimated by the MCMC)
#'
#' @return The prior log density for beta
#' @export
#'
#' @examples
#' beta <- c(-6.5188650, 1.39889900,-0.2601626,
#'           -1.2675330, 0.04916187, 1.9759470,
#'            1.3406900,-2.23430300,-2.1876150,
#'            1.1648310, 2.01463700, 0.9876818,
#'           -0.4144843,-1.35885900,-4.1775160)
#' tau <- 0.01
#' aao_dprior_beta(beta, tau)
aao_dprior_beta <- function(beta, tau) {
  if (tau <= 0) return(-Inf)
  n <- length(beta)
  -1 * beta %*% buildmat(n) %*% beta / (2*tau)
}


#' (Log) density of the prior distribution of the variance scaling constant
#'
#' The covariance scalar needs to be estimated as well. It can be interpreted as a
#' penalty term. The currently implemented prior is a half-normal distribution
#' with variance hyperparameter v.
#'
#' @param tau The covariance scalar
#' @param v The variance hyperparameter
#'
#' @return The prior log density for tau
#' @export
#'
#' @examples
#' aao_dprior_tau(1, 0.1)
#' aao_dprior_tau(0, 0.1) # returns -Inf
aao_dprior_tau <- function(tau, v) {
  if (tau <= 0) return(-Inf)
  -0.5*log(v) - (tau * tau)/(2*v)
}



#' (Log) density of the prior joint distribution of the spline coefficients and
#' covariance scalar
#'
#' The priors for beta and tau are independent, so their joint distribution is
#' simply their product.
#'
#' @param beta The spline coefficients
#' @param tau The covariance scalar
#' @param v_tau The variance hyperparameter
#'
#' @return The sum of \code{\link{aao_dprior_beta}} and \code{\link{aao_dprior_tau}}
#'   with the appropriate arguments.
#' @export
#'
#' @examples
#' beta <- c(-6.5188650, 1.39889900,-0.2601626,
#'           -1.2675330, 0.04916187, 1.9759470,
#'            1.3406900,-2.23430300,-2.1876150,
#'            1.1648310, 2.01463700, 0.9876818,
#'           -0.4144843,-1.35885900,-4.1775160)
#' tau <- 0.1
#' v_tau <- 0.01
#' aao_dprior(beta, tau, v_tau)
aao_dprior <- function(beta, tau, v_tau) {
  aao_dprior_beta(beta, tau) + aao_dprior_tau(tau, v_tau)
}



#' The proposal distribution for the spline coefficients
#'
#' Currently this is implemented as a \code{k}-dimensional multivariate normal
#' distribution. The covariance matrix is that of a
#' penalty matrix of second-order differences, as described in Lang & Brezger
#' (2001).
#'
#' @param beta The spline coefficients
#' @param center The mean of the proposal distribution
#' @param tau The variance hyperparameter
#'
#' @return For \code{aao_rq_beta}: A new set of proposed spline coefficients
#' @export
#'
#' @examples
#' # proposal function
#' beta <- c(-6.5188650, 1.39889900,-0.2601626,
#'           -1.2675330, 0.04916187, 1.9759470,
#'            1.3406900,-2.23430300,-2.1876150,
#'            1.1648310, 2.01463700, 0.9876818,
#'           -0.4144843,-1.35885900,-4.1775160)
#' tau <- 0.1
#' aao_rq_beta(beta, tau)
aao_rq_beta <- function(center, tau) {
  n <- length(center)
  drop(mvtnorm::rmvnorm(1, mean = center, sigma = tau*buildmat(n)))
}

#' @rdname aao_rq_beta
#' @return For \code{aao_dq_beta}: The log density of beta, with respect to a
#'   distribution with mean \code{center}
#' @export
#'
#' @examples
#'
#'
#' # density function
#' beta <- c(-6.5188650, 1.39889900,-0.2601626,
#'           -1.2675330, 0.04916187, 1.9759470,
#'            1.3406900,-2.23430300,-2.1876150,
#'            1.1648310, 2.01463700, 0.9876818,
#'           -0.4144843,-1.35885900,-4.1775160)
#' tau <- 0.1
#' beta2 <- aao_rq_beta(beta, tau)
#' aao_dq_beta(beta, beta2, tau)
aao_dq_beta <- function(beta, center, tau) {
  n <- length(beta)
  mvtnorm::dmvnorm(beta, mean = center, sigma = tau*buildmat(n), log = TRUE)
}



#' The proposal distribution for the covariance scalar
#'
#' @param tau The covariance scalar
#' @param center The mean of the proposal distribution
#' @param v The variance hyperparameter
#'
#' @return For \code{aao_rq_tau}: A new proposed covariance scalar
#' @export
#'
#' @examples
#' # proposal function
#' aao_rq_tau(0.5, 0.1)
aao_rq_tau <- function(center, v) {
  stats::rnorm(1, center, sqrt(v))
}

#' @rdname aao_rq_tau
#' @return For \code{aao_dq_tau}: The log density of tau, with respect to a
#'   distribution with mean \code{center}
#' @export
#'
#' @examples
#'
#'
#' # density function
#' aao_dq_tau(0.5, 0.4, 0.1)
aao_dq_tau <- function(tau, center, v) {
  stats::dnorm(tau, mean = center, sd = sqrt(v), log = TRUE)
}
