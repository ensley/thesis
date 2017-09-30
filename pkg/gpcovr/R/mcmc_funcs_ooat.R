#' (Log) density of the prior distribution on each spline coefficient
#'
#' This is a normal prior on each of the coefficients. The means are second-order differences, see Lang & Brezger (2001).
#'
#' @param beta The spline coefficients
#' @param tau The prior variance parameter
#' @param i The position of the coefficient within \code{beta}
#'
#' @return The prior log density for \code{beta[i]}
#' @export
#'
#' @examples
#' beta <- c(-6.5188650, 1.39889900,-0.2601626,
#'           -1.2675330, 0.04916187, 1.9759470,
#'            1.3406900,-2.23430300,-2.1876150,
#'            1.1648310, 2.01463700, 0.9876818,
#'           -0.4144843,-1.35885900,-4.1775160)
#' tau <- 0.01
#' ooat_dprior_beta(beta, tau, 4)
ooat_dprior_beta <- function(beta, tau, i) {
  if (tau <= 0) return(-Inf)
  if (i <= 2) result <- stats::dnorm(beta[i], mean = 0, sd = sqrt(tau), log = T)
  else result <- stats::dnorm(beta[i], mean = 2*beta[i-1] - beta[i-2], sd = sqrt(tau), log = T)
  result
}


#' (Log) density of the prior distribution of the variance scaling constant
#'
#' The variance scalar needs to be estimated as well. It can be interpreted as a
#' penalty term. The currently implemented prior is a half-normal distribution
#' with variance hyperparameter v.
#'
#' @param tau The variance scalar
#' @param v The variance hyperparameter
#'
#' @return The prior log density for tau
#' @export
#'
#' @examples
#' ooat_dprior_tau(1, 0.1)
#' ooat_dprior_tau(0, 0.1) # returns -Inf
ooat_dprior_tau <- function(tau, v) {
  if (tau <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v) - tau^2/(2*v)
}

#' (Log) density of the prior joint distribution of the spline coefficients and
#' variance scalar
#'
#' The priors for beta and tau are independent, so their joint distribution is
#' simply their product.
#'
#' @param beta The spline coefficients
#' @param tau The variance scalar
#' @param v_tau The variance hyperparameter
#'
#' @return The sum of \code{\link{ooat_dprior_beta}} for all \code{i} and \code{\link{ooat_dprior_tau}},
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
#' ooat_dprior(beta, tau, v_tau)
ooat_dprior <- function(beta, tau, v_tau) {
  ooat_dprior_tau(tau, v_tau) + sum(sapply(2:length(beta), function(i) ooat_dprior_beta(beta, tau, i)))
}

#' The proposal distribution for each of the spline coefficients
#'
#' This is just a normal distribution with mean \code{center} and variance
#' \code{v_beta}.
#'
#' @param beta The current spline coefficient
#' @param center Mean of the proposal distribution
#' @param v_beta Variance of the proposal distribution
#'
#' @return For \code{ooat_rq_beta}: A new proposed spline coefficient
#' @export
#'
#' @examples
#' # proposal function
#' ooat_rq_beta(2.2, 0.1)
ooat_rq_beta <- function(center, v_beta) {
  stats::rnorm(1, mean = center, sd = sqrt(v_beta))
}

#' @rdname ooat_rq_beta
#'
#' @return For \code{ooat_dq_beta}: The log density of beta, with respect to a
#'   distribution with mean \code{center}.
#' @export
#'
#' @examples
#'
#'
#' # density function
#' ooat_dq_beta(2, 2.2, 0.1)
ooat_dq_beta <- function(beta, center, v_beta) {
  stats::dnorm(beta, mean = center, sd = sqrt(v_beta), log = TRUE)
}

#' The proposal distribution for the variance scalar
#'
#' This is just a normal distribution with mean \code{center} and variance
#' \code{v}.
#'
#' @param tau The current variance scalar
#' @param center The mean of the proposal distribution
#' @param v The variance of the proposal distribution
#'
#' @return For \code{ooat_rq_tau}: A new proposed variance scalar
#' @export
#'
#' @examples
#' # proposal function
#' ooat_rq_tau(0.2, 0.05)
ooat_rq_tau <- function(center, v) {
  stats::rnorm(1, center, sqrt(v))
}

#' @rdname ooat_rq_tau
#'
#' @return For \code{ooat_dq_tau}: The log density of tau, with respect to a
#'   distribution with mean \code{center}
#' @export
#'
#' @examples
#'
#'
#' # density function
#' ooat_dq_tau(0.3, 0.2, 0.05)
ooat_dq_tau <- function(tau, center, v) {
  stats::dnorm(tau, mean = center, sd = sqrt(v), log = TRUE)
}
