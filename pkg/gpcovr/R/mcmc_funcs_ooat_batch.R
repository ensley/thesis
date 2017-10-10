#' Estimate spline coefficients (initializer)
#'
#' This and \code{\link{estbeta_ooat}} probably should not be two separate functions
#' but oh well. This function contains the main MCMC loop.
#'
#' @param filename Name of the file to store the betas in
#' @param errname Name of the file to store the covariance functions in
#' @param nbatch Number of iterations to run this function for. See "Details"
#' @param H The distance matrix between the observations. \code{M} x \code{M}
#'   matrix, where \code{M} is the number of observations.
#' @param Y The vector of observations. Length \code{M}
#' @param B The number of random draws to take from the spline-approximated
#'   spectral density
#' @param init_beta Initial values for the spline coefficients
#' @param init_tau Initial value for the variance scalar
#' @param knots The vector of knot locations
#' @param v_beta The \code{beta} proposal variance hyperparameter
#' @param v_tau The \code{tau} proposal variance hyperparameter
#' @param t_v The prior variance hyperparameter (not being adapted)
#' @param r.opt The optimal acceptance rate (used for adaptive tuning)
#' @param r The number of previous iterations to use for adaptive tuning
#' @param eps Not used I think
#' @param nugget Nugget effect
#' @param errthin Thinning parameter for calculating error bounds. E.g. if
#'   \code{errthin} = 100, every 100th sample will be written out to the error
#'   file
#'
#' @details If \code{nbatch} is too big, the program will crash eventually
#'   because R is dumb and won't release the memory after each expensive
#'   likelihood calculation. A value of around 100 is recommended.
#'
#' @return Returns a list of information to pass to \code{\link{estbeta_ooat}}. Also writes output to some csv files.
#' @export
#'
estbeta_ooat_initialize <- function(filename, errname, nbatch, H, Y, B, init_beta, init_tau, knots, v_beta, v_tau, t_v, r.opt = 0.23, r = 5, eps = 1e-3, nugget = 1e-3, errthin = 100) {
  if (file.exists(filename)) file.remove(filename)
  if (file.exists(errname)) file.remove(errname)
  filecon <- file(filename, 'a+')
  errcon <- file(errname, 'a+')
  # define constants
  M <- length(Y)
  K <- length(knots)
  c0 <- 10; c1 <- 0.8; k <- 3
  # input checking
  if(nrow(H) != M || ncol(H) != M)  stop('H must be square with dimensions equal to the length of Y')
  if(length(init_beta) != K)        stop('length of initial beta guess must equal the number of knots')
  # storage
  beta_samps <- matrix(NA, nrow = nbatch, ncol = K)
  tau_samps <- rep(NA, length = nbatch)

  beta_samps[1, ] <- beta <- init_beta
  tau_samps[1] <- tau <- init_tau
  cat('i = 1\n')
  prev_lik <- lik_star <- likelihood(beta, knots, Y, H, B, nugget, debug = FALSE)
  last_accept <- TRUE

  curve <- gpumc::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))

  write(beta_samps[1, , drop = FALSE], filecon, ncolumns = K, append = T, sep = ',')
  write(curve, errcon, ncolumns = 500, append = T, sep = ',')

  for(i in 2:nbatch) {
    cat('\ni =', i, '\n')
    if(i >= 2) {
      # adaptive tuning
      gamma1 <- c0 / (i + k)^c1
      if(i == 2) {
        r.hat <- rep(1, K)
        r.hat.tau <- 1
      } else if(i < r+1) {
        r.hat <- apply(beta_samps[1:(i-1), ], 2, function(col) mean(diff(col) != 0))
        r.hat.tau <- mean(diff(tau_samps[1:(i-1)]) != 0)
      } else {
        r.hat <- apply(beta_samps[(i-r):(i-1), ], 2, function(col) mean(diff(col) != 0))
        r.hat.tau <- mean(diff(tau_samps[(i-r):(i-1)]) != 0)
      }
      cat('r.hat =', r.hat, '\n')
      # v_beta <- exp(log(v_beta) + gamma1*(r.hat - r.opt))
      v_beta <- rep(1e-4, K)
      v_tau <- exp(log(v_tau) + gamma1*(r.hat.tau - r.opt))
    }

    # propose betas
    for (j in 2:K) {
      cat('\tKNOT', j, '\n')
      cat('variance constant =', v_beta[j], '\n')
      beta_star <- beta
      beta_star[j] <- ooat_rq_beta(beta[j], v_beta[j])
      cat('new beta is', beta_star[j] - beta[j], 'larger than old beta\n')
      lik <- ifelse(last_accept, lik_star, prev_lik)
      cat('last draw was', ifelse(last_accept, 'ACCEPTED', 'REJECTED'), '\n')
      cat('last likelihood =', prev_lik, '\n')
      lik_star <- likelihood(beta_star, knots, Y, H, B, nugget, debug = FALSE)
      cat('proposed likelihood =', lik_star, '\n')
      # calculate acceptance probability
      numprior <- ooat_dprior_beta(beta_star, tau, j)
      denomprior <- ooat_dprior_beta(beta, tau, j)
      num <- lik_star + numprior
      denom <- lik + denomprior
      logr <- num - denom
      # accept or reject
      if(!is.nan(logr) & log(stats::runif(1)) < logr) {
        beta <- beta_star
        prev_lik <- lik_star
        last_accept <- TRUE
        cat('ACCEPTED\n')
      } else {
        prev_lik <- lik
        last_accept <- FALSE
        cat('REJECTED\n')
      }
    }


    # propose tau
    tau_star <- ooat_rq_tau(tau, v_tau)
    # calculate acceptance probability
    num <- ooat_dprior(beta, tau_star, t_v)
    denom <- ooat_dprior(beta, tau, t_v)
    logr <- num - denom
    # accept or reject
    if(log(stats::runif(1)) < logr) {
      tau <- tau_star
    }

    # save current values
    beta_samps[i, ] <- beta
    tau_samps[i] <- tau

    # write information out to files
    write(beta_samps[i, , drop = FALSE], file = filecon, ncolumns = K, append = T, sep = ',')

    if (i %% errthin == 0) {
      # get covariance function error bounds for this set of betas and write it out
      curve <- gpumc::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))
      write(curve, errcon, ncolumns = 500, append = T, sep = ',')
    }
  }

  flush(filecon)
  flush(errcon)
  close(filecon)
  close(errcon)

  args <- list(filename = filename,
               errname = errname,
               H = H,
               Y = Y,
               B = B,
               knots = knots,
               r.opt = r.opt,
               r = r,
               eps = eps,
               t_v = t_v,
               nugget = nugget)

  output <- list(args = args,
                 last_few_beta = utils::tail(beta_samps, r),
                 last_few_tau = utils::tail(tau_samps, r),
                 v_beta = v_beta,
                 v_tau = v_tau,
                 prev_lik = prev_lik,
                 last_accept = last_accept)

  return(output)
}


#' Estimate spline coefficients
#'
#' This and \code{\link{estbeta_ooat_initialize}} probably should not be two separate
#' functions but oh well. This function contains the main MCMC loop.
#'
#' @param filename Name of the file to store the betas in
#' @param errname Name of the file to store the covariance functions in
#' @param nbatch Number of iterations to run this function for. See "Details"
#' @param batchidx Index of the current batch
#' @param H The distance matrix between the observations. \code{M} x \code{M}
#'   matrix, where \code{M} is the number of observations.
#' @param Y The vector of observations. Length \code{M}
#' @param B The number of random draws to take from the spline-approximated
#'   spectral density
#' @param last_few_beta Last \code{r} values for the spline coefficients
#' @param last_few_tau Last \code{r} for the variance scalar
#' @param prev_lik The most recent likelihood value
#' @param last_accept Whether \code{prev_lik} was accepted or not
#' @param knots The vector of knot locations
#' @param v_beta The \code{beta} proposal variance hyperparameter
#' @param v_tau The \code{tau} proposal variance hyperparameter
#' @param t_v The prior variance hyperparameter (not being adapted)
#' @param r.opt The optimal acceptance rate (used for adaptive tuning)
#' @param eps Not used I think
#' @param nugget Nugget effect
#' @param errthin Thinning parameter for calculating error bounds. E.g. if
#'   \code{errthin} = 100, every 100th sample will be written out to the error
#'   file
#'
#' @details If \code{nbatch} is too big, the program will crash eventually
#'   because R is dumb and won't release the memory after each expensive
#'   likelihood calculation. A value of around 100 is recommended.
#'
#' @return Returns a list of information to pass to \code{\link{estbeta_ooat}}. Also
#'   writes output to some csv files.
#' @export
#'
estbeta_ooat <- function(filename, errname, nbatch, batchidx, H, Y, B, last_few_beta, last_few_tau, prev_lik, last_accept, knots, v_beta, v_tau, t_v, r.opt = 0.23, eps, nugget = 1e-3, errthin = 100) {
  filecon <- file(filename, 'a+')
  errcon <- file(errname, 'a+')
  # define constants
  M <- length(Y)
  K <- length(knots)
  r <- length(last_few_tau)
  c0 <- 10; c1 <- 0.8; k <- 3
  # input checking
  if(nrow(H) != M | ncol(H) != M) stop('H must be square with dimensions equal to the length of Y')
  if(ncol(last_few_beta) != K) stop('length of initial beta guess must equal the number of knots')
  # storage
  beta_samps <- rbind(last_few_beta, matrix(NA, nrow = nbatch, ncol = K))
  tau_samps <- c(last_few_tau, rep(NA, length = nbatch))

  beta <- last_few_beta[nrow(last_few_beta), ]
  tau <- last_few_tau[length(last_few_tau)]
  lik_star <- prev_lik

  for(i in 1:nbatch) {
    ii <- nbatch * batchidx + i
    cat('\ni =', ii, '\n')
    # adaptive tuning
    # scaling factor
    gamma1 <- c0 / (nbatch * batchidx + i + k)^c1
    r.hat <- apply(beta_samps[i:(i+r-1), ], 2, function(col) mean(diff(col) != 0))
    cat('r.hat =', r.hat, '\n')
    r.hat.tau <- mean(diff(tau_samps[i:(i+r-1)]) != 0)
    # v_beta <- exp(log(v_beta) + gamma1*(r.hat - r.opt))
    v_beta <- rep(1e-4, K)
    v_tau <- exp(log(v_tau) + gamma1*(r.hat.tau - r.opt))

    # propose betas
    for(j in 2:K) {
      cat('\tKNOT', j, '\n')
      cat('proposal variance =', v_beta[j], '\n')
      beta_star <- beta
      beta_star[j] <- ooat_rq_beta(beta[j], v_beta[j])
      cat('new beta is', beta_star[j] - beta[j], 'larger than old beta\n')
      if(last_accept) {
        lik <- lik_star
      } else {
        lik <- prev_lik
      }
      cat('last draw was', ifelse(last_accept, 'ACCEPTED', 'REJECTED'), '\n')
      cat('last likelihood =', prev_lik, '\n')
      lik_star <- likelihood(beta_star, knots, Y, H, B, nugget, debug = FALSE)
      cat('proposed likelihood =', lik_star, '\n')
      # calculate acceptance probability
      numprior <- ooat_dprior_beta(beta_star, tau, j)
      denomprior <- ooat_dprior_beta(beta, tau, j)
      num <- lik_star + numprior
      denom <- lik + denomprior
      logr <- num - denom
      # accept or reject
      if(!is.nan(logr) & log(stats::runif(1)) < logr) {
        beta <- beta_star
        prev_lik <- lik_star
        last_accept <- TRUE
        cat('ACCEPTED\n')
      } else {
        prev_lik <- lik
        last_accept <- FALSE
        cat('REJECTED\n')
      }
    }

    # propose tau
    tau_star <- ooat_rq_tau(tau, v_tau)
    # calculate acceptance probability
    num <- ooat_dprior(beta, tau_star, t_v)
    denom <- ooat_dprior(beta, tau, t_v)
    logr <- num - denom
    # accept or reject
    if(log(stats::runif(1)) < logr) {
      tau <- tau_star
    }

    # save current values
    beta_samps[i+r, ] <- beta
    tau_samps[i+r] <- tau

    # write information out to files
    write(beta_samps[i+r, , drop = FALSE], file = filecon, ncolumns = K, append = T, sep = ',')

    if (i %% errthin == 0) {
      # get covariance function error bounds for this set of betas
      curve <- gpumc::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))
      write(curve, errcon, ncolumns = 500, append = T, sep = ',')
    }
  }

  flush(filecon)
  flush(errcon)
  close(filecon)
  close(errcon)

  args <- list(filename = filename,
               errname = errname,
               H = H,
               Y = Y,
               B = B,
               knots = knots,
               r.opt = r.opt,
               r = r,
               eps = eps,
               t_v = t_v,
               nugget = nugget)

  output <- list(args = args,
                 last_few_beta = utils::tail(beta_samps, r),
                 last_few_tau = utils::tail(tau_samps, r),
                 v_beta = v_beta,
                 v_tau = v_tau,
                 prev_lik = prev_lik,
                 last_accept = last_accept)

  return(output)
}
