#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

args <- commandArgs(TRUE)
if(length(args) != 4) {
  stop("incorrect arguments")
}

filepath <- args[1]
filename <- args[2]
outpath <- args[3]
N <- as.numeric(args[4])

cat('running ', N, ' iterations\n')

args <- readRDS(file.path(filepath, filename))


# CONSTANTS ---------------------------------------------------------------

alpha <- 5
sigma <- 1
nu <- 1.5
nugget <- 0.01

init_nu <- 1
init_alpha <- 10
init_sigma <- 3

hyper_nu <- 1
hyper_alpha <- 100
hyper_sigma <- 1


# FUNCTIONS ---------------------------------------------------------------


##### Priors #####


# v_alpha around 100 or so?
d_prior_alpha <- function(alpha, v_alpha) {
  dunif(alpha, 0, v_alpha, log = T)
}

# v_nu should be pretty small
d_prior_nu <- function(nu, v_nu) {
  if (nu <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v_nu) - nu^2/(2*v_nu)
}

# v_sigma should probably be pretty small
d_prior_sigma <- function(sigma, v_sigma) {
  if (sigma <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v_sigma) - sigma^2/(2*v_sigma)
}


##### Proposals #####


r_q_alpha <- function(alpha, v_alpha) {
  rnorm(1, mean = alpha, sd = sqrt(v_alpha))
}

d_q_alpha <- function(alpha, center, v_alpha) {
  dnorm(alpha, mean = center, sd = sqrt(v_alpha), log = TRUE)
}

r_q_nu <- function(nu, v_nu) {
  rnorm(1, mean = nu, sd = sqrt(v_nu))
}

d_q_nu <- function(nu, center, v_nu) {
  dnorm(nu, mean = center, sd = sqrt(v_nu), log = TRUE)
}

r_q_sigma <- function(sigma, v_sigma) {
  rnorm(1, mean = sigma, sd = sqrt(v_sigma))
}

d_q_sigma <- function(sigma, center, v_sigma) {
  dnorm(sigma, mean = center, sd = sqrt(v_sigma), log = TRUE)
}

normal_ll_exact <- function(x, h, nu, alpha, sigma, tau) {
  lik <- foreach (i = 1:ncol(x), .combine = '+', .packages = 'gpcovr') %dopar% {
    covmat <- matern_cor(h, nu, alpha, sigma)
    diag(covmat) <- sigma + tau
    cholmat <- chol(covmat)
    log_det <- sum(log(diag(cholmat)))
    quadform <- drop(crossprod(backsolve(cholmat, x[ ,i], transpose = TRUE)))
    -0.5 * (length(x) * log(2*pi) + 2 * log_det + quadform)
  }
  lik
}


##### MCMC #####

mcmc_true <- function(numit, H, Y, nugget, init_nu, init_alpha, init_sigma, v_nu, v_alpha, v_sigma, hyper_nu, hyper_alpha, hyper_sigma, r.opt = 0.23, r = 5) {
  # constants
  M <- length(Y)
  c0 <- 10; c1 <- 0.8; k <- 3
  # storage and initialization
  nu_samps <-    nu_props <- c(init_nu, rep(NA, numit-1))
  alpha_samps <- alpha_props <- c(init_alpha, rep(NA, numit-1))
  sigma_samps <- sigma_props <- c(init_sigma, rep(NA, numit-1))
  
  nu_accepts <- c(1, rep(0, numit-1))
  alpha_accepts <- c(1, rep(0, numit-1))
  sigma_accepts <- c(1, rep(0, numit-1))
  
  gamma1 <- gamma2 <- rep(0, numit)
  
  v_nu_vec <- c(v_nu, rep(0, numit-1))
  v_alpha_vec <- c(v_alpha, rep(0, numit-1))
  v_sigma_vec <- c(v_sigma, rep(0, numit-1))
  # initialization
  nu <- init_nu
  alpha <- init_alpha
  sigma <- init_sigma
  cat('i = 1\n')
  lik_last <- normal_ll_exact(Y, H, nu, alpha, sigma, nugget)
  
  # loop
  for (i in 2:numit) {
    cat('i =', i, '\n')
    # adaptive tuning
    gamma1 <- c0 / (i + k)^c1
    if(i < r+1) {
      nu_r_hat <- mean(nu_accepts[1:(i-1)])
      alpha_r_hat <- mean(alpha_accepts[1:(i-1)])
      sigma_r_hat <- mean(sigma_accepts[1:(i-1)])
    } else {
      nu_r_hat <- mean(nu_accepts[(i-r):(i-1)])
      alpha_r_hat <- mean(alpha_accepts[(i-r):(i-1)])
      sigma_r_hat <- mean(sigma_accepts[(i-r):(i-1)])
    }
    v_nu_vec[i] <- exp(log(v_nu_vec[i-1]) + gamma1*(nu_r_hat - r.opt))
    v_alpha_vec[i] <- exp(log(v_alpha_vec[i-1]) + gamma1*(alpha_r_hat - r.opt))
    v_sigma_vec[i] <- exp(log(v_sigma_vec[i-1]) + gamma1*(sigma_r_hat - r.opt))
    ### propose parameters one at a time
    ## nu
    # propose
    nu_star <- r_q_nu(nu, v_nu_vec[i])
    # cat('proposed nu:', nu_star, '\n')
    # calculate likelihood
    if (sigma_accepts[i-1] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (nu_star < 0) -Inf else normal_ll_exact(Y, H, nu_star, alpha, sigma, nugget)
    # calculate acceptance probability
    numprior <- d_prior_nu(nu_star, hyper_nu)
    denomprior <- d_prior_nu(nu, hyper_nu)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      nu <- nu_star
      nu_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ## alpha
    # propose
    alpha_star <- r_q_alpha(alpha, v_alpha_vec[i])
    # cat('proposed alpha:', alpha_star, '\n')
    # calculate likelihood
    if (nu_accepts[i] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (alpha_star < 0) -Inf else normal_ll_exact(Y, H, nu, alpha_star, sigma, nugget)
    # calculate acceptance probability
    numprior <- d_prior_alpha(alpha_star, hyper_alpha)
    denomprior <- d_prior_alpha(alpha, hyper_alpha)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      alpha <- alpha_star
      alpha_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ## sigma
    # propose
    sigma_star <- r_q_sigma(sigma, v_sigma_vec[i])
    # cat('proposed sigma:', sigma_star, '\n')
    # calculate likelihood
    if (alpha_accepts[i] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (sigma_star < 0) -Inf else normal_ll_exact(Y, H, nu, alpha, sigma_star, nugget)
    # calculate acceptance probability
    numprior <- d_prior_sigma(sigma_star, hyper_sigma)
    denomprior <- d_prior_sigma(sigma, hyper_sigma)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      sigma <- sigma_star
      sigma_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ### save current values
    nu_samps[i] <- nu
    alpha_samps[i] <- alpha
    sigma_samps[i] <- sigma
  }
  
  result <- cbind(nu = nu_samps,
                  alpha = alpha_samps,
                  sigma = sigma_samps)
  colnames(result) <- c('nu', 'alpha', 'sigma')
  
  v <- cbind(nu = v_nu_vec,
             alpha = v_alpha_vec,
             sigma = v_sigma_vec)
  colnames(v) <- c('nu', 'alpha', 'sigma')
  return(list(samps = result, tuning = v))
}

result_true <- mcmc_true(N, args$dist_obs, args$Y, nugget, init_nu, init_alpha, init_sigma, 0.1, 0.1, 0.1, 1, 100, 5, r.opt = 0.23, r = 5)
saveRDS(result_true, file.path(outpath, 'result_true.rds'))