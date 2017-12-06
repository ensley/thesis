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

args <- readRDS(file.path(filepath, filename))


# CONSTANTS ---------------------------------------------------------------

lambda <- 5
var <- 1
scale <- 1.5
nugget <- 0.01

init_lambda <- 2
init_var <- 1
init_scale <- 1

hyper_lambda <- 100
hyper_var <- 1
hyper_scale <- 1


# FUNCTIONS ---------------------------------------------------------------

damped_cov <- function(h, lambda, var, scale) {
  hh <- h / scale
  var * exp(-lambda * hh) * cos(hh)
}

##### Priors #####


# v_lambda around 100 or so?
d_prior_lambda <- function(lambda, v_lambda) {
  if (lambda <= 1) return(-Inf)
  dunif(lambda, 1, v_lambda, log = T)
}

# v_var should be pretty small
d_prior_var <- function(var, v_var) {
  if (var <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v_var) - var^2/(2*v_var)
}

# v_scale should probably be pretty small
d_prior_scale <- function(scale, v_scale) {
  if (scale <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v_scale) - scale^2/(2*v_scale)
}


##### Proposals #####


r_q_lambda <- function(lambda, v_lambda) {
  rnorm(1, mean = lambda, sd = sqrt(v_lambda))
}

d_q_lambda <- function(lambda, center, v_lambda) {
  dnorm(lambda, mean = center, sd = sqrt(v_lambda), log = TRUE)
}

r_q_var <- function(var, v_var) {
  rnorm(1, mean = var, sd = sqrt(v_var))
}

d_q_var <- function(var, center, v_var) {
  dnorm(var, mean = center, sd = sqrt(v_var), log = TRUE)
}

r_q_scale <- function(scale, v_scale) {
  rnorm(1, mean = scale, sd = sqrt(v_scale))
}

d_q_scale <- function(scale, center, v_scale) {
  dnorm(scale, mean = center, sd = sqrt(v_scale), log = TRUE)
}

normal_ll_exact <- function(x, h, lambda, var, scale, tau) {
  cat('\tlambda =', lambda, '\n\tvar =', var, '\n\tscale =', scale, '\n')
  lik <- foreach (i = 1:ncol(x), .combine = '+', .export = 'damped_cov') %dopar% {
    covmat <- damped_cov(h, lambda, var, scale)
    diag(covmat) <- var + tau
    cholmat <- chol(covmat)
    log_det <- sum(log(diag(cholmat)))
    quadform <- drop(crossprod(backsolve(cholmat, x[ ,i], transpose = TRUE)))
    -0.5 * (length(x) * log(2*pi) + 2 * log_det + quadform)
  }
  lik
}


##### MCMC #####

mcmc_true <- function(numit, H, Y, nugget, init_lambda, init_var, init_scale, v_lambda, v_var, v_scale, hyper_lambda, hyper_var, hyper_scale, r.opt = 0.23, r = 5) {
  # constants
  M <- length(Y)
  c0 <- 10; c1 <- 0.8; k <- 3
  # storage and initialization
  lambda_samps    <-    lambda_props <- c(init_lambda, rep(NA, numit-1))
  var_samps <- var_props <- c(init_var, rep(NA, numit-1))
  scale_samps <- scale_props <- c(init_scale, rep(NA, numit-1))
  
  lambda_accepts <- c(1, rep(0, numit-1))
  var_accepts <- c(1, rep(0, numit-1))
  scale_accepts <- c(1, rep(0, numit-1))
  
  gamma1 <- gamma2 <- rep(0, numit)
  
  v_lambda_vec <- c(v_lambda, rep(0, numit-1))
  v_var_vec <- c(v_var, rep(0, numit-1))
  v_scale_vec <- c(v_scale, rep(0, numit-1))
  # initialization
  lambda <- init_lambda
  var <- init_var
  scale <- init_scale
  cat('i = 1\n')
  lik_last <- normal_ll_exact(Y, H, lambda, var, scale, nugget)
  
  # loop
  for (i in 2:numit) {
    cat('i =', i, '\n')
    # adaptive tuning
    gamma1 <- c0 / (i + k)^c1
    if(i < r+1) {
      lambda_r_hat <- mean(lambda_accepts[1:(i-1)])
      var_r_hat <- mean(var_accepts[1:(i-1)])
      scale_r_hat <- mean(scale_accepts[1:(i-1)])
    } else {
      lambda_r_hat <- mean(lambda_accepts[(i-r):(i-1)])
      var_r_hat <- mean(var_accepts[(i-r):(i-1)])
      scale_r_hat <- mean(scale_accepts[(i-r):(i-1)])
    }
    v_lambda_vec[i] <- exp(log(v_lambda_vec[i-1]) + gamma1*(lambda_r_hat - r.opt))
    v_var_vec[i] <- exp(log(v_var_vec[i-1]) + gamma1*(var_r_hat - r.opt))
    v_scale_vec[i] <- exp(log(v_scale_vec[i-1]) + gamma1*(scale_r_hat - r.opt))
    ### propose parameters one at a time
    ## lambda
    # propose
    lambda_star <- r_q_lambda(lambda, v_lambda_vec[i])
    # cat('proposed lambda:', lambda_star, '\n')
    # calculate likelihood
    if (scale_accepts[i-1] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (lambda_star < 0) -Inf else normal_ll_exact(Y, H, lambda_star, var, scale, nugget)
    # calculate acceptance probability
    numprior <- d_prior_lambda(lambda_star, hyper_lambda)
    denomprior <- d_prior_lambda(lambda, hyper_lambda)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      lambda <- lambda_star
      lambda_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ## var
    # propose
    var_star <- r_q_var(var, v_var_vec[i])
    # cat('proposed var:', var_star, '\n')
    # calculate likelihood
    if (lambda_accepts[i] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (var_star < 0) -Inf else normal_ll_exact(Y, H, lambda, var_star, scale, nugget)
    # calculate acceptance probability
    numprior <- d_prior_var(var_star, hyper_var)
    denomprior <- d_prior_var(var, hyper_var)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      var <- var_star
      var_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ## scale
    # propose
    scale_star <- r_q_scale(scale, v_scale_vec[i])
    # cat('proposed scale:', scale_star, '\n')
    # calculate likelihood
    if (var_accepts[i] == 1) {
      # if we accepted the last parameter, the likelihood is equal to the previous
      # iteration's proposal likelihood
      likelihood <- lik_last
    }
    # otherwise, the likelihood is the same as the likelihood from the previous iteration.
    # now calculate the new proposal likelihood
    likelihood_star <- if (scale_star < 0) -Inf else normal_ll_exact(Y, H, lambda, var, scale_star, nugget)
    # calculate acceptance probability
    numprior <- d_prior_scale(scale_star, hyper_scale)
    denomprior <- d_prior_scale(scale, hyper_scale)
    num <- likelihood_star + numprior
    denom <- likelihood + denomprior
    logr <- num - denom
    # accept or reject
    if(!is.nan(logr) & log(runif(1)) < logr) {
      scale <- scale_star
      scale_accepts[i] <- 1
      lik_last <- likelihood_star
      # if(debug) cat('ACCEPTED\n')
    } else {
      lik_last <- likelihood
      # if(debug) cat('REJECTED\n')
    }
    ### save current values
    lambda_samps[i] <- lambda
    var_samps[i] <- var
    scale_samps[i] <- scale
  }
  
  result <- cbind(lambda = lambda_samps,
                  var = var_samps,
                  scale = scale_samps)
  colnames(result) <- c('lambda', 'var', 'scale')
  
  v <- cbind(lambda = v_lambda_vec,
             var = v_var_vec,
             scale = v_scale_vec)
  colnames(v) <- c('lambda', 'var', 'scale')
  return(list(samps = result, tuning = v))
}

result_true <- mcmc_true(N, args$dist_obs, args$Y, nugget, init_lambda, init_var, init_scale, 0.1, 0.1, 0.1, 100, 1, 5, r.opt = 0.23, r = 5)
saveRDS(result_true, file.path(outpath, 'result_true.rds'))