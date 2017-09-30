#! /usr/bin --vanilla --default-packages=utils

library(cudatest)
library(openacc)
library(coda)
library(splines)
library(mvtnorm)
library(methods)
library(RColorBrewer)


d_prior_beta <- function(beta, tau, i) {
  if (tau <= 0) return(-Inf)
  if (i <= 2) result <- dnorm(beta[i], mean = 0, sd = sqrt(tau), log = T)
  else result <- dnorm(beta[i], mean = 2*beta[i-1] - beta[i-2], sd = sqrt(tau), log = T)
  result
}


d_prior_tau <- function(tau, v) {
  if (tau <= 0) return(-Inf)
  # half-normal prior
  -0.5*log(v) - tau^2/(2*v)
}

d_prior <- function(beta, tau, v_tau) {
  d_prior_tau(tau, v_tau) + sum(sapply(2:length(beta), function(i) d_prior_beta(beta, tau, i)))
}

r_q_beta <- function(beta, v_beta) {
  rnorm(1, mean = beta, sd = sqrt(v_beta))
}

d_q_beta <- function(beta, center, v_beta) {
  dnorm(beta, mean = center, sd = sqrt(v_beta), log = TRUE)
}

r_q_tau <- function(tau, v) {
  rnorm(1, tau, sqrt(v))
}

d_q_tau <- function(tau_star, tau, v) {
  dnorm(tau_star, mean = tau, sd = sqrt(v), log = TRUE)
}



likelihood <- function(beta, knots, Y, H, B, nugget, debug) {
  M <- length(Y)

  slopes <- get_slopes(beta, knots)

  if(debug) cat('slopes:', slopes, '\n')
  if(slopes[2] >= 0 | slopes[1] <= 0) return(-Inf)
  mineigen <- -Inf
  cnt <- 0
  while(mineigen < 0) {
    cnt <- cnt + 1
    if(cnt != 1) cat('Redrawing...\n')
    if(cnt > 2) {
      cat('Likelihood calculation failed. Returning -Inf\n')
      return(-Inf)
    }
    # take random samples from spectral density
    wtilde <- rlogspline(B, beta, knots)
    if(any(is.infinite(exp(wtilde)))) return(-Inf)
    if(debug) {
      ### plot current density estimate
      # dev.off()
      par(mfrow = c(2,1))
      x <- seq(-20, 5, length = 200)
      plot(x, predict_natspl(nsbasis(x, knots), beta), type = 'l', ylim = c(-15, 5))
      lines(xx, yy, col = 'orange')
      abline(v = knots, col = 'grey', lty = 2)
      hist(wtilde, probability = T, breaks = 50); lines(x, dlogspline(x, beta, knots))
    }
    # estimate covariance matrix via MC integration
    Sigma <- diag(1 + nugget, nrow = M, ncol = M)
    Sigma[lower.tri(Sigma)] <- cudatest::mc(H[lower.tri(H)], exp(wtilde))
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
    mineigen <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$value)
    if(debug) cat('min eigval is', mineigen, '\n')
    rm(wtilde)
  }
  cholmat <- chol(Sigma)
  logdet <- sum(log(diag(cholmat)))
  quadform <- drop(crossprod(backsolve(cholmat, Y, transpose = TRUE)))
  l <- -0.5 * (M * log(2*pi) + 2*logdet + quadform)
  rm(Sigma, cholmat)
  gc()
  return(l)
}


estbeta_initialize <- function(filename, errname, nbatch, H, Y, B, init_beta, init_tau, knots, v_beta, v_tau, t_v, r.opt = 0.23, r = 5, eps = 1e-3, nugget = 1e-3) {
  filecon <- file(filename, 'a+')
  errcon <- file(errname, 'a+')
  # define constants
  M <- length(Y)
  K <- length(knots)
  c0 <- 10; c1 <- 0.8; k <- 3
  # input checking
  if(nrow(H) != M || ncol(H) != M) stop('H must be square with dimensions equal to the length of Y')
  if(length(init_beta) != K) stop('length of initial beta guess must equal the number of knots')
  # storage
  beta_samps <- matrix(NA, nrow = nbatch, ncol = K)
  tau_samps <- rep(NA, length = nbatch)

  beta_samps[1, ] <- beta <- init_beta
  tau_samps[1] <- tau <- init_tau
  cat('i = 1\n')
  prev_lik <- lik_star <- likelihood(beta, knots, Y, H, B, nugget, debug = FALSE)
  last_accept <- TRUE

  curve <- cudatest::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))
  
  write(beta_samps[1, , drop = FALSE], filecon, ncolumns = K, append = T, sep = ',')
  write(curve, errcon, ncolumns = 500, append = T, sep = ',')

  for(i in 2:nbatch) {
    cat('\ni =', i, '\n')
    if(i > 2) {
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
      beta_star[j] <- r_q_beta(beta[j], v_beta[j])
      cat('new beta is', beta_star[j] - beta[j], 'larger than old beta\n')
      lik <- ifelse(last_accept, lik_star, prev_lik)
      cat('last draw was', ifelse(last_accept, 'ACCEPTED', 'REJECTED'), '\n')
      cat('last likelihood =', prev_lik, '\n')
      lik_star <- likelihood(beta_star, knots, Y, H, B, nugget, debug = FALSE)
      cat('proposed likelihood =', lik_star, '\n')
      # calculate acceptance probability
      numprior <- d_prior_beta(beta_star, tau, j)
      denomprior <- d_prior_beta(beta, tau, j)
      num <- lik_star + numprior
      denom <- lik + denomprior
      logr <- num - denom
      # accept or reject
      if(!is.nan(logr) & log(runif(1)) < logr) {
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
    tau_star <- r_q_tau(tau, v_tau)
    # calculate acceptance probability
    num <- d_prior(beta, tau_star, t_v)
    denom <- d_prior(beta, tau, t_v)
    logr <- num - denom
    # accept or reject
    if(log(runif(1)) < logr) {
      tau <- tau_star
    }

    # save current values
    beta_samps[i, ] <- beta
    tau_samps[i] <- tau

    # get covariance function error bounds for this set of betas
    curve <- cudatest::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))
    
    # write information out to files
    write(beta_samps[i, , drop = FALSE], filecon, ncolumns = K, append = T, sep = ',')
    write(curve, errcon, ncolumns = 500, append = T, sep = ',')
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
                 last_few_beta = tail(beta_samps, r),
                 last_few_tau = tail(tau_samps, r),
                 v_beta = v_beta,
                 v_tau = v_tau,
                 prev_lik = prev_lik,
                 last_accept = last_accept)

  return(output)
}


estbeta <- function(filename, errname, nbatch, batchidx, H, Y, B, last_few_beta, last_few_tau, prev_lik, last_accept, knots, v_beta, v_tau, t_v, r.opt = 0.23, eps, nugget = 1e-3) {
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
      beta_star[j] <- r_q_beta(beta[j], v_beta[j])
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
      numprior <- d_prior_beta(beta_star, tau, j)
      denomprior <- d_prior_beta(beta, tau, j)
      num <- lik_star + numprior
      denom <- lik + denomprior
      logr <- num - denom
      # accept or reject
      if(!is.nan(logr) & log(runif(1)) < logr) {
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
    tau_star <- r_q_tau(tau, v_tau)
    # calculate acceptance probability
    num <- d_prior(beta, tau_star, t_v)
    denom <- d_prior(beta, tau, t_v)
    logr <- num - denom
    # accept or reject
    if(log(runif(1)) < logr) {
      tau <- tau_star
    }

    # save current values
    beta_samps[i+r, ] <- beta
    tau_samps[i+r] <- tau

    # get covariance function error bounds for this set of betas
    curve <- cudatest::mc(seq(1e-5, 5, length = 500), exp(rlogspline(B, beta, knots)))
    
    # write information out to files
    write(beta_samps[i+r, , drop = FALSE], ncolumns = K, filecon, append = T, sep = ',')
    write(curve, errcon, ncolumns = 500, append = T, sep = ',')
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
                 last_few_beta = tail(beta_samps, r),
                 last_few_tau = tail(tau_samps, r),
                 v_beta = v_beta,
                 v_tau = v_tau,
                 prev_lik = prev_lik,
                 last_accept = last_accept)

  return(output)
}


create_locations <- function(M, ngrid, mindist = 0.005) {
  # create grid of prediction observations
  xpred <- seq(0, 1, length = ngrid)
  grid_pred <- expand.grid(x = xpred, y = xpred)
  # create initial grid for non-gridded observations over unit square.
  # points 3*mindist apart.
  xobs <- seq(mindist, 1-mindist, by = 3*mindist)
  grid_obs <- expand.grid(x = xobs, y = xobs)
  # perturb grid points by no more than mindist in any direction
  grid_obs <- grid_obs + runif(2 * length(xobs) * length(xobs), -mindist, mindist)
  # sample M of these points
  grid_obs <- grid_obs[sample(nrow(grid_obs), size = M), ]
  # concatenate these points and indicate whether they are prediction or
  # observation locations
  grid <- rbind(grid_obs, grid_pred)
  h <- as.matrix(dist(grid))
  dimnames(h) <- NULL
  grid$type <- factor(c(rep('obs', M), rep('pred', ngrid * ngrid)))
  # split the distance matrix into prediction/observation locations only
  h_obs <- h[grid$type == 'obs',grid$type == 'obs']
  h_pred <- h[grid$type == 'pred',grid$type == 'pred']
  # create output
  out <- list(locs = grid, dist_obs = h_obs, dist_pred = h_pred, mindist = mindist)
  class(out) <- 'GPlocations'
  return(out)
}


prepare_model <- function(family, params) {
  if (family == 'matern') {
    names(params) <- c('nu', 'rho', 'sigma', 'nugget')
    model <- RMhandcock(params['nu'], notinvnu = TRUE, scale = params['rho'], var = params['sigma']) + 
      RMnugget(var = params['nugget'])
    specdens <- function(w) dmatern(w, nu = params['nu'], alpha = 1/params['rho'], sigma = params['sigma'])
  } else if (family == 'dampedcos') {
    names(params) <- c('lambda', 'theta', 'sigma', 'nugget')
    model <- RMdampedcos(lambda = params['lambda'], scale = params['theta'], var = params['sigma']) +
      RMnugget(var = params['nugget'])
    integrand <- function(h, w) {
      besselJ(h*w, 0) * h * RFcov(model, h)
    }
    specdens <- Vectorize(function(w) {
      int <- integrate(integrand, w = w, 0, Inf)$value
      1/(2*pi) * int
    }, vectorize.args = 'w')
  }
  
  out <- list(model = model, params = params, specdens = specdens)
  class(out) <- 'GPmodel'
  return(out)
}


simulate_gp <- function(locations, family, params) {
  locations$m <- prepare_model(family, params)
  locations$Y <- as.vector(RFsimulate(locations$m$model, locations$locs$x, locations$locs$y))
  class(locations) <- 'GPsimulated'
  return(locations)
}


plot.GPlocations <- function(locobj, ...) {
  pred_idx <- which(locobj$locs$type == 'pred')
  grid_pred <- locobj$locs[pred_idx, ]
  grid_obs <- locobj$locs[-pred_idx, ]
  plot.default(grid_pred$x, grid_pred$y, col = 'grey')
  points(grid_obs$x, grid_obs$y, pch = 19)
}

plot.GPsimulated <- function(obj, ...) {
  pred_idx <- which(obj$locs$type == 'pred')
  grid_pred <- obj$locs[pred_idx, ]
  grid_obs <- obj$locs[-pred_idx, ]
  Y_pred <- obj$Y[pred_idx]
  Y_obs <- obj$Y[-pred_idx]
  image(unique(grid_pred$x),
        unique(grid_pred$y),
        matrix(Y_pred, nrow = length(unique(grid_pred$x)), byrow = T),
        col = rev(brewer.pal(7, 'Blues')),
        xlab = '', ylab = '')
  points(grid_obs$x, grid_obs$y, col = 'white', pch = 20)
  points(grid_obs$x, grid_obs$y, col = 'black', pch = 1)
}