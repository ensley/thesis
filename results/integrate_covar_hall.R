#! /usr/bin --vanilla --default-packages=utils,gpcovr
library(gpcovr)
Rcpp::sourceCpp('~/Dropbox/psu/research/specdens/hallpaper.cpp')

cmd_args <- commandArgs(TRUE)
if (length(cmd_args) != 4) stop('incorrect number of arguments')

datadir <- cmd_args[1]
outdir <- cmd_args[2]
start <- as.numeric(cmd_args[3])
end <- as.numeric(cmd_args[4])



for (i in start:end) {
  cat('\n\ni=', i, '\n\n')
  if (i < 10) {
    dataset_folder <- paste0('dataset00', i)
  } else if (i < 100) {
    dataset_folder <- paste0('dataset0', i)
  } else {
    dataset_folder <- paste0('dataset', i)
  }
  
  args <- readRDS(file.path(datadir, paste0('init_args_', i, '.rds')))
  
  t <- theta <- seq(0, 20, length = 200)
  X <- args$gp$Y[args$gp$locs$type == 'obs']
  d <- as.matrix(args$gp$locs[args$gp$locs$type == 'obs',1:2])
  h <- 0.1
  d_dist <- as.matrix(dist(d))
  dd <- as.numeric(d_dist)
  
  rho_hat(0.01, X, dd, h)
  y <- sapply(t, function(tt) rho_hat(tt, X, dd, h))
  # plot(t, y, type = 'l', main = 'rho hat', xlab = 't', ylab = 'rho_hat(t)')
  
  yy <- sapply(t, function(tt) rho1_hat(tt, X, dd, h, 5, 5.1))
  # lines(t, yy, col = 'blue')
  
  yyy <- th <- c(); nextth <- 0; th_step <- 0.01
  while( (nexty <- rho_cross(nextth, X, dd, h, 5, 5.1)) >= 0) {
    cat('th =', nextth, '\n')
    yyy <- c(yyy, nexty)
    th <- c(th, nextth)
    nextth <- nextth + th_step
  }
  # plot(th, yyy, type = 'l')
  
  theta_new <- c(-rev(th), th)
  yyyy <- c(rev(yyy), yyy)
  trap <- function(t, theta, y) {
    yyy <- y * cos(theta * t) / (2*pi)
    as.double((tail(theta, -1) - head(theta, -1)) %*% (tail(yyy, -1) + head(yyy, -1))) / 2
  }
  z <- sapply(theta, function(th) trap(th, theta_new, yyy))
  # plot(theta, z / z[1], type = 'l', ylim = c(-1, 1), lwd = 2)
  # lines(theta, RandomFields::RFcov(args$gp$m$model, theta), lty = 2)
  # abline(h = 0)
  
  cov_hall <- stats::splinefun(theta, z)
  cov_hall_norm <- stats::splinefun(theta, z / z[1])
  
  cov <- read.csv(file.path(outdir, dataset_folder, 'covariance.csv'))
  cov$cov_hall <- cov_hall(cov$h)
  cov$cov_hall_norm <- cov_hall_norm(cov$h)
  write.csv(cov, file.path(outdir, dataset_folder, 'covariance.csv'), row.names = FALSE)
}