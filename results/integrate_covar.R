#! /usr/bin --vanilla --default-packages=utils,gpcovr
library(gpcovr)



cmd_args <- commandArgs(TRUE)
if (length(cmd_args) != 4) stop('incorrect number of arguments')

# datadir <- '~/git/thesis/data/matern'
# outdir <- '~/git/thesis/results/matern/aao'

datadir <- cmd_args[1]
outdir <- cmd_args[2]
start <- as.numeric(cmd_args[3])
end <- as.numeric(cmd_args[4])

burnin <- 1
h <- seq(0, 5, length = 500)

for (i in start:end) {
  if (i < 10) {
    dataset_folder <- paste0('dataset00', i)
  } else if (i < 100) {
    dataset_folder <- paste0('dataset0', i)
  } else {
    dataset_folder <- paste0('dataset', i)
  }
  
  args <- readRDS(file.path(datadir, paste0('init_args_', i, '.rds')))
  b <- read.csv(file.path(outdir, dataset_folder, 'allbetas.csv'), header = FALSE)
  bf <- colMeans(tail(b, -burnin))
  
  fit <- fit_matern(args$gp, n_samps = 5000)
  
  cov_true <- RandomFields::RFcov(args$gp$m$model, h)
  cov_spline <- gpumc::mc(h, exp(rlogspline(50000, bf, args$knots)))
  cov_bestmatern <- RandomFields::RFcov(fit$model, h)
  
  cov <-  data.frame(h, cov_true, cov_spline, cov_bestmatern)
  
  saveRDS(fit, file.path(outdir, dataset_folder, 'bestmatern.rds'))
  write.csv(cov, file.path(outdir, dataset_folder, 'covariance.csv'), row.names = FALSE)
}
