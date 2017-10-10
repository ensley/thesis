#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)

args <- commandArgs(TRUE)
if(length(args) != 4) {
  stop("incorrect arguments")
}

filepath <- args[1]
filename <- args[2]
outpath <- args[3]
nbatch <- as.numeric(args[4])

args <- readRDS(file.path(filepath, filename))

batch <- estbeta_aao_initialize(file.path(outpath, 'allbetas.csv'),
                            file.path(outpath, 'errbounds.csv'),
                            nbatch,
                            args$gp$dist_obs,
                            args$gp$Y[args$gp$locs$type == 'obs'],
                            args$B,
                            args$beta_init,
                            args$tau_init,
                            args$knots,
                            args$sigma.m,
                            args$v_tau,
                            args$t_v,
                            r.opt = 0.23,
                            r = 6,
                            eps = 1e-3,
                            args$nugget)

saveRDS(batch, file.path(outpath, 'batch.rds'))

