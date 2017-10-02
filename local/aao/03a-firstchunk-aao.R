#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)

args <- commandArgs(TRUE)
if(length(args) != 2) {
  stop("incorrect arguments")
}

filepath <- args[1]
nbatch <- as.numeric(args[2])

args <- readRDS(file.path(filepath, 'init_args.rds'))

batch <- estbeta_aao_initialize(file.path(filepath, 'allbetas.csv'),
                            file.path(filepath, 'errbounds.csv'),
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

saveRDS(batch, file.path(filepath, 'batch.rds'))
