#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)

args <- commandArgs(TRUE)
if(length(args) != 3) {
  stop("incorrect arguments")
}

filepath <- args[1]
outpath <- args[2]
nbatch <- as.numeric(args[3])

args <- readRDS(file.path(filepath))

batch <- estbeta_aao_initialize(file.path(filepath, 'allbetas.csv'),
                            file.path(filepath, 'errbounds.csv'),
                            nbatch,
                            args$locations$dist_obs,
                            args$gp$Y[args$gp$locs$type == 'obs'],
                            args$B,
                            args$beta_init,
                            args$tau_init,
                            args$knots,
                            args$v_beta,
                            args$v_tau,
                            args$t_v,
                            r.opt = 0.23,
                            r = 6,
                            eps = 1e-3,
                            args$nugget)

saveRDS(batch, file.path(outpath, 'batch.rds'))

