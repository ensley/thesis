#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)

args <- commandArgs(TRUE)
if(length(args) != 3) {
  stop("incorrect arguments")
}

filepath <- args[1]
nbatch <- as.numeric(args[2])
batchidx <- as.numeric(args[3])

b <- readRDS(file.path(filepath, 'batch.rds'))

batch <- estbeta(file.path(filepath, 'allbetas.csv'),
                 file.path(filepath, 'errbounds.csv'),
                 nbatch,
                 batchidx,
                 b$args$H,
                 b$args$Y,
                 b$args$B,
                 b$last_few_beta,
                 b$last_few_tau,
                 b$prev_lik,
                 b$last_accept,
                 b$args$knots,
                 b$v_beta,
                 b$v_tau,
                 b$args$t_v,
                 b$args$r.opt,
                 b$args$eps,
                 b$args$nugget)

saveRDS(batch, file.path(filepath, 'batch.rds'))
