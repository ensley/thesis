#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)

args <- commandArgs(TRUE)

if (length(args) < 8) {
  stop('Not enough arguments to 02a')
}

filepath <- args[1]
family <- args[2]
M <- as.numeric(args[3])          # number of observations in the Gaussian process
ngrid <- as.numeric(args[4])      # number of observations along one edge of prediction grid
K <- as.numeric(args[5])          # number of knots (and therefore number of betas)
B <- as.numeric(args[6])          # number of samples to take from the spectral density
seed <- args[7]                   # will the random seed be set?
params <- as.numeric(args[8:length(args)])

if (seed == 'y') {
  set.seed(1234)
  cat('SEED SET.\n')
}



cat('CREATING LOCATIONS...')



locations <- create_locations(M, ngrid)



cat('DONE.\nSIMULATING OBSERVATIONS...')



gp <- simulate_gp(locations, family, params)



cat('DONE.\n')



knots <- seq(gp$m$knotlocs[1], gp$m$knotlocs[2], length = K)
x <- exp(seq(gp$m$logx[1], gp$m$logx[2], length = 500))
y <- gp$m$specdens(x)
logx <- log(x)
logy <- logx + log(y)
w <- round(seq(1, length(x), length = 50))
xx <- logx[w]
yy <- logy[w]
spl <- nsbasis(xx, knots)
cat('KNOT LOCATIONS: ', knots)


cat('FINDING TRUE BETA VECTOR...')



f <- fitspline(500000, yy, spl)
beta_true <- colMeans(f$b)
rm(f)



cat('DONE.\nDETERMINING INITIAL BETA VECTOR...\n')



invalid <- TRUE
while(invalid) {
  beta_init <- beta_true + runif(K, -0.5, 0.5)
  slopes <- get_slopes(beta_init, knots)
  invalid <- slopes[1] <= 0 || slopes[2] >= 0
  cat('\tANOTHER TRY\n')
}



cat('INITIAL BETA DETERMINED.\n')



tau_init <- 1
v_beta <- rep(1e-4, K)
v_tau <- 1
t_v <- 5
nugget <- 1e-2



cat('INITIAL VALUES: ', beta_init, '\n')



args <- list(locations = locations,
             gp = gp,
             B = B,
             knots = knots,
             beta_true = beta_true,
             beta_init = beta_init,
             tau_init = tau_init,
             v_beta = v_beta,
             v_tau = v_tau,
             t_v = t_v,
             nugget = nugget,
             x = x,
             y = y)


saveRDS(args, file.path(filepath, 'init_args.rds'))
print('initial data saved to init_args.rds')
