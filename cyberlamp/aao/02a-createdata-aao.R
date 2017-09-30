#! /usr/bin --vanilla --default-packages=utils,openacc

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
locations <- create_locations(M, ngrid)
cat('LOCATIONS CREATED.\n')
gp <- simulate_gp(locations, family, params)
cat('OBSERVATIONS SIMULATED.\n')

knots <- seq(-1, 2, length = K)
x <- exp(seq(-5, 3, length = 500))
y <- gp$m$specdens(x)
logx <- log(x)
logy <- logx + log(y)
w <- round(seq(1, length(x), length = 50))
xx <- logx[w]
yy <- logy[w]

spl <- nsbasis(xx, knots)
f <- fitspline(500000, yy, spl)
beta_true <- colMeans(f$b)
rm(f)

cat('TRUE BETA VECTOR FOUND.\n')


# spl <- nsbasis(logx, knots)
# plot(xx, yy, pch = 20)
# abline(v = knots)
# lines(logx, predict_natspl(nsbasis(logx, knots), beta_true))
# legend('topleft', c('final estimate', '"Truth"'), col = c('black', 'blue'), lty = 1)
# abline(v = args$knots, lty = 3, col = 'grey')


invalid <- TRUE
while(invalid) {
  beta_init <- beta_true + runif(K, -0.5, 0.5)
  slopes <- get_slopes(beta_init, knots)
  invalid <- slopes[1] <= 0 || slopes[2] >= 0
  cat('ANOTHER TRY\n')
}

# lines(logx, predict_natspl(nsbasis(logx, knots), beta_init))

cat('INITIAL BETA DETERMINED.\n')

tau_init <- 1
sigma.m <- 1e-3
v_tau <- 1
t_v <- 5
nugget <- 1e-2

cat('INITIAL VALUES: ', beta_init, '\n')

####### THIS IS FOR A MATERN COVARIANCE FUNCTION

# set.seed(1234)
# # STEP 1: create observation points
# mindist <- 0.005
# x <- seq(mindist, 1-mindist, by = 3 * mindist)
# xpred <- seq(0, 1, length = 50)
# grid_all <- expand.grid(x = x, y = x) + runif(2 * length(x)^2, -mindist, mindist)
# grid <- grid_all[sample(1:nrow(grid_all), size = 400), ]
# grid <- rbind(grid, expand.grid(x = xpred, y = xpred))
# grid$xtype <- factor(c(rep('obs', 400), rep('pred', 50*50)))
# h <- as.matrix(dist(grid[ ,-3])); dimnames(h) <- NULL
# hobs <- h[grid$xtype == 'obs',grid$xtype == 'obs']
# hpred <- h[grid$xtype == 'pred',grid$xtype == 'pred']
# # STEP 2: simulate from model
# covmat <- matern_cor(h, 1.5, 4.76, 1)
# diag(covmat) <- 1.01
# R <- chol(covmat)
# u <- rnorm(nrow(grid))
# z <- crossprod(R, u)
# Yobs <- z[grid$xtype == 'obs']
# Ypred <- z[grid$xtype == 'pred']
# # STEP 3: place the knots
# knots <- seq(0, 3, length = K)
#
# msg <- 'Gaussian process generated with generateGP(). Matern covariance function. mindist set to 0.002, tuning variance sigma.m fixed at 1e-4.'
#
# x <- exp(seq(-5, 5, length = 500))
# y <- dmatern(x, 1.5, 4.76, 1)
# logx <- log(x)
# logy <- logx + log(y)
# w <- round(seq(1, length(x), length = 50))
# xx <- logx[w]
# yy <- logy[w]
#
# spl <- nsbasis(xx, knots)
# f <- fitspline(500000, yy, spl)
# beta_true <- colMeans(f$b)
# rm(f)
#
#
# beta_init <- beta_true + runif(K, -0.5, 0.5)
# tau_init <- 1
# sigma.m <- 1e-3
# v_tau <- 1
# t_v <- 5
# nugget <- 1e-2
#
# beta_true
# beta_init

####### END MATERN

####### THIS IS FOR A DAMPED COSINE COVARIANCE FUNCTION
# set.seed(1234)
# # create observation points
# mindist <- 0.002
# x <- seq(mindist, 1-mindist, by = 3 * mindist)
# xpred <- seq(0, 1, length = 50)
# grid_all <- expand.grid(x = x, y = x) + runif(2 * length(x)^2, -mindist, mindist)
# grid <- grid_all[sample(1:nrow(grid_all), size = 400), ]
# grid <- rbind(grid, expand.grid(x = xpred, y = xpred))
# grid$xtype <- factor(c(rep('obs', 400), rep('pred', 50*50)))
# h <- as.matrix(dist(grid[ ,-3])); dimnames(h) <- NULL
# hobs <- h[grid$xtype == 'obs',grid$xtype == 'obs']
# hpred <- h[grid$xtype == 'pred',grid$xtype == 'pred']
# # simulate from model
# model_sim <- RMdampedcos(lambda = 1, var = 2, scale = 0.7) + RMnugget(var = 0.01)
# Y <- RFsimulate(model_sim, grid$x, grid$y)@data$variable1
# Yobs <- Y[grid$xtype == 'obs']
# Ypred <- Y[grid$xtype == 'pred']
# # z <- RFsimulate(model_sim, seq(0.005, 0.995, length = 100), seq(0.005, 0.995, length = 100))
#
#
# integrand <- function(h, w, lambda) {
#   besselJ(h*w, 0) * h * exp(-lambda * h) * cos(h)
# }
#
# specdens <- Vectorize(function(w, lambda) {
#   int <- integrate(integrand, w = w, lambda = lambda, 0, Inf)$value
#   1/(2*pi) * int
# }, vectorize.args = 'w')
#
#
# logx <- seq(-5, 3, length = 1000)
# knots <- seq(-0.5, 1.5, length = 15)
# x <- exp(logx)
# y <- specdens(x, 1)
# # plot(x, s, type = 'l')
#
# logy <- log(y) + logx
#
# # plot(logx, logy, type = 'l')
#
# w <- floor(seq(1, length(x), length = 50))
# xx <- logx[w]
# yy <- logy[w]
#
# plot(xx, yy, pch = 20)
# abline(v = knots, lty = 3, col = 'grey')
#
# spl <- nsbasis(xx, knots)
# f <- fitspline(100000, yy, spl)
# beta_true <- colMeans(f$b)
#
# # set.seed(123456)
# beta_init <- beta_true + runif(K, -0.5, 0.5)
# tau_init <- 1
# sigma.m <- 1e-4
# prop.Sigma <- 1e-3 * diag(K)
# v_tau = 1
# t_v = 5
# nugget = 1e-2
#
# beta_true
# beta_init
#
# msg <- 'Exponentially damped cosine covariance function. Tuning variance sigma.m fixed at 1e-4.'
#
# lines(logx, predict_natspl(nsbasis(logx, knots), betatrue))


####### END DAMPED COSINE

args <- list(locations = locations,
             gp = gp,
             B = B,
             knots = knots,
             beta_true = beta_true,
             beta_init = beta_init,
             tau_init = tau_init,
             sigma.m = sigma.m,
             v_tau = v_tau,
             t_v = t_v,
             nugget = nugget,
             x = x,
             y = y)

# args <- list(hobs = locations$dist_obs,
#              hpred = locations$dist_pred,
#              grid = locations$locs,
#              Yobs = gp$Y[gp$locs$type == 'obs'],
#              Ypred = gp$Y[gp$locs$type == 'pred'],
#              B = B,
#              knots = knots,
#              beta_true = beta_true,
#              beta_init = beta_init,
#              tau_init = tau_init,
#              sigma.m = sigma.m,
#              v_tau = v_tau,
#              t_v = t_v,
#              nugget = nugget,
#              x = x,
#              y = y,
#              mindist = locations$mindist)

saveRDS(args, file.path(filepath, 'init_args.rds'))
print('initial data saved to init_args.rds')
