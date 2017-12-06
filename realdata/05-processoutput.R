#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(tidyverse)
library(gpcovr)

args <- commandArgs(TRUE)

if (length(args) < 3) {
  stop('Not enough arguments to 05a')
}

fileinpath <- args[1]
fileinname <- args[2]
outpath <- args[3]

setwd(outpath)
args <- readRDS(file.path(fileinpath, fileinname))
b <- read.csv(file.path(outpath, 'allbetas.csv'), header = FALSE)
x <- seq(1e-5, 5, length = 500)
# err_bounds <- get_error_bounds(b, args$gp, args$knots, thin = 1000)
cov_err_bounds <- read.csv(file.path(outpath, 'errbounds.csv'), header = FALSE)
cov_err_bounds2 <- t(apply(cov_err_bounds, 2, function(x) stats::quantile(x, c(0.025, 0.975))))


N <- nrow(b)
# burnin <- floor(0.2 * N)
burnin <- 1
bf <- colMeans(tail(b, -burnin))
apply(b, 2, function(x) mean(diff(x) != 0))






# LOCATIONS PLOT ----------------------------------------------------------

# pdf('locations.pdf')
# plot(args$gp)
# dev.off()

### END LOCATIONS PLOT




# TRACE PLOT --------------------------------------------------------------

pdf('trace.pdf')
trace_plot(b, 4)
dev.off()

### END TRACE PLOT




# LOG-LOG DENSITY PLOT ----------------------------------------------------

pdf('loglogdensity.pdf')
loglogdensity_plot(bf,
                   args$knots,
                   show_knots = TRUE,
                   show_points = TRUE,
                   beta_init = args$beta_init)
dev.off()

### END LOG-LOG DENSITY PLOT



# LOG DENSITY PLOT --------------------------------------------------------

pdf('logdensity.pdf')
logdensity_plot(bf,
                args$knots,
                show_knots = TRUE,
                show_points = TRUE,
                beta_init = args$beta_init)
dev.off()

### END LOG DENSITY PLOT


gpobj <- readRDS(file.path(fileinpath, 'gpobj.rds'))
# predict(gpobj, beta = bf, knots = args$knots)
grid <- gpobj$locs[ ,1:2]
is_pred <- which(gpobj$locs$type == 'pred')
pts_obs <- gpobj$Y[-is_pred]
pts_pred <- gpobj$Y[is_pred]
dist_pred2obs <- as.matrix(stats::dist(grid))[is_pred,-is_pred]
draws <- exp(rlogspline(50000, bf, args$knots))
gamma1 <- matrix(gpumc::mc(dist_pred2obs, draws), nrow = length(is_pred))
G <- matrix(gpumc::mc(gpobj$dist_obs, draws), nrow = length(pts_obs))
mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
if (mineigen < 0) diag(G) <- diag(G) + abs(mineigen) + 1e-4
Ginv <- solve(G)
preds_spline <- drop(gamma1 %*% Ginv %*% pts_obs)
gpcovr::mse(preds_spline, pts_pred)

# FIT BEST MATERN ---------------------------------------------------------

args <- readRDS(file.path(fileinpath, 'init_args.rds'))

# starting = list('sigma.sq' = 1, 'tau.sq' = 0.1, 'phi' = 0.5, 'nu' = 1)
# tuning = list('sigma.sq' = 0.01, 'tau.sq' = 0.01, 'phi' = 0.005, 'nu' = 0.005)
# priors = list('sigma.sq.ig' = c(2, 1), 'tau.sq.ig' = c(2, 1), 'phi.unif' = c(1e-5, 20), 'nu.unif' = c(1e-5, 5))
# n_samps <- 5000
# is_obs <- which(gpobj$locs$type == 'obs')
# coords <- as.matrix(gpobj$locs[is_obs,1:2])
# y <- gpobj$Y[is_obs]
# burn.in <- floor(0.75*n_samps)
# 
# m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'matern', n.samples = n_samps, n.report = 200)
# 
# pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
# model <- RandomFields::RMwhittle(nu = pars['nu'],
#                                  notinvnu = TRUE,
#                                  var = pars['sigma.sq'],
#                                  scale = 1/pars['phi']) +
#   RandomFields::RMnugget(var = pars['tau.sq'])
# 
# params <- c(pars['nu'],
#             (2 * sqrt(pars['nu'])) / pars['phi'],
#             pars['sigma.sq'],
#             pars['tau.sq'])
# names(params) <- c('nu', 'rho', 'sigma', 'nugget')
# 
# fit <- list(model = model, fitobj = m_i, params = params)

results <- readRDS('~/git/thesis/results/realdata/try4/result_true.rds')
burnin <- 300
m <- coda::as.mcmc(results$samps)
params <- colMeans(window(m, start = burnin))
model_matern <- RandomFields::RMwhittle(nu = params['nu'],
                                 notinvnu = TRUE,
                                 var = params['sigma'],
                                 scale = 1/params['alpha']) +
  RandomFields::RMnugget(var = 0.01)


gamma1 <- matrix(RandomFields::RFcov(model_matern, as.numeric(dist_pred2obs)),
                 nrow = length(is_pred))
G <- matrix(RandomFields::RFcov(model_matern, as.numeric(gpobj$dist_obs)),
            nrow = length(pts_obs))
mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
if (mineigen < 0) diag(G) <- diag(G) + mineigen + 1e-4
Ginv <- solve(G)
preds_bestmatern <- drop(gamma1 %*% Ginv %*% pts_obs)
lines(x, RandomFields::RFcov(model_matern, x))



results <- readRDS('~/git/thesis/results/realdata/trydamped1/result_true.rds')
burnin <- 300
m <- coda::as.mcmc(results$samps)
params <- colMeans(window(m, start = burnin))
model_damped <- RandomFields::RMdampedcos(lambda = params['lambda'],
                                   var = params['var'],
                                   scale = params['scale']) +
  RandomFields::RMnugget(var = 0.01)


gamma1 <- matrix(RandomFields::RFcov(model_damped, as.numeric(dist_pred2obs)),
                 nrow = length(is_pred))
G <- matrix(RandomFields::RFcov(model_damped, as.numeric(gpobj$dist_obs)),
            nrow = length(pts_obs))
mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
if (mineigen < 0) diag(G) <- diag(G) + mineigen + 1e-4
Ginv <- solve(G)
preds_bestdamped <- drop(gamma1 %*% Ginv %*% pts_obs)
lines(x, RandomFields::RFcov(model_damped, x))


# COVARIANCE PLOT ---------------------------------------------------------

covar_plot <- function(h, beta, knots, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  l1 <- gpumc::mc(h, exp(rlogspline(50000, beta, knots)))
  l2 <- if (!is.null(beta_true)) gpumc::mc(h, exp(rlogspline(50000, beta_true, knots))) else NULL
  l3 <- if (!is.null(gpModel)) RandomFields::RFcov(gpModel$model, h)
  
  if (!is.null(err_bounds)) {
    min_lower_bd <- min(err_bounds[ ,1])
    min_upper_bd <- min(err_bounds[ ,2])
    max_lower_bd <- max(err_bounds[ ,1])
    max_upper_bd <- max(err_bounds[ ,2])
  } else {
    min_lower_bd <- min_upper_bd <- max_lower_bd <- max_upper_bd <- NULL
  }
  
  plot_ymin <- min(min_lower_bd, min_upper_bd, min(l1), suppressWarnings(min(l2)), suppressWarnings(min(l3)))
  plot_ymax <- max(max_lower_bd, max_upper_bd, max(l1), suppressWarnings(max(l2)), suppressWarnings(max(l3)))
  
  graphics::plot(range(h), c(plot_ymin, plot_ymax), type = 'n',
                 xlab = 'h',
                 ylab = 'C(h)',
                 main = 'Covariance function')
  
  if (!is.null(err_bounds)) graphics::polygon(c(h, rev(h)), c(err_bounds[ ,1], rev(err_bounds[ ,2])), col = 'grey90', border = NA)
  graphics::lines(h, l1, lwd = 2, col = 'blue')
  if (!is.null(beta_true)) graphics::lines(h, l2, lwd = 2, lty = 2)
  if (!is.null(gpModel)) graphics::lines(h, l3, col = 'orange', lwd = 2, lty = 2)
  graphics::abline(h = 0, lty = 3)
  # graphics::legend('topright', c('estimated', 'optimal', 'true'), col = c('blue', 'black', 'orange'), lwd = 2, lty = c(1, 2, 2))
  
  invisible(gpModel)
}

pdf('covariance_with_others.pdf')
covar_plot(x,
           bf,
           args$knots,
           gpModel = args$gp$m,
           beta_true = args$beta_true,
           err_bounds = cov_err_bounds2)
lines(x, RandomFields::RFcov(model_matern, x))
lines(x, RandomFields::RFcov(model_damped, x))
dev.off()

### END COVARIANCE PLOT


# m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'exponential', n.samples = n_samps, n.report = 200)
# 
# pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
# model <- RandomFields::RMexp(var = pars['sigma.sq'],
#                              scale = 1/pars['phi']) +
#   RandomFields::RMnugget(var = pars['tau.sq'])
# 
# gamma1 <- matrix(RandomFields::RFcov(model, as.numeric(dist_pred2obs)),
#                  nrow = length(is_pred))
# G <- matrix(RandomFields::RFcov(model, as.numeric(gpobj$dist_obs)),
#             nrow = length(pts_obs))
# mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
# if (mineigen < 0) diag(G) <- diag(G) + mineigen + 1e-4
# Ginv <- solve(G)
# preds_bestexp <- drop(gamma1 %*% Ginv %*% pts_obs)
# 
# 
# 
# m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'spherical', n.samples = n_samps, n.report = 200)
# 
# pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
# model <- RandomFields::RMexp(var = pars['sigma.sq'],
#                              scale = 1/pars['phi']) +
#   RandomFields::RMnugget(var = pars['tau.sq'])
# 
# gamma1 <- matrix(RandomFields::RFcov(model, as.numeric(dist_pred2obs)),
#                  nrow = length(is_pred))
# G <- matrix(RandomFields::RFcov(model, as.numeric(gpobj$dist_obs)),
#             nrow = length(pts_obs))
# mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
# if (mineigen < 0) diag(G) <- diag(G) + mineigen + 1e-4
# Ginv <- solve(G)
# preds_bestspherical <- drop(gamma1 %*% Ginv %*% pts_obs)
# 
# 
# 
# 
# m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'gaussian', n.samples = n_samps, n.report = 200)
# 
# pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
# model <- RandomFields::RMexp(var = pars['sigma.sq'],
#                              scale = 1/pars['phi']) +
#   RandomFields::RMnugget(var = pars['tau.sq'])
# 
# gamma1 <- matrix(RandomFields::RFcov(model, as.numeric(dist_pred2obs)),
#                  nrow = length(is_pred))
# G <- matrix(RandomFields::RFcov(model, as.numeric(gpobj$dist_obs)),
#             nrow = length(pts_obs))
# mineigen <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
# if (mineigen < 0) diag(G) <- diag(G) + mineigen + 1e-4
# Ginv <- solve(G)
# preds_bestgaussian <- drop(gamma1 %*% Ginv %*% pts_obs)



# COMPUTE PREDICTIONS -----------------------------------------------------

preds <- data.frame(actual = predYpred,
                    spline = preds_spline,
                    bestmatern = preds_bestmatern,
                    bestexp = preds_bestexp,
                    bestspherical = preds_bestspherical,
                    bestgaussian = preds_bestgaussian)

write.csv(preds, 'preds.csv')

# COMPUTE MSE -------------------------------------------------------------

mse_spline <- mse(preds$spline, preds$actual)
mse_bestmatern <- mse(preds$bestmatern, preds$actual)
mse_bestexponential <- mse(preds_bestexp, preds$actual)
mse_bestgaussian <- mse(preds_bestgaussian, preds$actual)
mse_bestspherical <- mse(preds_bestspherical, preds$actual)

mses <- data.frame(type = c('spline', 'best matern', 'best exponential', 'best gaussian', 'best spherical'),
                   mse = c(mse_spline, mse_bestmatern, mse_bestexponential, mse_bestgaussian, mse_bestspherical))
mses %>% mutate(mse_rel = mse / mse_spline) %>% 
  arrange(mse_rel) %>% knitr::kable(format = 'latex', booktabs = TRUE)


# PLOT PREDICTIONS --------------------------------------------------------

pdf('pairs.pdf')
pairs(preds, pch = 20, cex = 0.5)
dev.off()


preds %>% 
  select(-X) %>% 
  gather(key = 'type', value = 'value', spline:bestgaussian) %>% 
  mutate(residual = value - actual) %>% 
  ggplot(aes(actual, residual, color = type)) +
  geom_point()
