#! /usr/bin --vanilla --default-packages=utils,gpcovr


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



# COVARIANCE PLOT ---------------------------------------------------------

pdf('covariance.pdf')
covar_plot(x,
           bf,
           args$knots,
           gpModel = args$gp$m,
           beta_true = args$beta_true,
           err_bounds = cov_err_bounds2)
dev.off()

### END COVARIANCE PLOT


gpobj <- readRDS(file.path(fileinpath, 'gpobj.rds'))
# predict(gpobj, beta = bf, knots = args$knots)
grid <- gpobj$locs[ ,1:2]
is_pred <- which(gpobj$locs$type == 'pred')
pts_obs <- gpobj$Y[-is_pred]
dist_pred2obs <- as.matrix(stats::dist(grid))[is_pred,-is_pred]
draws <- exp(rlogspline(50000, bf, args$knots))
gamma1 <- matrix(gpumc::mc(dist_pred2obs, draws), nrow = length(is_pred))
G <- matrix(gpumc::mc(gpobj$dist_obs, draws), nrow = length(pts_obs))
diag(G) <- 1 + 1e-5
Ginv <- solve(G)
preds_spline <- drop(gamma1 %*% Ginv %*% pts_obs)
gpcovr::mse(preds, predYpred)

# FIT BEST MATERN ---------------------------------------------------------

starting = list('sigma.sq' = 1, 'tau.sq' = 0.1, 'phi' = 0.5, 'nu' = 1)
tuning = list('sigma.sq' = 0.01, 'tau.sq' = 0.01, 'phi' = 0.005, 'nu' = 0.005)
priors = list('sigma.sq.ig' = c(2, 1), 'tau.sq.ig' = c(2, 1), 'phi.unif' = c(1e-5, 20), 'nu.unif' = c(1e-5, 5))
n_samps <- 5000
is_obs <- which(gpobj$locs$type == 'obs')
coords <- as.matrix(gpobj$locs[is_obs,1:2])
y <- gpobj$Y[is_obs]
burn.in <- floor(0.75*n_samps)

m_i <- spBayes::spLM(y ~ 1, coords = coords, starting = starting, tuning = tuning, priors = priors, cov.model = 'matern', n.samples = n_samps, n.report = 200)

pars <- apply(m_i$p.theta.samples[burn.in:nrow(m_i$p.theta.samples), ], 2, stats::median)
model <- RandomFields::RMwhittle(nu = pars['nu'],
                                 notinvnu = TRUE,
                                 var = pars['sigma.sq'],
                                 scale = 1/pars['phi']) +
  RandomFields::RMnugget(var = pars['tau.sq'])

params <- c(pars['nu'],
            (2 * sqrt(pars['nu'])) / pars['phi'],
            pars['sigma.sq'],
            pars['tau.sq'])
names(params) <- c('nu', 'rho', 'sigma', 'nugget')

fit <- list(model = model, fitobj = m_i, params = params)

# COMPUTE PREDICTIONS -----------------------------------------------------

preds <- data.frame(actual = predYpred,
                    spline = preds_spline,
                    bestmatern = othermodel_predict(gpobj, fit$model))

write.csv(preds, 'preds.csv')

# COMPUTE MSE -------------------------------------------------------------

mse_spline <- mse(preds$spline, preds$actual)
mse_bestmatern <- mse(preds$bestmatern, preds$actual)

# PLOT PREDICTIONS --------------------------------------------------------

pdf('pairs.pdf')
pairs(preds, pch = 20, cex = 0.5)
dev.off()
