#! /usr/bin --vanilla --default-packages=utils,gpcovr


library(gpcovr)

args <- commandArgs(TRUE)

if (length(args) < 2) {
  stop('Not enough arguments to 05a')
}

filepath <- args[1]
family <- args[2]

setwd(filepath)
args <- readRDS('init_args.rds')
b <- read.csv('allbetas.csv', header = FALSE)
x <- seq(1e-5, 5, length = 500)
err_bounds <- get_error_bounds(b, args$gp, args$knots, thin = 1000)
cov_err_bounds <- read.csv('errbounds.csv', header = FALSE)
cov_err_bounds2 <- t(apply(cov_err_bounds, 2, function(x) stats::quantile(x, c(0.025, 0.975))))


N <- nrow(b)
burnin <- floor(0.5 * N)
bf <- colMeans(tail(b, -burnin))
apply(tail(b, -burnin), 2, function(x) mean(diff(x) != 0))






# LOCATIONS PLOT ----------------------------------------------------------

pdf('locations.pdf')
plot(args$gp)
dev.off()

### END LOCATIONS PLOT




# TRACE PLOT --------------------------------------------------------------

pdf('trace.pdf')
trace_plot(tail(b, -burnin), 4)
dev.off()

### END TRACE PLOT




# LOG-LOG DENSITY PLOT ----------------------------------------------------

pdf('loglogdensity.pdf')
loglogdensity_plot(bf,
                   args$knots,
                   show_knots = TRUE,
                   show_points = TRUE,
                   gpModel = args$gp$m,
                   beta_true = args$beta_true,
                   beta_init = args$beta_init,
                   err_bounds = err_bounds)
dev.off()

### END LOG-LOG DENSITY PLOT



# LOG DENSITY PLOT --------------------------------------------------------

pdf('logdensity.pdf')
logdensity_plot(bf,
                args$knots,
                show_knots = TRUE,
                show_points = TRUE,
                gpModel = args$gp$m,
                beta_true = args$beta_true,
                err_bounds = err_bounds)
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




# FIT BEST MATERN ---------------------------------------------------------

fit <- fit_matern(args$gp, 5000)

# COMPUTE PREDICTIONS -----------------------------------------------------

preds <- data.frame(actual = args$gp$Y[args$gp$locs$type == 'pred'],
                    spline = predict(args$gp, beta = bf, knots = args$knots),
                    bestmatern = othermodel_predict(args$gp, fit$model),
                    optimal = optimal_predict(args$gp))

write.csv(preds, 'preds.csv')

# COMPUTE MSE -------------------------------------------------------------

mse_spline <- mse(preds$spline, preds$actual)
mse_bestmatern <- mse(preds$bestmatern, preds$actual)
mse_optimal <- mse(preds$optimal, preds$actual)

# PLOT PREDICTIONS --------------------------------------------------------

pdf('pairs.pdf')
pairs(preds, pch = 20, cex = 0.5)
dev.off()
