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
err_bounds <- get_error_bounds(b, args$gp, args$knots, thin = 1000)
cov_err_bounds <- read.csv(file.path(outpath, 'errbounds.csv'), header = FALSE)
cov_err_bounds2 <- t(apply(cov_err_bounds, 2, function(x) stats::quantile(x, c(0.025, 0.975))))


N <- nrow(b)
# burnin <- floor(0.2 * N)
burnin <- 1
bf <- colMeans(tail(b, -burnin))
apply(b, 2, function(x) mean(diff(x) != 0))






# LOCATIONS PLOT ----------------------------------------------------------

plot.GPsimulated <- function(x, ..., bw = FALSE) {
  pred_idx <- which(x$locs$type == 'pred')
  grid_pred <- x$locs[pred_idx, ]
  grid_obs <- x$locs[-pred_idx, ]
  Y_pred <- x$Y[pred_idx]
  Y_obs <- x$Y[-pred_idx]
  pal <- ifelse(bw, 'Greys', 'Blues')
  graphics::plot.default(c(0, 1), c(0, 1), type = 'n', ...)
  graphics::image(unique(grid_pred$x),
                  unique(grid_pred$y),
                  matrix(Y_pred, nrow = length(unique(grid_pred$x)), byrow = F),
                  col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(20)),
                  xlab = '', ylab = '', add = TRUE)
  graphics::points(grid_obs$x, grid_obs$y, col = 'white', pch = 20)
  graphics::points(grid_obs$x, grid_obs$y, col = 'black', pch = 1)
}

pdf('locations.pdf')
plot(args$gp, xlab = '', ylab = '')
dev.off()

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
                   gpModel = args$gp$m,
                   beta_true = args$beta_true,
                   beta_init = args$beta_init,
                   err_bounds = err_bounds)
dev.off()

### END LOG-LOG DENSITY PLOT



# LOG DENSITY PLOT --------------------------------------------------------

logdensity_plot <- function(beta, knots, show_knots = TRUE, show_points = FALSE, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL)
{
  seq_length <- if (!is.null(err_bounds)) nrow(err_bounds) else 500
  if (!is.null(gpModel)) {
    x <- exp(seq(gpModel$logx[1], gpModel$logx[2], length = seq_length))
    y <- gpModel$specdens(x)
    logx <- log(x)
    logy <- logx + log(y)
    w <- round(seq(1, length(x), length = 50))
    xx <- logx[w]
    yy <- logy[w]
    
    spl <- nsbasis(logx, knots)
    pts <- exp(yy) / normalize_dlogspline(spl, beta_true, knots)
  } else {
    show_points <- FALSE
    logx <- seq(-5, 5, length = seq_length)
    pts <- NULL
  }
  
  l1 <- dlogspline(logx, beta, knots)
  l2 <- if (!is.null(beta_true)) dlogspline(logx, beta_true, knots) else NULL
  l3 <- if (!is.null(beta_init)) dlogspline(logx, beta_init, knots) else NULL
  
  if (!is.null(err_bounds)) {
    min_lower_bd <- min(exp(err_bounds[ ,1]))
    min_upper_bd <- min(exp(err_bounds[ ,2]))
    max_lower_bd <- max(exp(err_bounds[ ,1]))
    max_upper_bd <- max(exp(err_bounds[ ,2]))
  } else {
    min_lower_bd <- min_upper_bd <- max_lower_bd <- max_upper_bd <- NULL
  }
  
  plot_ymin <- min(suppressWarnings(min(pts)), min_lower_bd, min_upper_bd, min(l1), suppressWarnings(min(l2)), suppressWarnings(min(l3)))
  plot_ymax <- max(suppressWarnings(max(pts)), max_lower_bd, max_upper_bd, max(l1), suppressWarnings(max(l2)), suppressWarnings(max(l3)))
  
  graphics::plot(range(logx), c(plot_ymin, plot_ymax), type = 'n',
                 xlab = expression(paste(log, ' ', omega)),
                 ylab = expression(paste('f(', log, ' ', omega, ')')),
                 main = 'log density')
  if (!is.null(err_bounds)) {
    graphics::polygon(c(logx, rev(logx)),
                      exp(c(err_bounds[ ,1], rev(err_bounds[ ,2]))),
                      col = 'grey90',
                      border = NA)
  }
  if (show_points) graphics::points(xx, pts, pch = 20, cex = 1.5, col = 'grey50')
  graphics::lines(logx, l1, lwd = 2, col = 'blue')
  # if (!is.null(beta_true)) graphics::lines(logx, l2, lty = 2, lwd = 2)
  # if (!is.null(beta_init)) graphics::lines(logx, l3, lty = 2, lwd = 2, col = 'orange')
  if (show_knots) graphics::abline(v = knots, col = 'grey', lty = 3)
  graphics::legend('topleft',
                   c('true points', 'estimated'),
                   col = c('grey50', 'blue'),
                   lty = c(NA, 1),
                   lwd = c(NA, 2),
                   pch = c(16, NA))
  
  invisible(gpModel)
}

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

covar_plot <- function(h, beta, knots, gpModel = NULL, beta_true = NULL, beta_init = NULL, err_bounds = NULL, ...)
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
                 main = 'Covariance function', ...)
  
  if (!is.null(err_bounds)) graphics::polygon(c(h, rev(h)), c(err_bounds[ ,1], rev(err_bounds[ ,2])), col = 'grey90', border = NA)
  graphics::lines(h, l1, lwd = 2, col = 'blue')
  if (!is.null(beta_true)) graphics::lines(h, l3, lwd = 2, lty = 2)
  # if (!is.null(gpModel)) graphics::lines(h, l3, col = 'orange', lwd = 2, lty = 2)
  graphics::abline(h = 0, lty = 3)
  graphics::legend('topright', c('estimated', 'true'), col = c('blue', 'black'), lwd = 2, lty = c(1, 2))
  
  invisible(gpModel)
}

pdf('covariance.pdf')
covar_plot(x,
           bf,
           args$knots,
           gpModel = args$gp$m,
           beta_true = args$beta_true,
           err_bounds = cov_err_bounds2,
           xlim = c(0, 2))
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
