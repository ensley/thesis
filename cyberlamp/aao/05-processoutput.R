#! /usr/bin --vanilla --default-packages=utils,openacc


library(openacc)
library(RColorBrewer)
library(RandomFields)
source('01a-functions-aao.R')


args <- commandArgs(TRUE)

if (length(args) < 2) {
  stop('Not enough arguments to 05a')
}

filepath <- args[1]
family <- args[2]

setwd(filepath)

plot.new()

args <- readRDS('init_args.rds')
b <- read.csv('allbetas.csv', header = FALSE)
N <- nrow(b)
# burnin <- floor(0.2 * N)
burnin <- 1
bf <- colMeans(tail(b, -burnin))
apply(b, 2, function(x) mean(diff(x) != 0))


plot(args$gp)

pdf('trace.pdf')
bmcmc <- coda::mcmc(b)
plot(bmcmc[ ,1:4])
plot(bmcmc[ ,2:5])
dev.off()

logx <- log(args$x)
logy <- logx + log(args$y)
w <- round(seq(1, length(args$x), length = 50))
xx <- logx[w]
yy <- logy[w]

spl <- nsbasis(logx, args$knots)
plot(logx, predict_natspl(nsbasis(logx, args$knots), bf), type = 'l')
lines(logx, predict_natspl(nsbasis(logx, args$knots), args$beta_true), col = 'blue')
legend('topleft', c('final estimate', '"Truth"'), col = c('black', 'blue'), lty = 1)
abline(v = args$knots, lty = 3, col = 'grey')

plot(logx, dloglogspline(logx, bf, args$knots), type = 'l')
lines(logx, dloglogspline(logx, args$beta_true, args$knots), col = 'blue')
legend('topleft', c('final estimate', '"Truth"'), col = c('black', 'blue'), lty = 1)
abline(v = args$knots, lty = 3, col = 'grey')


plot(logx, dlogspline(logx, args$beta_true, args$knots), type = 'l', col = 'blue')
lines(logx, dlogspline(logx, bf, args$knots))
legend('topleft', c('final estimate', '"Truth"'), col = c('black', 'blue'), lty = 1)
abline(v = args$knots, lty = 3, col = 'grey')


err_bounds <- get_error_bounds(tail(b, -burnin), logx, args$knots)
hgrid <- seq(1e-5, 5, length = 500)
cov_curves <- read.csv('errbounds.csv', header = FALSE)
cov_err_bounds <- apply(cov_curves, 2, function(co) quantile(co, c(0.025, 0.975)))


pdf('loglogdensity.pdf')
pts <- yy - log(normalize_dlogspline(nsbasis(logx, args$knots), args$beta_true, args$knots))
l1 <- dloglogspline(logx, args$beta_true, args$knots)
l2 <- dloglogspline(logx, bf, args$knots)
l3 <- dloglogspline(logx, args$beta_init, args$knots)
plot_ymin <- min(c(min(pts), min(err_bounds[ ,1]), min(err_bounds[ ,2]), min(l1), min(l2), min(l3)))
plot_ymax <- max(c(max(pts), max(err_bounds[ ,1]), max(err_bounds[ ,2]), max(l1), max(l2), max(l3)))
plot(range(logx), c(plot_ymin, plot_ymax), type = 'n',
        xlab = expression(paste(log, ' ', omega)),
        ylab = expression(paste(log, ' f(', log, ' ', omega, ')')),
        main = 'log-log density')
polygon(c(logx, rev(logx)), c(err_bounds[ ,1], rev(err_bounds[ ,2])), col = 'grey90', border = NA)
points(xx, pts, pch = 20, cex = 1.5, col = 'grey50')
lines(logx, l1, lwd = 2)
lines(logx, l2, col = 'blue', lwd = 2)
lines(logx, l3, lty = 2, lwd = 2, col = 'orange')
abline(v = args$knots, col = 'grey', lty = 3)
legend('topleft',
       c('true points', 'closest fit', 'estimated', 'initial'),
       col = c('grey50', 'black', 'blue', 'orange'),
       lty = c(NA, 1, 1, 2),
       lwd = c(NA, 2, 2, 2),
       pch = c(16, NA, NA, NA))
dev.off()



pdf('logdensity.pdf')
pts <- exp(yy)/normalize_dlogspline(nsbasis(logx, args$knots), args$beta_true, args$knots)
l1 <- dlogspline(logx, args$beta_true, args$knots)
l2 <- dlogspline(logx, bf, args$knots)
l3 <- dlogspline(logx, args$beta_init, args$knots)
plot_ymin <- min(c(min(pts), min(exp(err_bounds[ ,1])), min(exp(err_bounds[ ,2])), min(l1), min(l2), min(l3)))
plot_ymax <- max(c(max(pts), max(exp(err_bounds[ ,1])), max(exp(err_bounds[ ,2])), max(l1), max(l2), max(l3)))
plot(range(logx), c(plot_ymin, plot_ymax), type = 'n',
        xlab = expression(paste(log, ' ', omega)),
        ylab = expression(paste(' f(', log, ' ', omega, ')')),
        main = 'log density')
polygon(c(logx, rev(logx)), exp(c(err_bounds[ ,1], rev(err_bounds[ ,2]))), col = 'grey90', border = NA)
points(xx, exp(yy)/normalize_dlogspline(nsbasis(logx, args$knots), args$beta_true, args$knots), pch = 20, cex = 1.5, col = 'grey50')
lines(logx, dlogspline(logx, args$beta_true, args$knots), lwd = 2)
lines(logx, dlogspline(logx, bf, args$knots), col = 'blue', lwd = 2)
lines(logx, dlogspline(logx, args$beta_init, args$knots), lty = 2, lwd = 2, col = 'orange')
abline(v = args$knots, col = 'grey', lty = 3)
legend('topleft',
       c('true points', 'closest fit', 'estimated', 'initial'),
       col = c('grey50', 'black', 'blue', 'orange'),
       lty = c(NA, 1, 1, 2),
       lwd = c(NA, 2, 2, 2),
       pch = c(16, NA, NA, NA))
dev.off()

pdf('covar.pdf')
l1 <- cudatest::mc(hgrid, exp(rlogspline(50000, args$beta_true, args$knots)))
l2 <- cudatest::mc(hgrid, exp(rlogspline(50000, bf, args$knots)))
l3 <- RFcov(args$gp$m$model, hgrid)
plot_ymin <- min(c(min(cov_err_bounds[1, ]), min(cov_err_bounds[2, ]), min(l1), min(l2), min(l3)))
plot_ymax <- max(c(max(cov_err_bounds[1, ]), max(cov_err_bounds[2, ]), max(l1), max(l2), min(l3)))
plot(range(hgrid), c(plot_ymin, plot_ymax), type = 'n',
        xlab = 'h',
        ylab = 'C(h)',
        main = 'Covariance function')
polygon(c(hgrid, rev(hgrid)), c(cov_err_bounds[1, ], rev(cov_err_bounds[2, ])), col = 'grey90', border = NA)
lines(hgrid, l3, lwd = 2)
# lines(model_sim - RMnugget(var = 0.01), xlim = c(0, 5), lwd = 2)
lines(hgrid, l2, col = 'blue', lwd = 2)
abline(h = 0, lty = 3)
legend('topright',
       c('true', 'estimated'),
       col = c('black', 'blue'),
       lwd = 2)
dev.off()

#####



grid <- args$locations$locs
grid_pred <- grid[grid$type == 'pred',1:2]
grid_obs <- grid[grid$type == 'obs',1:2]
dist <- as.matrix(dist(grid[ ,1:2]))[grid$type == 'pred',grid$type == 'obs']
gamma1 <- matrix(cudatest::mc(dist, exp(rlogspline(50000, bf, args$knots))), nrow = length(args$gp$Y[grid$type == 'pred']))
G <- matrix(cudatest::mc(args$locations$dist_obs, exp(rlogspline(50000, bf, args$knots))), nrow = length(args$gp$Y[grid$type == 'obs']))
Ginv <- solve(G)

pred_ests <- drop(gamma1 %*% Ginv %*% args$gp$Y[grid$type == 'obs'])
mse <- mean((pred_ests - args$gp$Y[grid$type == 'pred'])^2)



# gamma1_best <- matrix(matern_cor(dist, 1.5, 4.76, 1), nrow = length(args$Ypred))
# G_best <- matrix(matern_cor(args$hobs, 1.5, 4.76, 1), nrow = length(args$Yobs))
gamma1_best <- matrix(RFcov(args$gp$m$model, as.numeric(dist)), nrow = length(args$gp$Y[grid$type == 'pred']))
G_best <- matrix(RFcov(args$gp$m$model, as.numeric(args$locations$dist_obs)), nrow = length(args$gp$Y[grid$type == 'obs']))
# diag(G_best) <- diag(G_best) + 0.01
Ginv_best <- solve(G_best)

pred_best <- drop(gamma1_best %*% Ginv_best %*% args$gp$Y[grid$type == 'obs'])
mse_best <- mean((pred_best - args$gp$Y[grid$type == 'pred'])^2)

( incr <- mse / mse_best - 1 )

write(incr, file = 'mse_perc_incr.txt')
