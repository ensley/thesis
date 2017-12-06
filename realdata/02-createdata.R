#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)
library(imager, quietly = TRUE)

args <- commandArgs(TRUE)

if (length(args) < 3) {
  stop('Not enough arguments to 02a')
}

filepath <- args[1]
# family <- args[2]
# M <- as.numeric(args[3])          # number of observations in the Gaussian process
# ngrid <- as.numeric(args[4])      # number of observations along one edge of prediction grid
K <- as.numeric(args[2])          # number of knots (and therefore number of betas)
B <- as.numeric(args[3])          # number of samples to take from the spectral density
# seed <- args[4]                   # will the random seed be set?
# 
# if (seed == 'y') {
#   set.seed(1234)
#   cat('SEED SET.\n')
# }



im <- load.image('~/git/thesis/film.png')
# demean
im <- im - mean(im)
N <- 10 # number of squares
M <- 25 # length of one side of each square
ncols <- dim(im)[1]
nrows <- dim(im)[2]
limit_start <- seq(1, by = M, length = N)
limit_end <- seq(M, by = M, length = N)

reps_list <- purrr::map2(limit_start, limit_end, ~ im[.x:.y,.y:.x,1,1])

purrr::iwalk(reps_list, function(v, i) {
  pdf(paste0('col_image_', i, '.pdf'))
  image(v, col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(20)))
  dev.off()
  pdf(paste0('bw_image_', i, '.pdf'))
  image(v, col = gray(seq(0, 1, length = 9)))
  dev.off()
})

pred_square <- im[50:74,180:156,1,1]
# pred_rows <- 18:21
# pred_cols <- 9:12
pred_rows <- 9:12
pred_cols <- 5:8
rmax <- max(pred_rows)
rmin <- min(pred_rows) - 1
cmax <- max(pred_cols)
cmin <- min(pred_cols) - 1

pdf('col_pred.pdf')
image(pred_square, col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(20)))
polygon((c(rmin, rmax, rmax, rmin) - 0.5)/24, 
        (c(cmin, cmin, cmax, cmax) - 0.5)/24, 
        border = 'red', lwd = 2)
dev.off()
pdf('bw_pred.pdf')
image(pred_square, col = gray(seq(0, 1, length = 20)))
polygon((c(rmin, rmax, rmax, rmin) - 0.5)/24, 
        (c(cmin, cmin, cmax, cmax) - 0.5)/24, 
        border = 'white', lwd = 2)
dev.off()

Y <- do.call('cbind', purrr::map(reps_list, as.numeric))

xgrid <- seq(0, 1, length = M)
grid <- expand.grid(x = xgrid, y = xgrid)
dist_obs <- as.matrix(dist(grid))


knots <- seq(-1, 2, length = K)
x <- exp(seq(-5, 3, length = 500))

beta_init <- c(-4.037633, 3.470971, -0.8238171, 1.340247, -2.783525, -0.2526895, -0.6673028, 1.366069, -1.687069, 1.228752, 3.486575, 1.057369, -1.554608, -1.636577, 1.288975)



tau_init <- 1
sigma.m <- 1e-4
v_tau <- 1
t_v <- 5
nugget <- 1e-2



cat('INITIAL VALUES: ', beta_init, '\n')



args <- list(dist_obs = dist_obs,
             Y = Y,
             B = B,
             knots = knots,
             beta_init = beta_init,
             tau_init = tau_init,
             sigma.m = sigma.m,
             v_tau = v_tau,
             t_v = t_v,
             nugget = nugget,
             x = x)


saveRDS(args, file.path(filepath, 'init_args.rds'))
print('initial data saved to init_args.rds')



subset_inds <- matrix(FALSE, nrow = M, ncol = M)
subset_inds[pred_rows,pred_cols] <- TRUE
subset_inds <- as.logical(subset_inds)

predYobs <- pred_square[!subset_inds]
predYpred <- pred_square[subset_inds]

pred_grid_obs <- grid[!subset_inds, ]
pred_grid_pred <- grid[subset_inds, ]

Npred <- length(predYpred)
Nobs <- length(predYobs)

pred2obs <- as.matrix(dist(rbind(pred_grid_obs, pred_grid_pred)))[Nobs + 1:Npred,1:Nobs]

pred_grid <- data.frame(cbind(rbind(pred_grid_obs, pred_grid_pred), 
                              type = c(rep('obs', Nobs), rep('pred', Npred))))


gpobj <- list(locs = pred_grid,
              dist_obs = as.matrix(dist(pred_grid_obs)),
              dist_pred = as.matrix(dist(pred_grid_pred)),
              mindist = NULL)
class(gpobj) <- 'GPlocations'
gpobj$m <- list(family = NULL,
                model = NULL,
                params = NULL,
                specdens = NULL,
                knotlocs = NULL,
                logx = NULL)
gpobj$Y <- predYobs
class(gpobj) <- 'GPsimulated'

saveRDS(gpobj, file.path(filepath, 'gpobj.rds'))
print('GPsimulated object containing film data saved to gpobj.rds')
