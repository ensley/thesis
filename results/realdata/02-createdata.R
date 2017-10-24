#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)
library(imager, quietly = TRUE)

args <- commandArgs(TRUE)

if (length(args) < 4) {
  stop('Not enough arguments to 02a')
}

filepath <- args[1]
# family <- args[2]
# M <- as.numeric(args[3])          # number of observations in the Gaussian process
# ngrid <- as.numeric(args[4])      # number of observations along one edge of prediction grid
K <- as.numeric(args[2])          # number of knots (and therefore number of betas)
B <- as.numeric(args[3])          # number of samples to take from the spectral density
seed <- args[4]                   # will the random seed be set?

if (seed == 'y') {
  set.seed(1234)
  cat('SEED SET.\n')
}



# cat('CREATING LOCATIONS...')
# 
# 
# 
# locations <- create_locations(M, ngrid)
# 
# 
# 
# cat('DONE.\nSIMULATING OBSERVATIONS...')
# 
# 
# 
# gp <- simulate_gp(locations, family, params)
# 
# 
# 
# cat('DONE.\n')

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
  image(v, col = rev(RColorBrewer::brewer.pal(9, 'Blues')))
  dev.off()
  pdf(paste0('bw_image_', i, '.pdf'))
  image(v, col = gray(seq(0, 1, length = 9)))
  dev.off()
})

pred_square <- im[50:75,180:155,1,1]

pdf('col_pred.pdf')
image(pred_square, col = rev(RColorBrewer::brewer.pal(9, 'Blues')))
dev.off()
pdf('bw_pred.pdf')
image(pred_square, col = gray(seq(0, 1, length = 9)))
dev.off()

Y <- do.call('cbind', purrr::map(reps_list, as.numeric))

xgrid <- seq(0, 1, length = M)
grid <- expand.grid(x = xgrid, y = xgrid)
dist_obs <- as.matrix(dist(grid))


knots <- seq(-1, 2, length = K)
x <- exp(seq(-5, 3, length = 500))
# y <- gp$m$specdens(x)
# logx <- log(x)
# logy <- logx + log(y)
# w <- round(seq(1, length(x), length = 50))
# xx <- logx[w]
# yy <- logy[w]
# spl <- nsbasis(xx, knots)



# cat('FINDING TRUE BETA VECTOR...')
# 
# 
# 
# f <- fitspline(500000, yy, spl)
# beta_true <- colMeans(f$b)
# rm(f)



# cat('DONE.\nDETERMINING INITIAL BETA VECTOR...\n')


# 
# invalid <- TRUE
# while(invalid) {
#   beta_init <- beta_true + runif(K, -0.5, 0.5)
#   slopes <- get_slopes(beta_init, knots)
#   invalid <- slopes[1] <= 0 || slopes[2] >= 0
#   cat('\tANOTHER TRY\n')
# }
# 

# beta_init <- c(-4.0119576, 3.3824580, -0.7388719, 1.2512943, -2.6219876, -0.3484350, -0.6633029, 1.3271634, -1.5979129, 1.0988002, 3.5645099, 1.1255109, -1.7148295, -1.4804619, 1.2237576)
beta_init <- c(-4.037633, 3.470971, -0.8238171, 1.340247, -2.783525, -0.2526895, -0.6673028, 1.366069, -1.687069, 1.228752, 3.486575, 1.057369, -1.554608, -1.636577, 1.288975)


# cat('INITIAL BETA DETERMINED.\n')



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
