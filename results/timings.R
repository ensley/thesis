#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)
library(tictoc)


eff_range <- function(dist, nu) {
  f <- function(alpha, h, nu) {
    matern_cor(h, nu, alpha, sigma = 1) - 0.05
  }
  uniroot(f, interval = c(1e-6, 100), h = dist, nu = nu)$root
}

nu <- 1.5
alpha <- eff_range(0.3, nu)
tau <- 0.01
nsamp <- c(5000, 10000, 20000, 50000, 100000)
N <- 100
K <- 15

df <- data.frame(nsamp = rep(nsamp, each = N), true = NA, est = NA)



tic.clearlog()
for (i in 1:length(nsamp)) {
  
  locations <- create_locations(400, 50)
  gp <- simulate_gp(locations, 'matern', c(nu, 1/alpha, 1, tau))
  
  knots <- seq(gp$m$knotlocs[1], gp$m$knotlocs[2], length = K)
  x <- exp(seq(gp$m$logx[1], gp$m$logx[2], length = 500))
  y <- gp$m$specdens(x)
  logx <- log(x)
  logy <- logx + log(y)
  w <- round(seq(1, length(x), length = 50))
  xx <- logx[w]
  yy <- logy[w]
  spl <- nsbasis(xx, knots)
  
  spl <- nsbasis(xx, knots)
  f <- fitspline(10000, yy, spl)
  beta_true <- colMeans(f$b[-(1:3000), ])
  rm(f)
  
  exact <- normal_ll_exact(gp$Y[gp$locs$type == 'obs'], gp$dist_obs, nu, alpha, sigma = 1, tau = tau)
  
  for (j in 1:N) {
    cat('i =', i, '\tj =', j, '\n')
    df$true[N*(i-1)+j] <- exact
    tic(i)
    df$est[N*(i-1)+j] <- likelihood(beta_true, knots, gp$Y[gp$locs$type == 'obs'], gp$dist_obs, B = nsamp[i], nugget = tau, debug = FALSE)
    toc(log = TRUE, quiet = TRUE)
  }
}
log_lst <- tic.log(format = FALSE)
timings <- unlist(lapply(log_lst, function(x) x$toc - x$tic))
names(timings) <- NULL
df$time <- timings


write.csv(df, 'timings.csv', row.names = FALSE)
