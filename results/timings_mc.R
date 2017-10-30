#! /usr/bin --vanilla --default-packages=utils,gpcovr

library(gpcovr)
library(tictoc)

B <- 20000
M <- 4*10^(seq(1, by = 0.5, length = 6))
N <- 1

locations <- create_locations(M[2], 1, mindist = 1e-4)
gp <- simulate_gp(locations, 'matern', c(1.5, 0.1548, 1, 0.01))

knots <- seq(gp$m$knotlocs[1], gp$m$knotlocs[2], length = 15)
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

tic.clearlog()
for (i in 1:length(M)) {
  for (j in 1:N) {
    cat('i =', i, '\tj =', j, '\n')
    
    H <- as.numeric(dist(cbind(runif(M[i]), runif(M[i]))))
    wtilde <- rlogspline(B, beta_true, knots)
    
    tic(paste('exact', M[i], j, sep = '_'))
    cov_exact <- RandomFields::RFcovmatrix(gp$m$model, distances = H, dim = 2)
    toc(log = TRUE, quiet = TRUE)
    
    tic(paste('mc', M[i], j, sep = '_'))
    cov_mc <- gpumc::mc(H, wtilde)
    toc(log = TRUE, quiet = TRUE)
    
    rm(H, wtilde, cov_exact, cov_mc)
    gc()
    gc(reset = TRUE)
  }
  gc()
  gc(reset = TRUE)
}

log_lst <- tic.log(format = FALSE)


library(tidyverse)
library(ggthemes)
theme_set(theme_few(base_size = 18))
timings <- log_lst %>% 
  transpose() %>% 
  as_tibble() %>% 
  unnest() %>% 
  separate(msg, c('method', 'size', 'id'), sep = '_', convert = TRUE) %>% 
  mutate(elapsed = toc - tic) 

timings %>% 
  ggplot(aes(size, elapsed, color = method)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_log10() +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = prettyNum) +
  labs(x = 'Number of observations',
       y = 'Elapsed time (s)')

