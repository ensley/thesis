library(tidyverse)
library(ggthemes)
theme_set(theme_few())

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

integrate_trap <- function(x, y) {
  if (!zero_range(diff(x))) stop('irregular grid')
  delta <- x[2] - x[1]
  delta * sum((tail(y, -1) + head(y, -1)) / 2)
}


directory <- '~/git/thesis/results'
files <- list.files(directory, pattern = 'covariance.csv$', recursive = TRUE, full.names = TRUE)


cov <- map_dfr(files, read_csv, .id = 'id')
integrals <- cov %>% 
  mutate_if(is.character, as.numeric) %>% 
  group_by(id) %>% 
  mutate(splinediff = (cov_true - cov_spline)^2,
         bestdiff = (cov_true - cov_bestmatern)^2) %>% 
  summarise(spline = integrate_trap(h, splinediff),
            bestmatern = integrate_trap(h, bestdiff)) %>% 
  gather(type, integral, spline:bestmatern)

integrals %>% 
  ggplot(aes(type, log(integral))) +
  geom_boxplot() +
  coord_flip()

cov1 <- cov %>% filter(id == 1)
with(cov1, integrate_trap(h, (cov_true-cov_spline)^2))
with(cov1, integrate_trap(h, (cov_true-cov_bestmatern)^2))

cov %>% 
  mutate_if(is.character, as.numeric) %>% 
  group_by(id) %>% 
  mutate(splinediff = (cov_true - cov_spline)^2,
         bestdiff = (cov_true - cov_bestmatern)^2) %>% 
  gather(type, integral, splinediff:bestdiff) %>% 
  ggplot(aes(h, integral, group = id)) +
  geom_line(color = 'gray50', alpha = 0.5) +
  facet_grid(type ~ .) +
  scale_x_continuous(limits = c(0, 2))

plot(h, Ch, type = 'l')
lines(h, Ctrue, col = 'blue')

plot(h, Cbestfit / Cbestfit[1], type = 'l')
lines(h, Ctrue, col = 'blue')

plot(h, (Ch-Ctrue)^2, type = 'l')