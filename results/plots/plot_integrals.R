library(tidyverse)
library(ggthemes)
theme_set(theme_few())

width = 12
height = 7


# HELPER FUNCTIONS --------------------------------------------------------

# boolean function to make sure a vector has a range of zero
# (i.e. all elements are identical)
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

# integrate via trapezoidal rule given vectors of x and y coordinates
integrate_trap <- function(x, y) {
  if (!zero_range(diff(x))) stop('irregular grid')
  delta <- x[2] - x[1]
  delta * sum((tail(y, -1) + head(y, -1)) / 2)
}

# labeller function for facet_grid
lblr <- function(value) {
  name <- list(
    'cov_bestmatern' = 'Best fit Matern',
    'cov_spline' = 'Spline fit',
    'cov_hall' = 'Hall method',
    'cov_hall_norm' = 'Hall method\n(normalized)',
    'matern' = 'Matern',
    'dampedcos' = 'Damped cosine'
  )
  name[value]
}


# FIND CSV FILES ----------------------------------------------------------


directory <- '~/git/thesis/results'
files <- list.files(directory, pattern = 'covariance.csv$', recursive = TRUE, full.names = TRUE)


# READ CSV FILES ----------------------------------------------------------


cov <- map_dfr(files, function(f) {
  df <- read_csv(f, col_types = cols(
    h = col_double(),
    cov_true = col_double(),
    cov_spline = col_double(),
    cov_bestmatern = col_double(),
    cov_hall = col_double(),
    cov_hall_norm = col_skip()
  ))
  matches <- stringr::str_split(f, '/')[[1]]
  df$family <- matches[7]
  df
}, .id = 'id')

n_covs <- 3

# CREATE VARIOUS DATA FRAMES AND CALCULATE INTEGRALS ----------------------


# cov with differences included
diffs <- cov %>% 
  mutate(id = as.numeric(id)) %>%
  mutate_at(vars(starts_with('cov_'), -cov_true), funs((cov_true - .)^2)) %>% 
  select(-cov_true)

# just integral values with id and family
integrals <- diffs %>% 
  group_by(family, id) %>% 
  summarise_at(vars(starts_with('cov_')), funs(integrate_trap(h, .))) %>% 
  gather(type, integral, -(family:id))

# first run (damped cosine, dataset001)
cov1 <- cov %>% 
  filter(id == 1)

# contains mean difference functions (for plot)
cov_meandiff <- diffs %>% 
  group_by(family, h) %>%
  summarise_at(vars(starts_with('cov_')), mean) %>% 
  gather(method, value, starts_with('cov_'))


# PLOTS -------------------------------------------------------------------


# boxplots of all integrals by family and method
integrals %>% 
  ggplot(aes(type, integral)) +
  geom_boxplot() +
  coord_flip() +
  facet_grid(family ~ .) +
  labs(title = '')
ggsave(file.path(directory, 'plots/boxplot_integrals.pdf'), width = width, height = height)

# true covariance vs estimates, one dataset
cov1 %>% 
  gather(method, value, starts_with('cov_')) %>%
  ggplot(aes(h, value, group = method)) +
  geom_line(aes(color = method, linetype = method), size = 1) +
  scale_color_manual(name = '', values = c(few_pal()(n_covs), 'gray')) +
  scale_linetype_manual(name = '', values = c(rep('solid', n_covs), 'dashed')) +
  labs(title = 'True vs Estimated Covariance Function',
       subtitle = 'Damped cosine, dataset 001',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/true_vs_est_one.pdf'), width = width, height = height)


cov %>% 
  filter(family == 'matern') %>%
  gather(method, value, starts_with('cov_')) %>% 
  ggplot(aes(h, value)) +
  geom_line(aes(group = interaction(id, method), color = method, size = method, alpha = method)) +
  scale_color_manual(name = '', values = c(few_pal()(n_covs), 'black')) +
  scale_size_manual(name = '', values = c(rep(0.5, n_covs), 1.2)) +
  scale_alpha_manual(name = '', values = c(rep(0.3, n_covs), 1)) +
  labs(title = 'True vs Estimated Covariance Function',
       subtitle = '100 Matern datasets',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/true_vs_est_matern.pdf'), width = width, height = height)

cov %>% 
  filter(family == 'dampedcos') %>%
  gather(method, value, starts_with('cov_')) %>% 
  ggplot(aes(h, value)) +
  geom_line(aes(group = interaction(id, method), color = method, size = method, alpha = method)) +
  scale_color_manual(name = '', values = c(few_pal()(n_covs), 'black')) +
  scale_size_manual(name = '', values = c(rep(0.5, n_covs), 1.2)) +
  scale_alpha_manual(name = '', values = c(rep(0.3, n_covs), 1)) +
  labs(title = 'True vs Estimated Covariance Function',
       subtitle = '100 damped cosine datasets',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/true_vs_est_dampedcos.pdf'), width = width, height = height)


cov %>% 
  gather(method, value, starts_with('cov_')) %>% 
  ggplot(aes(h, value)) +
  geom_line(aes(group = interaction(id, method), color = method, size = method, alpha = method)) +
  facet_grid(family ~ .) +
  scale_color_manual(name = '', values = c(few_pal()(n_covs), 'black')) +
  scale_size_manual(name = '', values = c(rep(0.5, n_covs), 1.2)) +
  scale_alpha_manual(name = '', values = c(rep(0.3, n_covs), 1)) +
  labs(title = 'True vs Estimated Covariance Function',
       subtitle = '100 datasets for each of Matern and damped cosine',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/true_vs_est_all.pdf'), width = width, height = height)


# squared differences, one dataset
diffs %>% 
  filter(id == 1) %>% 
  gather(method, value, starts_with('cov_')) %>% 
  ggplot(aes(h, value)) +
  geom_line(size = 1) +
  facet_grid(method ~ family, labeller = as_labeller(lblr)) +
  labs(title = 'Squared difference between true and estimated C(h)',
       subtitle = 'Damped cosine, dataset 001',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/sqdiff_one.pdf'), width = width, height = height)


# squared differences, all datasets (with mean)
diffs %>% 
  gather(method, value, starts_with('cov_')) %>% 
  ggplot(aes(h, value, group = id)) +
  geom_line(color = 'gray60', alpha = 0.5) +
  geom_line(data = cov_meandiff, aes(h, value, group = method), size = 1) +
  facet_grid(family ~ method, labeller = as_labeller(lblr)) +
  labs(title = 'Squared difference between true and estimated C(h)',
       subtitle = 'Comparing spline method to the best-fitting Matern model',
       x = 'h',
       y = 'C(h)')
ggsave(file.path(directory, 'plots/sqdiff.pdf'), width = width, height = height)

  