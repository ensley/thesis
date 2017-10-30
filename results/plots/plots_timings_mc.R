library(tidyverse)

timings <- read_csv('~/git/thesis/results/timings_mc.csv', col_names = TRUE)
names(timings) <- c('tic', 'toc', 'method', 'size', 'elapsed')

timings %>% 
  group_by(method, size) %>% 
  summarise(elapsed_avg = median(elapsed),
            elapsed_min = quantile(elapsed, 0.25),
            elapsed_max = quantile(elapsed, 0.75)) %>% 
  ggplot(aes(size, elapsed_avg, color = method)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = elapsed_min, ymax = elapsed_max), width = 0.1) +
  scale_x_log10() +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = prettyNum) +
  scale_color_few() +
  labs(title = 'Runtime of Exact Covariance Calculation vs MC Approximation',
       x = 'Number of observations',
       y = 'Elapsed time (s)')
ggsave('~/git/thesis/results/plots/timings_mc.pdf', width = 12, height = 7)