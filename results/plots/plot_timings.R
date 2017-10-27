library(tidyverse)

df <- read_csv('~/git/thesis/results/timings.csv')

df %>% 
  filter(time < 1) %>% 
  ggplot(aes(factor(nsamp), time)) +
  geom_boxplot()
ggsave('~/git/thesis/results/plots/time_boxplots.pdf')


df %>% 
  mutate(rel_est = - (est - true) / true) %>% 
  ggplot(aes(factor(nsamp), rel_est)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = 'grey70')
ggsave('~/git/thesis/results/plots/acc_boxplots.pdf')

df %>% 
  mutate(rel_est = - (est - true) / true) %>% 
  group_by(nsamp) %>% 
  summarise(time = mean(time),
            rel_est = mean(rel_est))
