library(tidyverse)
directory <- '~/git/thesis/results'

files <- list.files(directory, pattern = "\\preds.csv$", recursive = TRUE, full.names = TRUE)
dfs <- lapply(files, function(f) {
  df <- read_csv(f,
           col_types = cols(
             X1 = col_skip(),
             actual = col_double(),
             spline = col_double(),
             bestmatern = col_double(),
             optimal = col_double(),
             hall = col_double()
             )
           )
  matches <- str_split(f, '/')[[1]]
  df$cov <- matches[7]
  df$method <- matches[8]
  df$dataset <- matches[9]
  df
})
allpreds <- bind_rows(dfs, .id = 'id')
mse <- allpreds %>% 
  group_by(method, cov, id = factor(as.numeric(id))) %>% 
  summarise_at(c('spline', 'bestmatern', 'optimal', 'hall'), funs(mean((. - actual)^2)))
mse_rel <- mse %>% 
  mutate_if(is.numeric, funs(. / optimal))

# mse_rel %>%
#   gather(type, mse, spline:hall) %>% 
#   filter(method == 'aao', mse <= 10) %>%
#   ggplot(aes(id, mse, fill = type)) +
#   geom_col(position = 'dodge') +
#   facet_grid(cov ~ .)

mse_rel %>%
  select(-optimal) %>% 
  gather(type, mse, spline, bestmatern, hall) %>%
  filter(type == 'spline') %>% 
  arrange(mse)

mse_rel %>% 
  select(-optimal) %>% 
  gather(type, mse, spline, bestmatern, hall) %>% 
  filter(method == 'aao', mse <= 100) %>% 
  ggplot(aes(type, mse)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, color = 'grey') +
  facet_grid(cov ~ .) +
  coord_flip() +
  ggthemes::theme_few()

mse_rel %>% ungroup() %>% count(cov)
