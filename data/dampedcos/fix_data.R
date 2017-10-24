files <- list.files(pattern = '.rds$')
good_init <- c(-4.0119576, 3.3824580, -0.7388719, 1.2512943, -2.6219876, -0.3484350, -0.6633029, 1.3271634, -1.5979129, 1.0988002, 3.5645099, 1.1255109, -1.7148295, -1.4804619, 1.2237576)

for (f in files) {
  init_args <- readRDS(f)
  init_args$beta_init <- good_init
  saveRDS(init_args, f)
}