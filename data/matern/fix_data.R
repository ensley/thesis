files <- list.files(pattern = '.rds$')
good_init <- c(-6.50585647, 1.32818210, -2.60928725, 3.13439640, -0.96105717, 0.75017782, -0.64266722, -0.68445846, 2.08036337 -2.56218706, 0.58193003, 1.82498580, 0.09131705, -1.60918892, 0.20776664)

for (f in files) {
  init_args <- readRDS(f)
  init_args$beta_init <- good_init
  saveRDS(init_args, f)
}