## code to prepare `example1_sCA` dataset goes here

load("death.2019.rda")
example1_sCA <- death.2019[-c(3, 8, 15, 19), ]
example1_sCA <- t(as.matrix(example1_sCA))

usethis::use_data(example1_sCA, overwrite = TRUE)
