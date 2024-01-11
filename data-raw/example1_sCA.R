## code to prepare `example1_sCA` dataset goes here

load("death.2019.rda")
example1_sCA <- death.2019

usethis::use_data(example1_sCA, overwrite = TRUE)
