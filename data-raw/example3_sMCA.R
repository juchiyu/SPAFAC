## code to prepare `example3_sMCA` dataset goes here

IOP.cat <- readRDS("IOP.cat.rds")
Demo <- readRDS("Demo.rds")

example3_sMCA <- cbind(Demo$gender, IOP.cat)

usethis::use_data(example3_sMCA, overwrite = TRUE)
