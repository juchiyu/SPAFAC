## code to prepare `example2_sDiSCA` dataset goes here
load("DisCA_EnFrAuthors.rda")

example2_sDiSCA <- FrEnAuthors.unblnc

usethis::use_data(example2_sDiSCA, overwrite = TRUE)

