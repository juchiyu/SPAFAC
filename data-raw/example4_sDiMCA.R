## code to prepare `example4_sDiMCA` dataset goes here

macs <- readxl::read_excel("MACS-SpaFac-Data.xlsx")

example4_sDiMCA <- macs

usethis::use_data(example4_sDiMCA, overwrite = TRUE)
