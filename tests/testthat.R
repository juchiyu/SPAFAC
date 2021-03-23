library(testthat)
# detach("package:SPAFAC", unload = TRUE)
# devtools::install_git(
#   "https://github.com/juchiyu/SPAFAC.git",
#   credentials = git2r::cred_user_pass("juchiyu", getPass::getPass())
# )

library(SPAFAC)

test_check("SPAFAC")
