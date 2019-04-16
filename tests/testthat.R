Sys.setenv("R_TESTS" = "")
library(testthat)
library(rmdcev)
rstan::expose_stan_functions(stanmodels$SimulationFunctions)

test_check("rmdcev")
