Sys.setenv("R_TESTS" = "")
library(testthat, rmdcev)
expose_stan_functions(stanmodels$SimulationFunctions)

test_check("rmdcev")
