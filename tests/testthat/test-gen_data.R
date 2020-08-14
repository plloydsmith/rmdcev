context("Test Generate data")

#library(rmdcev)
tol <- 0.01

nobs <- 200
nalts <- 5
test_that("Generate gamma MDCEV data", {
	model = "gamma"
	sim.data <- GenerateMDCEVData(model = model, nobs = nobs, nalts = nalts)

	output <- mdcev(~ b1 + b2 + b3 + b4 + b5 + b6,
					   data = sim.data$data,
					   psi_ascs = 0,
					   model = model,
					   algorithm = "MLE",
					print_iterations = F)
	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 13)
})

test_that("Generate alpha MDCEV data", {
	model = "alpha"
	sim.data <- GenerateMDCEVData(model = model, psi_j = 0, nobs = nobs, nalts = nalts)

	output <- mdcev(~ b1 + b2 + b3,
					data = sim.data$data,
					psi_ascs = 0,
					model = model,
					algorithm = "MLE",
					print_iterations = F)
	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 10)
})

test_that("Generate hybrid0 MDCEV data", {
	model = "hybrid0"
	sim.data <- GenerateMDCEVData(model = model, psi_i_parms = 0, nobs = nobs, nalts = nalts)

	output <- mdcev(~ b1 + b2 + b3,
					data = sim.data$data,
					psi_ascs = 0,
					model = model,
					algorithm = "MLE",
					print_iterations = F)
	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 9)
})

test_that("Generate kt_ee MDCEV data", {
	model = "kt_ee"
	sim.data <- GenerateMDCEVData(model = model,
								  psi_i_parms = -1,
								  phi_parms = 0.5,
								  gamma_parms = c(2, 4, 6, 8, 10),
								  nobs = nobs, nalts = nalts)

	expect_equal(length(sim.data[["data"]][["quant"]]), nobs * nalts)
})
