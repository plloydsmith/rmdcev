# Unit tests for helper functions in helperFunctions.R and PrepareSimulationData.R

# ── CreateListsRow / CreateListsCol round-trip ────────────────────────────

test_that("CreateListsRow splits matrix into rows", {
	m <- matrix(1:6, nrow = 3, ncol = 2)
	lst <- rmdcev:::CreateListsRow(m)
	expect_equal(length(lst), 3)
	expect_equal(lst[[1]], m[1, ])
	expect_equal(lst[[3]], m[3, ])
})

test_that("CreateListsCol splits matrix into columns", {
	m <- matrix(1:6, nrow = 3, ncol = 2)
	lst <- rmdcev:::CreateListsCol(m)
	expect_equal(length(lst), 2)
	expect_equal(lst[[1]], m[, 1])
	expect_equal(lst[[2]], m[, 2])
})

test_that("CreateListsCol handles vector input", {
	v <- c(10, 20, 30)
	lst <- rmdcev:::CreateListsCol(v)
	expect_equal(length(lst), 3)
	expect_equal(lst[[2]], 20)
})

# ── GrabParms ─────────────────────────────────────────────────────────────

test_that("GrabParms returns correct dimensions", {
	# Build a minimal est_sim tibble (long format)
	est_sim <- tibble::tibble(
		sim_id = rep(1:3, each = 2),
		parms  = rep(c("psi.1", "psi.2"), 3),
		value  = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
	)
	result <- rmdcev:::GrabParms(est_sim, "psi")
	expect_equal(nrow(result), 3)   # nsims
	expect_equal(ncol(result), 2)   # n_psi
})

# ── CreateBlankPolicies error paths ──────────────────────────────────────

test_that("CreateBlankPolicies errors when price_change_only=FALSE and n_psi==0", {
	data(data_rec, package = "rmdcev")
	data_small <- mdcev.data(data_rec, subset = id <= 50,
							 id.var = "id", alt.var = "alt", choice = "quant")

	# Fit a model with psi_ascs=0 and no formula covariates -> n_psi == 0
	output <- mdcev(~ 0,
					data = data_small,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	expect_error(
		CreateBlankPolicies(npols = 1, output, price_change_only = FALSE),
		regexp = "psi|phi|price_change_only"
	)
})

# ── .extract_parameter_draws nsims-cap warning ────────────────────────────

test_that(".extract_parameter_draws warns when nsims > available draws", {
	data(data_rec, package = "rmdcev")
	data_small <- mdcev.data(data_rec, subset = id <= 50,
							 id.var = "id", alt.var = "alt", choice = "quant")

	# n_draws = 5 -> theta_tilde has 5 rows
	output <- mdcev(~ 0,
					data = data_small,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					std_errors = "mvn",
					n_draws = 5,
					print_iterations = FALSE,
					backend = "rstan")

	policies <- CreateBlankPolicies(1, output, price_change_only = TRUE)

	# nsims = 100 > n_draws = 5 -> should warn
	expect_warning(
		PrepareSimulationData(output, policies, nsims = 100),
		regexp = "capped|nsims"
	)
})
