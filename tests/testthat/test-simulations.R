context("Test Simulations")

tol <- 0.1
data(data_rec, package = "rmdcev")
data_rec

data_rec <- mdcev.data(data_rec, subset = id < 100,
				   id.var = "id",
				   alt.var = "alt",
				   choice = "quant")

result <- mdcev( ~ alt-1,
				   data = data_rec,
				   model = "hybrid0",
				 psi_ascs = 0,
				   algorithm = "MLE",
				 std_errors = "mvn",
				   print_iterations = FALSE,
				   backend = "rstan")

nalts <- result$stan_data[["J"]]
model_num <- result$stan_data[["model_num"]]
npols <- 2
policies<-	CreateBlankPolicies(npols, result, price_change_only = TRUE)

df_sim <- PrepareSimulationData(result, policies, nsims = 3)
income <- df_sim[["df_indiv"]][["income"]][[1]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[1]]
price <- df_sim[["df_indiv"]][["price"]][[1]]

price_j <- price[-1]


quant_num <- income - quant_j  %*% price[-1]
quant <- c(quant_num, quant_j)
psi_j <- result[["stan_data"]][["dat_psi"]][1:nalts,] %*% t(result[["stan_fit"]][["par"]][["psi"]])
phi_j <- rep(0,nalts)
gamma_j <- result[["stan_fit"]][["par"]][["gamma"]]
gamma <- c(1, gamma_j)
alpha <- rep(0, nalts+1)
#alpha <- rep(result[["stan_fit"]][["par"]][["alpha"]], nalts+1)
#alpha <- c(result[["stan_fit"]][["par"]][["alpha"]], rep(0,nalts))
scale <- result[["stan_fit"]][["par"]][["scale"]]

test_that("Conditional error hybrid0 draw", {
	skip_on_os("solaris")

tol_e <- 1e-20
tol_l <- 1e-20
max_loop = 999

PRNG <-rmdcev_get_rng(seed = 3)
o <- rmdcev_get_stream()

error <- DrawError_rng(quant_num, quant_j, price[-1],
			  psi_j, phi_j, gamma_j, alpha, scale,
			  model_num = model_num, nalts = nalts, nerrs = 2, cond_error = 1, draw_mlhs = 1,
			  PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err / price
#	 Test hybrid algo
	algo_gen = 0
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e, max_loop, o)

	expect_true(sum(abs(mdemand - quant)) < tol)
	# Test general algo
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen = 1, model_num, tol_e = tol_e, max_loop = max_loop, o)

	expect_true(sum(abs(mdemand - quant)) < tol)

error <- c(0.0000000000, -0.77612399, -0.55169780, -0.02143232, -0.81241994, -0.19768961,
		 -0.85621517, -0.82511555, -0.60241553, -1.14157948, -0.62773715, 0.69217576,
		0.61801650, 0.26789713, -0.66958419, -0.35102966, -0.90048354, 0.32263623)

psi_b_err <- exp(c(0, psi_j) + error)
MUzero_b <- psi_b_err / price
mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
							 nalts, algo_gen = 0, model_num, tol_e = tol_e, max_loop = max_loop, o)

	util <- ComputeUtilJ(income, mdemand[-1], price[-1],
							 psi_b_err, phi_j, gamma[-1], alpha,
							 nalts, model_num, o)

	expect_true(abs(util - log(income)) < tol)

	price_p <- price + c(.001,rep(1,nalts))
 	MUzero_p <- psi_b_err / price_p

	hdemand <- HicksianDemand(util, price_p, MUzero_p,  c(1, phi_j), gamma, alpha,
							nalts, algo_gen = 0, model_num, tol_l = tol_l, max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_true(abs(wtp_err - (-62.4995)) < tol)

})

test_that("Test demand simulation", {
	skip_on_os("solaris")

	# Test conditional errors
	demand <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err =1, nerrs = 3, sim_type = "demand")
	expect_true(abs(demand[[1]][[1]][1,1] - 62499.5) < .01)
	expect_true(abs(demand[[1]][[1]][1,2] - 0) < .01)

	# Test unconditional errors (currently returns -Inf for gamam0)
	#	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
	#						 cond_err = 0, nerrs = 3, sim_type = "welfare")
	#	sum_wtp <- SummaryWelfare(wtp)
	#	expect_true(sum(abs(sum_wtp$Mean)) < .01)

})


test_that("Test full simulation function", {

	# Test conditional errors
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err =1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < .01)

	# Test unconditional errors (currently returns -Inf for gamam0)
#	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
#						 cond_err = 0, nerrs = 3, sim_type = "welfare")
#	sum_wtp <- SummaryWelfare(wtp)
#	expect_true(sum(abs(sum_wtp$Mean)) < .01)

})


test_that("Conditional error hybrid draw", {
	skip_on_os("solaris")
	data_rec <- mdcev.data(data_rec, subset = id < 100,
						   id.var = "id",
						   alt.var = "alt",
						   choice = "quant")

	result <- mdcev( ~ alt-1,
					 data = data_rec,
					 model = "hybrid",
					 psi_ascs = 0,
					 algorithm = "MLE",
					 print_iterations = FALSE,
					 backend = "rstan")

	nalts <- result$stan_data[["J"]]
	model_num <- result$stan_data[["model_num"]]
	npols <- 2
	policies<-	CreateBlankPolicies(npols, result, price_change_only = TRUE)

	df_sim <- PrepareSimulationData(result, policies, nsims = 1)
	income <- df_sim[["df_indiv"]][["income"]][[2]]
	quant_j <- df_sim[["df_indiv"]][["quant_j"]][[2]]
	price <- df_sim[["df_indiv"]][["price"]][[2]]

	price_j <- price[-1]


	quant_num <- income - quant_j  %*% price[-1]
	quant <- c(quant_num, quant_j)
	psi_j <- result[["stan_data"]][["dat_psi"]][1:nalts,] %*% t(result[["stan_fit"]][["par"]][["psi"]])
	phi_j <- rep(0,nalts)
	gamma_j <- result[["stan_fit"]][["par"]][["gamma"]]
	gamma <- c(1, gamma_j)
	alpha <- rep(result[["stan_fit"]][["par"]][["alpha"]], nalts+1)
	scale <- result[["stan_fit"]][["par"]][["scale"]]
	print(scale, digits =10)
	expect_true(abs(scale - 0.6953965) < tol)

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop = 999

	PRNG <-rmdcev_get_rng(seed = 3)
	o <- rmdcev_get_stream()

	error <- DrawError_rng(quant_num, quant_j, price[-1],
						   psi_j, phi_j, gamma_j, alpha, scale,
						   model_num = model_num, nalts = nalts, nerrs = 2, cond_error = 1, draw_mlhs = 1,
						   PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err / price
	#	 Test hybrid algo
	algo_gen = 0
	mdemand_0 <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e, max_loop, o)

	expect_true(sum(abs(mdemand_0 - quant)) < tol)
	# Test general algo
	mdemand_1 <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen = 1, model_num, tol_e = tol_e, max_loop = max_loop, o)

	expect_true(sum(abs(mdemand_0 - mdemand_1)) < tol)

	error <- c(0.0000000000, 1.70860392,0.09740258,0.43580157,0.57046604,
			0.22554382, -0.95752753,0.79871627,0.64340736,0.38675191,
			-0.35739897,0.15933184, -0.40152362, -0.33500928,0.82087624,
			0.16345204,0.02210642, -0.04967792)

	psi_b_err <- exp(c(0, psi_j) + error)
	MUzero_b <- psi_b_err / price
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen = 0, model_num, tol_e = tol_e, max_loop = max_loop, o)

	util <- ComputeUtilJ(income, mdemand[-1], price[-1],
						 psi_b_err, phi_j, gamma[-1], alpha,
						 nalts, model_num, o)
	print(util, digits =10)
	expect_true(abs(util - 29.58674801) < tol)

	price_p <- price + c(.001,rep(1,nalts))
	MUzero_p <- psi_b_err / price_p

	hdemand <- HicksianDemand(util, price_p, MUzero_p,  c(1, phi_j), gamma, alpha,
							  nalts, algo_gen = 1, model_num, tol_l = tol_l, max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	print(wtp_err, digits =10)
	expect_true(abs(wtp_err - (-41.57857841)) < tol)

	# Test general form
	hdemand <- HicksianDemand(util, price_p, MUzero_p,  c(1, phi_j), gamma, alpha,
							  nalts, algo_gen = 0, model_num, tol_l = tol_l, max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	print(wtp_err, digits =10)
	expect_true(abs(wtp_err - (-41.57857841)) < tol)

})
