#' @title GenerateMDCEVDataRP
#' @description Simulate random parameter data for MDCEV model
#' @inheritParams FitMDCEV
#' @inheritParams GenerateMDCEVData
#' @importFrom stats runif
#' @return list with data for stan model and parms_true with parameter values
#' @export
GenerateMDCEVDataRP <- function(model,
								nobs = 1000,
								ngoods = 10,
								inc_lo = 100000,
								inc_hi = 150000,
								price_lo = 100,
								price_hi = 500,
								alpha_parms = 0.5,
								scale_parms = 1,
								gamma_parms = runif(ngoods, 1, 10),
								psi_j_parms = c(-5, 0.5, 2),
								nerrs = 1,
								tol = 1e-20,
								max_loop = 999){

	inc <-  runif(nobs, inc_lo, inc_hi)
	price <- matrix(runif(nobs*ngoods, price_lo, price_hi), nobs, ngoods)

	num_psi <- length(psi_j_parms)
	RP <- num_psi + ngoods + 1

	if (model == "les"){
		model_num <- 1
		algo_gen <- 1
	} else if (model == "alpha"){
		model_num <- 2
		alpha_parms <- 0 + runif(ngoods+1, 0.01, .98)
		algo_gen <- 1
	} else if (model == "gamma"){
		model_num <- 3
		algo_gen <- 0
	} else if (model == "gamma0"){
		RP <- RP - 1 # subtract one for fixed alpha
		model_num <- 4
		algo_gen <- 0
		alpha_parms <- NULL
	} else
		stop("No model specificied. Choose a model")

	beta <- c(psi_j_parms, gamma_parms, alpha_parms)

	# Random correlation matrix to construct a covariance matrix of the individual betas
	tau <- rgamma(RP, 4, 2)
	Sigma <- diag(tau) %*% cor(matrix(rnorm(RP*(RP+1)), RP+1, RP)) %*% diag(tau)

	a <- c(rep(-Inf, num_psi), rep(0.01, RP-num_psi))
	b <- c(rep(Inf, num_psi), rep(Inf, RP-num_psi))

	beta_individual <- rtmvnorm(n=nobs, mean=beta, sigma=Sigma, lower=a, upper=b)

	#	beta_individual <- MASS::mvrnorm(nobs, beta, Sigma)

	indexes <- tibble(individual = rep(1:nobs, each = ngoods),
					  task = rep(1:nobs, each = ngoods),
					  row = 1:(nobs*ngoods)) %>%
		group_by(task) %>%
		summarise(task_individual = first(individual),
				  start = first(row),
				  end = last(row))

	# Create psi variables that vary over alternatives
	psi_j <- cbind(rep(1,ngoods), # add constant term
				   matrix(runif(ngoods*(num_psi-1), 0 , 1), nrow = ngoods))
	psi_j <-  rep(1, nobs) %x% psi_j

	psi_beta <- beta_individual[,1:num_psi] %x% rep(1, ngoods)

	psi_sims <- matrix(rowSums(psi_j * psi_beta), ncol = ngoods, byrow = TRUE)
	psi_sims <- CreateListsRow(psi_sims)
	psi_sims <- list(psi_sims )
	names(psi_sims) <- "psi_sims"


	# Create gamma variables that vary over alternatives
	#	dat_gamma <- rep(1, nobs) %x% diag(ngoods)
	#	gamma_sims <- matrix(rowSums(dat_gamma * gamma_beta), ncol = ngoods, byrow = TRUE)

	gamma_sims <- beta_individual[,(num_psi+1):(num_psi+ngoods)] #%x% rep(1, ngoods)

	if (model == "les"){
		alpha_sims <- cbind(beta_individual[,RP], matrix(0, nobs, ngoods))
	} else if (model == "alpha"){
		alpha_sims <- beta_individual[,(num_psi+ngoods+1):RP]
		gamma_sims <- matrix(1, nobs, ngoods)
	} else if (model == "gamma"){
		alpha_sims <- matrix(rep(beta_individual[,(num_psi+ngoods+1):RP],ngoods+1), nrow = nobs, ncol = ngoods+1, byrow = F)
	} else if (model == "gamma0"){
		alpha_sims <- matrix(1e-6, nobs, ngoods+1)
	} else
		stop("No model specificied. Choose a model")


	gamma_sims <- CreateListsRow(gamma_sims)
	gamma_sims <- list(gamma_sims )
	names(gamma_sims) <- "gamma_sims"

	alpha_sims <- CreateListsRow(alpha_sims)
	alpha_sims <- list(alpha_sims )
	names(alpha_sims) <- "alpha_sims"

	inc_list <- list(as.list(inc))
	names(inc_list) <- "inc" # price normalized MU at zero

	price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
	price_list <- list(CreateListsRow(price_list))
	names(price_list) <- "price" # price normalized MU at zero

	df_indiv <- c(inc_list, price_list, psi_sims, gamma_sims, alpha_sims)

	expose_stan_functions(stanmodels$SimulationFunctions)
	#	model_src <- stanc_builder("src/stan_files/SimulationFunctions.stan")
	#	expose_stan_functions(model_src)

	quant <- pmap(df_indiv, CalcmdemandOne_rng,
				  scale_sim=scale_parms,
				  nerrs=nerrs, algo_gen = algo_gen,
				  tol = tol, max_loop = max_loop)

	# Convert simulated data into estimation data
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)
	quant <- quant[,2:(ncol(quant))]
	quant <- as.vector(t(quant))
	price <- as.vector(t(price))

	id <- rep(1:nobs, each = ngoods)
	good <- rep(1:ngoods, times = nobs)
	inc <- rep(inc, each = ngoods)

	data <- as.data.frame(cbind(id, good, quant, price, psi_j, inc))

	parms_true <- df_indiv
	parms_true$inc <- parms_true$price <-  NULL

	out <- list(data = data,
				parms_true = beta_individual,
				Sigma = Sigma,
				beta = beta)
	return(out)
}

