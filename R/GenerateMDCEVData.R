#' @title GenerateMDCEVData
#' @description Simulate data for MDCEV model
#' @inheritParams mdcev
#' @param nobs Number of individuals
#' @param ngoods Number of non-numeraire goods
#' @param inc_lo Low bound of income for uniform draw
#' @param inc_hi High bound of income for uniform draw
#' @param price_lo Low bound of price for uniform draw
#' @param price_hi High bound of price for uniform draw
#' @param alpha_parms Parameter value for alpha term
#' @param scale_parms Parameter value for scale term
#' @param gamma_parms Parameter value for gamma terms
#' @param psi_i_parms Parameter value for psi terms that vary by individual
#' @param psi_j_parms Parameter value for psi terms that vary by good
#' @param nerrs Number of error draws for demand simulation
#' @param tol Tolerance level for simulations if using general approach
#' @param max_loop maximum number of loops for simulations if using general approach
#' @return list with data for stan model and parms_true with parameter values
#' @export
#' @examples
#' \donttest{
#' data <- GenerateMDCEVData(model = "hybrid0")
#'}
GenerateMDCEVData <- function(model, nobs = 1000, ngoods = 10,
							  inc_lo = 100000, inc_hi = 150000,
							  price_lo = 100, price_hi = 500,
							  alpha_parms = 0.5,
							  scale_parms = 1,
							  gamma_parms = stats::runif(ngoods, 1, 2),
							  psi_i_parms = c(-1.5, 3, -2, 1, 2),
							  psi_j_parms = c(-5, 0.5, 2),
							  nerrs = 1,
							  tol = 1e-20,
							  max_loop = 999){

	income <-  stats::runif(nobs, inc_lo, inc_hi)
	price <- matrix(stats::runif(nobs*ngoods, price_lo, price_hi), nobs, ngoods)

	true <- gamma_parms
	parms <- c(paste(rep('gamma', ngoods), 1:ngoods, sep=""))
	gamma_true <- cbind(parms, true)

	true <- scale_parms
	parms <- 'scale'
	scale_true <- cbind(parms, true)

	# Create psi variables that vary over alternatives
	psi_j <- cbind(rep(1,ngoods), # add constant term
				   matrix(stats::runif(ngoods*(length(psi_j_parms)-1), 0 , 1), nrow = ngoods))
	psi_j <-  rep(1, nobs) %x% psi_j

	psi_i <- matrix(2 * stats::runif(nobs * length(psi_i_parms)), nobs,length(psi_i_parms))
	psi_i <- psi_i %x% rep(1, ngoods)

	dat_psi <- cbind(psi_j, psi_i)
	colnames(dat_psi) <- c(paste(rep('b', ncol(dat_psi)), 1:ncol(dat_psi), sep=""))

	true <- c(psi_j_parms, psi_i_parms)

	parms <- paste0(rep('psi', length(true)), sep="_",
					colnames(dat_psi))
	parms_true <- cbind(parms, true)

	if (model == "gamma"){
		model_num <- 1
		alpha_parms <- c(alpha_parms, rep(0, ngoods))
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "alpha"){
		model_num <- 2
		alpha_parms <- 0 + stats::runif(ngoods+1, 0.01, .98)
		gamma_parms <- rep(1, ngoods)
		true <- alpha_parms
		parms  <- c(paste(rep('alpha',ngoods),1:(ngoods+1),sep=""))
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "hybrid"){
		model_num <- 3
		alpha_parms <- rep(alpha_parms, ngoods+1)
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 0
	} else if (model == "hybrid0"){
		model_num <- 4
		alpha_parms <- rep(1e-3, ngoods+1)
		parms_true <- rbind(parms_true, gamma_true, scale_true)
		algo_gen <- 0
	} else
		stop("No model specificied. Choose a model")

	psi_parms <- c(psi_j_parms, psi_i_parms)

	psi_sims <- matrix(dat_psi %*% psi_parms, ncol = ngoods, byrow = TRUE)
	psi_sims <- CreateListsRow(psi_sims)
	psi_sims <- list(psi_sims )
	names(psi_sims) <- "psi_sims"

	income_list <- list(as.list(income))
	names(income_list) <- "income" # price normalized MU at zero

	price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
	price_list <- list(CreateListsRow(price_list))
	names(price_list) <- "price" # price normalized MU at zero

	df_indiv <- c(income_list, price_list, psi_sims)

	rstan::expose_stan_functions(stanmodels$SimulationFunctions)

	quant <- purrr::pmap(df_indiv, CalcmdemandOne_rng,
				  gamma_sim=gamma_parms,
				  alpha_sim=alpha_parms,
				  scale_sim=scale_parms,
				  nerrs=nerrs, algo_gen = algo_gen,
				  tol = tol, max_loop = max_loop)

	# Convert simulated data into data form for estimation
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)
	quant <- quant[,2:(ncol(quant))]
	quant <- as.vector(t(quant))
	price <- as.vector(t(price))

	id <- rep(1:nobs, each = ngoods)
	good <- rep(1:ngoods, times = nobs)
	income <- rep(income, each = ngoods)

	data <- as.data.frame(cbind(id, good, quant, price, dat_psi, income))

	out <- list(data = data,
				parms_true = parms_true)
	return(out)
}
