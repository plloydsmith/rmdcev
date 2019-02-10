#' @title SimulateMdcevData
#' @description Simulate data for MDCEV model
#' @inheritParams FitMDCEV
#' @param nobs Number of individuals
#' @param ngoods Number of non-nuemraire goods
#' @param inc_lo Mean income
#' @param inc_hi Standard deviation of income
#' @param price_lo Mean price for non-numeraire goods
#' @param price_hi Standard deviation of price
#' @param alpha_parms Parameter value for alpha term
#' @param scale_parms Parameter value for scale term
#' @param gamma_parms Parameter value for gamma terms
#' @param psi_i_parms Parameter value for psi terms that vary by individual
#' @param psi_j_parms Parameter value for psi terms that vary by good
#' @param nerrs Number of error draws for demand simulation
#' @importFrom stats runif
#' @return data for stan model
SimulateMdcevData <- function(model, nobs = 1000, ngoods = 10,
							  inc_lo = 100000, inc_hi = 50000,
							  price_lo = 5, price_hi = 100,
							  alpha_parms = 0.8,
							  scale_parms = 1,
							  gamma_parms = 1 + runif(ngoods, 0, 2),
							  psi_i_parms = c(-1.5, 3, -2, 1, 2),
							  psi_j_parms = c(-5, 0.5, 2),
					 			nerrs = 15){

#nobs = 100
#ngoods = 10
#inc_mean = 100000
#inc_sd = 50000
#price_mean = 5
#price_sd = 100
#alpha_parms = 0.8
#scale_parms = 1
#gamma_parms = 1 + runif(ngoods, 0, 2)
#psi_i_parms = c(-1.5, 3, -2, 1, 2)
#psi_j_parms = c(-5, 0.5, 2)
#nerrs = 15
	inc <- inc_lo + runif(nobs, 0, inc_hi) # budget
	price <- price_lo + matrix(runif(nobs*ngoods, 0, price_hi), nobs, ngoods)  # price of non-numeraire good

#	psi_names <- c(paste(rep('psi_j',length(psi_j_parms)),1:length(psi_j_parms),sep=""),
#				   paste(rep('psi_i',length(psi_i_parms)),1:length(psi_i_parms),sep=""))
#	scale_names <- c('scale')
	true <- gamma_parms
	parms <- c(paste(rep('gamma', ngoods), 1:ngoods, sep=""))
	gamma_true <- cbind(parms, true)

	true <- scale_parms
	parms <- 'scale1'
	scale_true <- cbind(parms, true)
	# depending on model specification
#	if (model_type == 2){
#		alpha_names <- c(paste(rep('alpha',ngoods),1:(ngoods+1),sep=""))
#	} else if (model_type == 4) {
#		alpha_names <- NULL
#	} else alpha_names <- "alpha"

	# order is important
#	names <- c(psi_i_names, psi_j_names, gamma_names, alpha_names, scale_names)

	# Create psi variables that vary over alternatives
	psi_j <- cbind(rep(1,ngoods),matrix(runif(ngoods*(length(psi_j_parms)-1), 0 , 1), nrow = ngoods))
	psi_j <-  rep(1, nobs) %x% psi_j

	psi_i <- 2 * matrix(runif(nobs * length(psi_i_parms)), nobs,length(psi_i_parms))
	psi_i <- psi_i %x% rep(1, ngoods)

	dat_psi = cbind(psi_j, psi_i)

	colnames(dat_psi) <- c(paste(rep('psi', ncol(dat_psi)), 1:ncol(dat_psi), sep=""))

	true <- c(psi_j_parms, psi_i_parms)
	parms <- colnames(dat_psi)
	parms_true <- cbind(parms, true)

	if (model == "les"){
		model_num <- 1
		alpha_parms <- c(alpha_parms, rep(0, ngoods))
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "alpha"){
		model_num <- 2
		alpha_parms <- 0 + runif(ngoods+1, 0.01, .98)
		gamma_parms <- rep(1, ngoods)
		true <- alpha_parms
		parms  <- c(paste(rep('alpha',ngoods),1:(ngoods+1),sep=""))
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "gamma"){
		model_num <- 3
		alpha_parms <- rep(alpha_parms, ngoods+1)
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 0
	} else if (model == "gamma0"){
		model_num <- 4
		alpha_parms <- rep(1e-6, ngoods+1)
		parms_true <- rbind(parms_true, gamma_true, scale_true)
		algo_gen <- 0
	} else
		stop("No model specificied. Choose a model")

	psi_parms <- c(psi_j_parms, psi_i_parms)
	psi_sims <- matrix(dat_psi %*% psi_parms, ncol = ngoods, byrow = TRUE)

	psi_sims <- CreateListsRow(psi_sims)
	psi_sims <- list(psi_sims )
	names(psi_sims) <- "psi_sims"

	inc_list <- list(as.list(inc))
	names(inc_list) <- "inc" # price normalized MU at zero

	price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
	price_list <- list(CreateListsRow(price_list))
	names(price_list) <- "price" # price normalized MU at zero

	df_indiv <- c(inc_list, price_list, psi_sims)

	expose_stan_functions(rmdcev:::stanmodels$SimulationFunctions)

	quant <- pmap(df_indiv, CalcmdemandOne_rng,
				  gamma_sim=gamma_parms,
				  alpha_sim=alpha_parms,
				  scale_sim=scale_parms,
				  nerrs=nerrs, algo_gen = algo_gen)

	# Convert simulated data into estimation data
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)

	quant <- quant[,2:(ncol(quant))]
	quant <- as.vector(t(quant))
	price <- as.vector(t(price))

	id <- rep(1:nobs, each = ngoods)
	good <- rep(1:ngoods, times = nobs)
	inc <- rep(inc, each = ngoods) # budget

	data <- as.data.frame(cbind(id, good, quant, price, dat_psi, inc))

	out <- list(data = data,
				parms_true = parms_true)
return(out)
}
