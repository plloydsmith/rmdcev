#' @title GenerateMDCEVData
#' @description Simulate data for KT models
#' @inheritParams mdcev
#' @param nobs Number of individuals
#' @param nalts Number of non-numeraire alts
#' @param income Vector of individual income
#' @param price Matrix of prices for non-numeraire alternatives.
#' @param alpha_parms Parameter value for alpha term
#' @param scale_parms Parameter value for scale term
#' @param gamma_parms Parameter value for gamma terms
#' @param psi_i_parms Parameter value for psi terms that vary by individual
#' @param psi_j_parms Parameter value for psi terms that vary by alt (all models except kt_ee)
#' @param phi_parms Parameter value for phi terms that vary by alt (kt_ee model only)
#' @param dat_psi_i (nobs X # psi_i_parms) matrix with individual-specific characteristics
#' @param dat_psi_j (nalts X # psi_j_parms) matrix with alternative-specific variables (all models except kt_ee)
#' @param dat_phi (nalts X # phi_parms) matrix with alternative-specific variables (kt_ee model only)
#' @param nerrs Number of error draws for demand simulation
#' @param tol Tolerance level for simulations if using general approach
#' @param max_loop maximum number of loops for simulations if using general approach
#' @return A `mdcev.data` object, which is a `data.frame` in long
#'     format. Also includes parms_true with parameter values
#' @export
#' @examples
#' \donttest{
#' data <- GenerateMDCEVData(model = "gamma")
#'}
GenerateMDCEVData <- function(model, nobs = 1000, nalts = 10,
							  income = stats::runif(nobs, 100000, 150000),
							  price = matrix(stats::runif(nobs*nalts, 100, 500), nobs, nalts),
							  alpha_parms = 0.5,
							  scale_parms = 1,
							  gamma_parms = stats::runif(nalts, 1, 10),
							  psi_i_parms = c(-1.5, 2, -1),
							  psi_j_parms = c(-5, 0.5, 2),
							  phi_parms = c(-5, 0.5, 2),
							  dat_psi_i = matrix(2 * stats::runif(nobs * length(psi_i_parms)), nobs, length(psi_i_parms)),
							  dat_psi_j = cbind(matrix(stats::runif(nalts*(length(psi_j_parms)), 0 , 1), nrow = nalts)),
							  dat_phi = cbind(matrix(stats::runif(nalts*(length(phi_parms)), 0 , 1), nrow = nalts)),
							  nerrs = 1,
							  tol = 1e-20,
							  max_loop = 999){

	true <- gamma_parms
	parms <- c(paste(rep('gamma', nalts), 1:nalts, sep = ""))
	gamma_true <- cbind(parms, true)

	true <- scale_parms
	parms <- 'scale'
	scale_true <- cbind(parms, true)

	if (model != "kt_ee"){
	# Create psi variables that vary over alternatives
		dat_psi_j <-  rep(1, nobs) %x% dat_psi_j
		dat_phi <- NULL
		phi_sims <- replicate(nobs, rep(0, nalts), simplify = FALSE)
	} else {
		psi_j_parms <- NULL
		dat_psi_j <- NULL
		dat_phi <-  rep(1, nobs) %x% dat_phi
		colnames(dat_phi) <- c(paste(rep('phi_', ncol(dat_phi)), 1:ncol(dat_phi), sep=""))
		phi_sims <- matrix(dat_phi %*% phi_parms, ncol = nalts, byrow = TRUE)
		phi_sims <- CreateListsRow(phi_sims)
		parms <- paste0(rep('phi', length(phi_parms)), sep="_", colnames(dat_phi))
		true <- c(phi_parms)
		phi_true <- cbind(parms, true)

	}

	dat_psi_i <- dat_psi_i %x% rep(1, nalts)
	dat_psi <- cbind(dat_psi_j, dat_psi_i)
	colnames(dat_psi) <- c(paste(rep('b', ncol(dat_psi)), 1:ncol(dat_psi), sep=""))

	# Create blank phi

	# Name all parameters true
	true <- c(psi_j_parms, psi_i_parms)

	parms <- paste0(rep('psi', length(true)), sep="_", colnames(dat_psi))

	parms_true <- cbind(parms, true)

	if (model == "gamma"){
		model_num <- 1
		alpha_parms <- c(alpha_parms, rep(0, nalts))
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "alpha"){
		model_num <- 2
		alpha_parms <- 0 + stats::runif(nalts+1, 0.01, .98)
		gamma_parms <- rep(1, nalts)
		true <- alpha_parms
		parms  <- c(paste(rep('alpha',nalts),1:(nalts+1),sep=""))
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, alpha_true, scale_true)
		algo_gen <- 1
	} else if (model == "hybrid"){
		model_num <- 3
		alpha_parms <- rep(alpha_parms, nalts+1)
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 0
	} else if (model == "hybrid0"){
		model_num <- 4
		alpha_parms <- rep(0, nalts+1)
		parms_true <- rbind(parms_true, gamma_true, scale_true)
		algo_gen <- 0
	} else if (model == "kt_ee"){
		model_num <- 5

		true <- c(psi_i_parms)

		parms <- paste0(rep('psi', length(true)), sep="_", colnames(dat_psi))

		parms_true <- cbind(parms, true)

		alpha_parms <- c(alpha_parms, rep(0, nalts))
		true <- alpha_parms[1]
		parms <- 'alpha1'
		alpha_true <- cbind(parms, true)
		parms_true <- rbind(parms_true, phi_true, gamma_true, alpha_true, scale_true)
		algo_gen <- 1
	} else
		stop("No model specificied. Choose a model")

	psi_parms <- c(psi_j_parms, psi_i_parms)

	psi_sims <- matrix(dat_psi %*% psi_parms, ncol = nalts, byrow = TRUE)
	psi_sims <- CreateListsRow(psi_sims)
	psi_j <- list(psi_sims )
	names(psi_j) <- "psi_j"

	phi_j <- list(phi_sims)
	names(phi_j) <- "phi_j"

	income_list <- list(as.list(income))
	names(income_list) <- "income" # price normalized MU at zero

	price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
	price_list <- list(CreateListsRow(price_list))
	names(price_list) <- "price" # price normalized MU at zero

	df_indiv <- c(income_list, price_list, psi_j, phi_j)

	PRNG <- rmdcev_get_rng(seed = 3L)
	o <- rmdcev_get_stream()

	quant <- purrr::pmap(df_indiv, CalcmdemandOne_rng,
						 gamma_j=gamma_parms,
						 alpha=alpha_parms,
						 scale=scale_parms,
						 nerrs=nerrs,
						 model_num = model_num,
						 algo_gen = algo_gen,
						 tol = tol, max_loop = max_loop,
						 PRNG, o)

	# Convert simulated data into data form for estimation
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)
	quant <- quant[,2:(ncol(quant))]
	quant <- as.vector(t(quant))
	price <- as.vector(t(price))

	id <- rep(1:nobs, each = nalts)
	alt <- rep(1:nalts, times = nobs)
	income <- rep(income, each = nalts)

	data <- as.data.frame(cbind(id, alt, quant, price, dat_psi, dat_phi, income))

	data <- mdcev.data(data,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

	out <- list(data = data,
				parms_true = parms_true)
	return(out)
}
