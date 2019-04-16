#' @title processMDCEVdata
#' @description Process MDCEV data
#' @inheritParams FitMDCEV
#' @param model_options list of model options
#' @export
processMDCEVdata <- function(data, psi_formula, lc_formula,
							 num_price,
							 model_options){

	dat_psi <- stats::model.matrix(psi_formula, data)

	NPsi <- ncol(dat_psi)

	J <- length(unique(data$good))
	I <- length(unique(data$id))

	if (model_options$fixed_scale == 1)
		NScale <- 0
	else
		NScale <- 1

	if (model_options$model == "les"){
		model_num <- 1
		NAlpha <- 1
		NGamma <- J
	} else if (model_options$model == "alpha"){
		model_num <- 2
		NAlpha <- J + 1
		NGamma <- 0
	} else if (model_options$model == "gamma"){
		model_num <- 3
		NAlpha <- 1
		NGamma <- J
	} else if (model_options$model == "gamma0"){
		model_num <- 4
		NAlpha <- 0
		NGamma <- J
	} else
		stop("No model specificied. Choose a model specification")

	n_parameters <- NPsi + NGamma + NAlpha + NScale

	if(is.null(num_price)) # default price numeraire is one
		num_price <- rep(1, I)

	# convert quant/price to matrices
	price <- matrix(data$price, ncol = J, byrow = TRUE)
	quant <- matrix(data$quant, ncol = J, byrow = TRUE)
	inc <- matrix(data$inc, ncol = J, byrow = TRUE)[,1]

	nonzero <- cbind(1, quant != 0)
	M <- rowSums(nonzero != 0)
	M_factorial = factorial(M-1)

	# Put data into one list for rstan
	stan_data =
		list(I = I, J = J, NPsi = NPsi,
			 dat_psi = as.matrix(dat_psi),
			 j_price = price,
			 j_quant = quant,
			 num_price = as.vector(num_price),
			 income = as.vector(inc),
			 M_factorial = M_factorial,
			 no_priors = model_options$no_priors,
			 prior_psi = model_options$prior_psi,
			 prior_gamma = model_options$prior_gamma,
			 prior_alpha = model_options$prior_alpha,
			 prior_scale = model_options$prior_scale,
			 model_num = model_num,
			 fixed_scale = model_options$fixed_scale,
			 n_parameters = n_parameters,
			 trunc_data = model_options$trunc_data)

	if (model_options$n_classes > 1){
		stan_data$prior_beta_m <- model_options$prior_beta_m

		data_class <- tbl_df(data) %>%
			distinct(id, .keep_all = T) %>%
			stats::model.matrix(lc_formula, .)
		stan_data$data_class <- as.matrix(data_class)
		stan_data$n_parameters <- n_parameters * model_options$n_classes + ncol(data_class) * (model_options$n_classes - 1)
		stan_data$L <- ncol(data_class) # number of membership variables
		stan_data$K <- model_options$n_classes
	}
return(stan_data)
}
