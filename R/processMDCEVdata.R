#' @title processMDCEVdata
#' @description Process MDCEV data
#' @inheritParams FitMDCEV
#' @param model_options list of model options
#' @keywords internal
processMDCEVdata <- function(data, psi_formula, lc_formula, model_options){

	dat_psi <- stats::model.matrix(psi_formula, data)

	NPsi <- ncol(dat_psi)

	J <- length(unique(data$good))
	I <- length(unique(data$id))

	if (model_options$model == "gamma"){
		model_num <- 1
	} else if (model_options$model == "alpha"){
		model_num <- 2
	} else if (model_options$model == "hybrid"){
		model_num <- 3
	} else if (model_options$model == "hybrid0"){
		model_num <- 4
	} else
		stop("No model specificied. Choose a model specification")

	# convert quant/price to matrices
	price <- matrix(data$price, ncol = J, byrow = TRUE)
	quant <- matrix(data$quant, ncol = J, byrow = TRUE)
	income <- matrix(data$income, ncol = J, byrow = TRUE)[,1]

	# Put data into one list for rstan
	stan_data =
		list(I = I, J = J, NPsi = NPsi,
			 K = model_options$n_classes,
			 dat_psi = as.matrix(dat_psi),
			 price_j = price,
			 quant_j = quant,
			 income = as.vector(income),
			 flat_priors = model_options$flat_priors,
			 prior_psi_sd = model_options$prior_psi_sd,
			 prior_gamma_sd = model_options$prior_gamma_sd,
			 prior_alpha_sd = model_options$prior_alpha_sd,
			 prior_scale_sd = model_options$prior_scale_sd,
			 prior_delta_sd = model_options$prior_delta_sd,
			 model_num = model_num,
			 fixed_scale1 = model_options$fixed_scale1,
			 trunc_data = model_options$trunc_data,
			 gamma_fixed = model_options$gamma_fixed,
			 alpha_fixed = model_options$alpha_fixed)

	if (model_options$n_classes > 1){
		data_class <- tbl_df(data) %>%
			dplyr::distinct(id, .keep_all = T) %>%
			stats::model.matrix(lc_formula, .)
		stan_data$data_class <- as.matrix(data_class)
		stan_data$L <- ncol(data_class) # number of membership variables
	}
return(stan_data)
}
