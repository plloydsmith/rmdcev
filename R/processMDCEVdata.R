#' @title processMDCEVdata
#' @description Process MDCEV data
#' @export

processMDCEVdata <- function(data,
							 data_class = NULL,
							 price_num = NULL,
							 model_options)
{

	NPsi <- length(data$dat_psi)
	dat_psi <- do.call(cbind, data$dat_psi)

	J <- ncol(data$quant)
	I <- nrow(data$quant)

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

	if(is.null(price_num)) # default price numeraire is one
		price_num <- rep(1, I)

	nonzero <- cbind(1, data$quant != 0)
	M <- rowSums(nonzero != 0)
	M_factorial = factorial(M-1)

	# Put data into one list for rstan
	stan_data =
		list(I = I, J = J, NPsi = NPsi,
			 dat_psi = as.matrix(dat_psi),
			 j_price = data$price,
			 j_quant = data$quant,
			 num_price = as.vector(price_num),
			 income = as.vector(data$inc),
			 M_factorial = M_factorial,
			 model_num = model_num,
			 fixed_scale = model_options$fixed_scale,
			 n_parameters = n_parameters,
			 trunc_data = model_options$trunc_data)

	if (model_options$n_classes > 1){
		if (is.null(data_class)) {# default constant if no membership variables
			stan_data$data_class <- as.matrix(rep(1, I))
			stan_data$L <- 1 # number of membership variables
		} else {
			stan_data$data_class <- data_class
			stan_data$L <- ncol(data_class) # number of membership variables
		}
		stan_data$K <- model_options$n_classes
	}
return(stan_data)
}
