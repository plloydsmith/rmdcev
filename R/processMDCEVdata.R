#' @title processMDCEVdata
#' @description Process MDCEV data
#' @export

processMDCEVdata <- function(dat,
							 dat_class = NULL,
							 n_classes,
							 price_num = NULL,
							 model_specification,
							 fixed_scale)
{

	NPsi <- length(dat$dat_psi)
	dat_psi <- do.call(cbind, dat$dat_psi)

	J <- ncol(dat$quant)
	I <- nrow(dat$quant)

	if (fixed_scale == 1)
		NScale <- 0
	else
		NScale <- 1



	if (model_specification == "les"){
		model_type <- 1
		NAlpha <- 1
		NGamma <- J
	} else if (model_specification == "alpha"){
		model_type <- 2
		NAlpha <- J + 1
		NGamma <- 0
	} else if (model_specification == "gamma"){
		model_type <- 3
		NAlpha <- 1
		NGamma <- J
	} else if (model_specification == "gamma0"){
		model_type <- 4
		NAlpha <- 0
		NGamma <- J
	} else
		stop("No model specificied. Choose a model_specification")

	n_parameters <- NPsi + NGamma + NAlpha + NScale

	if(is.null(price_num)) # default price numeraire is one
		price_num <- rep(1, I)

	exp_num <- dat$inc - rowSums(dat$price * dat$quant)

	nonzero <- cbind(1, dat$quant != 0)
	M <- rowSums(nonzero != 0)
	M_factorial = factorial(M-1)

	#------------------------------------#
	# Put data into one list for rstan
	stan_dat =
		list(I = I, J = J, NPsi = NPsi,
			 dat_psi = as.matrix(dat_psi),
			 j_price = dat$price,
			 j_quant = dat$quant,
			 num_price = as.vector(price_num),
			 num_quant = as.vector(exp_num),
			 M = M,
			 M_factorial = M_factorial,
			 nonzero = nonzero,
			 model_type = model_type,
			 fixed_scale = fixed_scale,
			 n_parameters = n_parameters)

	if (n_classes > 1){

		if (is.null(dat_class)) {# default constant if no membership variables
			stan_dat$dat_class <- as.matrix(rep(1, dat$I))
			stan_dat$L <- 1 # number of membership variables
		} else {
			stan_dat$dat_class <- dat_class
			stan_dat$L <- ncol(dat_class) # number of membership variables
		}
		stan_dat$K <- n_classes
	}
return(stan_dat)
}
