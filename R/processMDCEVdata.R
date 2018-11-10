#' @title processMDCEVdata
#' @description Process MDCEV data
#' @export

processMDCEVdata <- function(stan.dat,
							 dat_class = NULL, n_classes, price_num = NULL)
{
	NPsi <- length(stan.dat$dat_psi)
	dat_psi <- do.call(cbind, stan.dat$dat_psi)

	ngoods <- ncol(quant)
	nobs <- nrow(quant)

	if(is.null(price_num)) # default price numeraire is one
		price_num <- rep(1, nobs)

	exp_num <- stan.dat$inc - rowSums(stan.dat$price * stan.dat$quant)

	nonzero <- cbind(1, quant != 0)
	zero <- 1 - nonzero
	M <- rowSums(nonzero != 0)
	M_factorial = factorial(M-1)

	#------------------------------------#
	# Put data into one list for rstan
	dat_stan_j =
		list(I = nobs, J = ngoods, NPsi = NPsi,
			 dat_psi = as.matrix(dat_psi),
			 j_price = stan.dat$price,
			 j_quant = stan.dat$quant,
			 num_price = as.vector(price_num),
			 num_quant = as.vector(exp_num),
			 M = M,
			 M_factorial = M_factorial,
			 nonzero = nonzero,
			 zero = zero)

	if (n_classes > 1){

		if (is.null(dat_class)) {# default constant if no membership variables
			dat_stan_j$dat_class <- as.matrix(rep(1, nobs))
			dat_stan_j$L <- 1 # number of membership variables
		} else {
			dat_stan_j$dat_class <- dat_class
			dat_stan_j$L <- ncol(dat_class) # number of membership variables
		}
		dat_stan_j$K <- n_classes
	}
return(dat_stan_j)
}
