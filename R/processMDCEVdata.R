#' @title processMDCEVdata
#' @description Process MDCEV data
#' @inheritParams mdcev
#' @param alt_names name of alternatives.
#' @param model_options list of model options
#' @keywords internal
processMDCEVdata <- function(formula, data, model_options){

	formula <- Formula::Formula(formula)

	psi.vars <- stats::formula(formula, rhs = 1, lhs = 0)
	# Delete intercept term for psi variables.
	psi.vars <- stats::terms(psi.vars)
	attr(psi.vars, "intercept") <- 0
	dat_psi <- stats::model.matrix(psi.vars, data)
	NPsi_ij <- ncol(dat_psi)

	if (NPsi_ij == 0){
		dat_psi <- matrix(0, 0, NPsi_ij)
	}

	J <- nrow(unique(attr(data, "index")["alt"]))
	I <- nrow(unique(attr(data, "index")["id"]))

	if (model_options$model == "kt_ee"){
		if (is.null(model_options$psi_ascs)) model_options$psi_ascs = 0
		phi.vars <- stats::formula(formula, rhs = 3, lhs = 0)
		phi.vars <- stats::terms(phi.vars)
		attr(phi.vars, "intercept") <- 0
		dat_phi <- stats::model.matrix(phi.vars, data)
		NPhi <- ncol(dat_phi)
	} else if (model_options$model != "kt_ee"){
		if (is.null(model_options$psi_ascs)) model_options$psi_ascs = 1
		NPhi <- 0
		dat_phi <- matrix(0, 0, NPhi)
	}

	if (model_options$model == "gamma"){
		model_num <- 1
	} else if (model_options$model == "alpha"){
		model_num <- 2
	} else if (model_options$model == "hybrid"){
		model_num <- 3
	} else if (model_options$model == "hybrid0"){
		model_num <- 4
	} else if (model_options$model == "kt_ee"){
		model_num <- 5
	} else
		stop("No model specificied. Choose a model specification")

	# convert quant/price to matrices and income to vector
	price.name <- attr(data, "price")
	quant.name <- attr(data, "choice")
	income.name <- attr(data, "income")

	price <- matrix(data[[price.name]], ncol = J, byrow = TRUE)
	quant <- matrix(data[[quant.name]], ncol = J, byrow = TRUE)
	income <- as.vector(matrix(data[[income.name]], ncol = J, byrow = TRUE)[,1])

	# Put data into one list for rstan
	stan_data = list(I = I, J = J, NPsi_ij = NPsi_ij, NPhi = NPhi,
			 K = model_options$n_classes,
			 dat_psi = as.matrix(dat_psi),
			 dat_phi = as.matrix(dat_phi),
			 price_j = price,
			 quant_j = quant,
			 income = income,
			 model_num = model_num)
	stan_data = c(stan_data, model_options)
	stan_data$n_classes = NULL
	stan_data$model = NULL

	if (model_options$n_classes > 1){
		lc.vars <- formula(formula, rhs = 2, lhs = 0)

		data_class <- as_tibble(data) %>%
			dplyr::distinct(id, .keep_all = T) %>%
			stats::model.matrix(lc.vars, .)
		stan_data$data_class <- as.matrix(data_class)
		stan_data$L <- ncol(data_class) # number of membership variables
	}
return(stan_data)
}
