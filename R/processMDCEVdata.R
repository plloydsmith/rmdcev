#' @title .split_psi_for_rp
#' @description Split dat_psi into random and fixed submatrices when psi_random is specified.
#'   Updates NPsi_ij (random count), dat_psi (random columns), NPsi_ij_fixed (fixed count),
#'   and dat_psi_fixed (fixed columns) in place.
#' @noRd
.split_psi_for_rp <- function(stan_data, psi_random, data) {
	if (stan_data$NPsi_ij == 0L) return(stan_data)

	psi_rand_vars <- stats::terms(stats::formula(Formula::Formula(psi_random), rhs = 1, lhs = 0))
	attr(psi_rand_vars, "intercept") <- 0
	rand_cols <- colnames(stats::model.matrix(psi_rand_vars, data))
	if (is.null(rand_cols)) rand_cols <- character(0)

	all_cols <- colnames(stan_data$dat_psi)
	missing  <- setdiff(rand_cols, all_cols)
	if (length(missing) > 0)
		stop("psi_random refers to variables not in the psi formula: ",
		     paste(missing, collapse = ", "), call. = FALSE)

	fixed_cols <- setdiff(all_cols, rand_cols)
	if (length(fixed_cols) == 0L) return(stan_data)  # all terms already random

	I_J        <- nrow(stan_data$dat_psi)
	fixed_idx  <- match(fixed_cols, all_cols)
	rand_idx   <- match(rand_cols,  all_cols)

	stan_data$dat_psi_fixed <- stan_data$dat_psi[seq_len(I_J), fixed_idx, drop = FALSE]
	stan_data$NPsi_ij_fixed  <- length(fixed_idx)

	if (length(rand_idx) == 0L) {
		stan_data$dat_psi <- matrix(0, 0, 0)
		stan_data$NPsi_ij  <- 0L
	} else {
		stan_data$dat_psi <- stan_data$dat_psi[seq_len(I_J), rand_idx, drop = FALSE]
		stan_data$NPsi_ij  <- length(rand_idx)
	}
	stan_data
}

#' @title processMDCEVdata
#' @description Process MDCEV data
#' @inheritParams mdcev
#' @param alt_names name of alternatives.
#' @param model_options list of model options
#' @noRd
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
		if (is.null(model_options$psi_ascs)) model_options$psi_ascs <- 0
		phi.vars <- stats::formula(formula, rhs = 3, lhs = 0)
		phi.vars <- stats::terms(phi.vars)
		attr(phi.vars, "intercept") <- 0
		dat_phi <- stats::model.matrix(phi.vars, data)
		NPhi <- ncol(dat_phi)
	} else {
		if (is.null(model_options$psi_ascs)) model_options$psi_ascs <- 1
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
	} else if (model_options$model == "gamma1"){
		model_num <- 6
	} else
		stop("No model specificied. Choose a model specification")

	# convert quant/price to matrices and income to vector
	price.name <- attr(data, "price")
	id.name <- attr(data, "id")
	quant.name <- attr(data, "choice")
	income.name <- attr(data, "income")

	price <- matrix(data[[price.name]], ncol = J, byrow = TRUE)
	quant <- matrix(data[[quant.name]], ncol = J, byrow = TRUE)
	if (!is.null(income.name)) {
		income <- as.vector(matrix(data[[income.name]], ncol = J, byrow = TRUE)[,1])
	} else if (model_num == 6) {
		# gamma1 has no income effects; dummy income keeps quant_num = 1 in Stan,
		# which is multiplied by (alpha_0 - 1) = 0 so it never enters the likelihood or WTP
		income <- rowSums(price * quant) + 1
	} else {
		stop("income is required for model '", model_options$model,
			 "'. Provide an income column in mdcev.data().", call. = FALSE)
	}

	# Put data into one list for rstan
	stan_data <- list(I = I, J = J, NPsi_ij = NPsi_ij, NPhi = NPhi,
			 K = model_options$n_classes,
			 dat_psi = as.matrix(dat_psi),
			 dat_phi = as.matrix(dat_phi),
			 price_j = price,
			 quant_j = quant,
			 income = income,
			 model_num = model_num)

	# Copy exactly the Stan data-block keys from model_options (explicit to avoid
	# accidental name collisions if model_options gains new keys in the future).
	stan_keys <- c("fixed_scale1", "single_scale", "trunc_data", "psi_ascs",
	               "gamma_ascs", "jacobian_analytical_grad", "flat_priors",
	               "prior_psi_sd", "prior_phi_sd", "prior_gamma_sd",
	               "prior_alpha_shape", "prior_scale_sd", "prior_delta_sd",
	               "gamma_nonrandom", "alpha_nonrandom")
	for (k in stan_keys) stan_data[[k]] <- model_options[[k]]

	# Convert logical flags to integers for Stan (Stan data blocks use int, not logical).
	bool_keys <- c("fixed_scale1", "single_scale", "trunc_data", "gamma_ascs",
	               "jacobian_analytical_grad", "flat_priors",
	               "gamma_nonrandom", "alpha_nonrandom")
	for (k in bool_keys) {
		if (!is.null(stan_data[[k]]))
			stan_data[[k]] <- as.integer(stan_data[[k]])
	}
	if (!is.null(stan_data$psi_ascs))
		stan_data$psi_ascs <- as.integer(stan_data$psi_ascs)

	if (model_options$n_classes > 1){
		lc.vars <- formula(formula, rhs = 2, lhs = 0)

		data_class_df <- as_tibble(data) %>%
			dplyr::distinct(!!sym(id.name), .keep_all = TRUE)
		data_class <- stats::model.matrix(lc.vars, data_class_df)
		stan_data$data_class <- as.matrix(data_class)
		stan_data$L <- ncol(data_class) # number of membership variables
	}
	return(stan_data)
}
