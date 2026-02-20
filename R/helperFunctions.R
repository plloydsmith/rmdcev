#' @title CreateListsRow
#' @description Convert matrix to a list with each row as an element.
#' @param x matrix to convert
#' @return A list
#' @noRd
CreateListsRow <- function(x) {
	lapply(seq_len(nrow(x)), function(i) x[i, ])
}

#' @title CleanInit
#' @description Reshape initial values to have the dimensions Stan expects.
#' @param init_input raw initial values (list or scalar)
#' @noRd
CleanInit <- function(init_input) {
	if (!is.list(init_input))
		return(init_input)
	temp <- lapply(init_input, function(x) matrix(x, nrow = 1, length(x)))
	if (!is.null(init_input$scale))
		temp$scale <- array(init_input$scale, dim = 1)
	temp
}

#' @title CreateListsCol
#' @description Convert matrix to a list with each column as an element.
#' @param x matrix (or vector) to convert
#' @noRd
CreateListsCol <- function(x) {
	if (is.vector(x)) as.list(x)
	else lapply(seq_len(ncol(x)), function(i) x[, i])
}

#' @title MultiplyMatrix
#' @description Matrix multiply mat_temp %*% x, returned reshaped to n_rows rows.
#' @noRd
MultiplyMatrix <- function(mat_temp, x, n_rows) {
	matrix(mat_temp %*% x, nrow = n_rows, byrow = TRUE)
}

#' @title DoCbind
#' @description Bind a list of vectors/matrices column-wise into one matrix.
#' @noRd
DoCbind <- function(x) {
	do.call(cbind, x)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
#' @param npols Number of policies to simulate
#' @param model Estimated model from mdcev
#' @param price_change_only Logical value for whether to include policy changes to dat_psi. Defaults to TRUE.
#' TRUE implies that only price changes are used in simulation.
#' @export
#' @examples
#' \donttest{
#' data_rec <- mdcev.data(data_rec, subset = id <= 500, id.var = "id",
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 0, data = data_rec,
#'                model = "hybrid0", algorithm = "MLE",
#'                std_errors = "mvn")
#' CreateBlankPolicies(npols = 2, mdcev_est)
#'}
CreateBlankPolicies <- function(npols, model, price_change_only = TRUE) {
	price_p   <- CreateListsRow(matrix(0, nrow = npols, ncol = model$stan_data$J + 1))
	model_num <- model$stan_data$model_num
	dat_psi_p <- NULL
	dat_phi_p <- NULL

	if (!price_change_only) {
		if (model_num < 5 && model$parms_info$n_vars$n_psi == 0)
			stop("No psi variables to vary! Use price_change_only == TRUE option.")
		if (model_num == 5 && model$parms_info$n_vars$n_phi == 0)
			stop("No phi variables to vary! Use price_change_only == TRUE option.")

		if (model_num < 5)
			dat_psi_p <- lapply(seq_len(npols), function(X) model$stan_data$dat_psi)
		else
			dat_phi_p <- lapply(seq_len(npols), function(X) model$stan_data$dat_phi)
	}

	list(price_p = price_p, dat_psi_p = dat_psi_p, dat_phi_p = dat_phi_p,
	     price_change_only = price_change_only)
}

#' @title CreatePsiMatrix
#' @param psi_j matrix (J x n_psi_j) of alt-specific attributes
#' @param psi_i matrix (I x n_psi_i) of person-specific attributes
#' @description Creates the Psi data matrix for use in mdcev model
#' @noRd
CreatePsiMatrix <- function(psi_j = NULL, psi_i = NULL) {
	if (!is.null(psi_i))
		psi_i <- lapply(psi_i, function(x) rep(x, each = nrow(psi_j)))
	if (!is.null(psi_j))
		psi_j <- lapply(psi_j, function(x) rep(x, times = nrow(psi_i)))
	c(psi_j, psi_i)
}

#' @title GrabParms
#' @param data long-format parameter draws (cols: sim_id, parms, value)
#' @param parm_name regex to match parameter names
#' @description Extract and pivot a subset of draws into a wide nsims x k matrix.
#' @noRd
GrabParms <- function(data, parm_name) {
	data %>%
		dplyr::filter(grepl(parm_name, .data$parms)) %>%
		tidyr::pivot_wider(id_cols = "sim_id", names_from = "parms", values_from = "value") %>%
		dplyr::select(-"sim_id") %>%
		as.matrix()
}

#' @title GrabIndividualParms
#' @param est_sim individual-level beta draws (cols: id, sim_id, parm_id, beta, parms)
#' @param parm_name regex to match parameter names
#' @description Extract individual-level draws for one parameter group; returns a list
#'   of I matrices (nsims x k), one per individual.
#' @noRd
GrabIndividualParms <- function(est_sim, parm_name) {
	est_sim %>%
		dplyr::filter(grepl(parm_name, .data$parms)) %>%
		dplyr::select("id", "sim_id", "parm_id", "beta") %>%
		tidyr::pivot_wider(names_from = "parm_id", values_from = "beta") %>%
		dplyr::select(-"sim_id") %>%
		dplyr::group_split(.data$id, .keep = FALSE)
}

#' @title extract_bayes_draws
#' @description Extract posterior draws as a data.frame from a fitted Bayes
#'   mdcev object, handling both the rstan and cmdstanr backends.
#'   Column names are returned in rstan dot-notation (e.g. \code{psi.1})
#'   so downstream code that was written for rstan works unchanged.
#' @param object An mdcev object with \code{algorithm == "Bayes"}.
#' @return A data.frame of posterior draws (rows = iterations, cols = parameters).
#' @noRd
extract_bayes_draws <- function(object) {
	if (isTRUE(object$backend == "rstan")) {
		if (!requireNamespace("rstan", quietly = TRUE))
			stop("rstan is required when backend = 'rstan'.")
		as.data.frame(rstan::extract(object$stan_fit, permuted = TRUE, inc_warmup = FALSE))
	} else {
		draws_df <- as.data.frame(posterior::as_draws_df(object$stan_fit$draws()))
		draws_df <- draws_df[, !names(draws_df) %in% c(".chain", ".iteration", ".draw"),
		                     drop = FALSE]
		# Normalize bracket notation (mu[1,2]) to rstan dot notation (mu.1.2)
		names(draws_df) <- gsub("\\[", ".", gsub(",", ".", gsub("\\]", "", names(draws_df))))
		draws_df
	}
}

#' @title get_bayes_summary
#' @description Return a tibble of posterior summary statistics for a fitted
#'   Bayes mdcev object, with columns \code{parms}, \code{n_eff}, \code{Rhat}
#'   (and others), working for both the rstan and cmdstanr backends.
#' @param object An mdcev object with \code{algorithm == "Bayes"}.
#' @return A tibble.
#' @noRd
get_bayes_summary <- function(object) {
	if (isTRUE(object$backend == "rstan")) {
		if (!requireNamespace("rstan", quietly = TRUE))
			stop("rstan is required when backend = 'rstan'.")
		summ_mat <- rstan::summary(object$stan_fit)$summary
		tibble::as_tibble(summ_mat) %>%
			dplyr::mutate(parms = row.names(summ_mat))
	} else {
		object$stan_fit$summary() %>%
			dplyr::rename(parms = "variable", n_eff = "ess_bulk", Rhat = "rhat")
	}
}

#' @title get_bayes_chain_info
#' @description Return a list with \code{chains}, \code{warmup}, and
#'   \code{total} post-warmup draws for a fitted Bayes mdcev object,
#'   working for both the rstan and cmdstanr backends.
#' @param object An mdcev object with \code{algorithm == "Bayes"}.
#' @return A named list.
#' @noRd
get_bayes_chain_info <- function(object) {
	if (isTRUE(object$backend == "rstan")) {
		if (!requireNamespace("rstan", quietly = TRUE))
			stop("rstan is required when backend = 'rstan'.")
		list(
			chains = object$stan_fit@sim[["chains"]],
			warmup = object$stan_fit@sim[["warmup"]],
			total  = object$stan_fit@sim[["chains"]] *
			         (object$stan_fit@sim[["iter"]] - object$stan_fit@sim[["warmup"]])
		)
	} else {
		meta    <- object$stan_fit$metadata()
		nchains <- object$stan_fit$num_chains()
		list(
			chains = nchains,
			warmup = meta$iter_warmup,
			total  = nchains * meta$iter_sampling
		)
	}
}

#' @title get_rstan_model
#' @description Return the pre-compiled rstan stanmodel object for a named model.
#' @param model_name Character string, either "mdcev" or "mdcev_rp".
#' @return A stanmodel object.
#' @noRd
get_rstan_model <- function(model_name) {
	if (!requireNamespace("rstan", quietly = TRUE))
		stop("Package 'rstan' is required for the rstan backend and hessian computation.")
	if (!is.list(stanmodels) || is.null(stanmodels[[model_name]]))
		stop("Stan model '", model_name, "' is not available. ",
		     "Ensure rstan is installed and the package was loaded correctly.")
	stanmodels[[model_name]]
}

#' @title rmdcev_get_rng
#' @description Create a boost::ecuyer1988 RNG seeded with \code{seed},
#'   returned as an external pointer suitable for passing to the compiled
#'   Stan simulation functions.  Defined here (not in RcppExports.R) so it
#'   survives \code{Rcpp::compileAttributes()} regeneration.
#' @noRd
rmdcev_get_rng <- function(seed = 0L) {
	if (!requireNamespace("rstan", quietly = TRUE))
		stop("Package 'rstan' is required for simulation. Please install it.")
	rstan::get_rng(seed = seed)
}

#' @title rmdcev_get_stream
#' @description Return an external pointer to \code{Rcpp::Rcout} for use as
#'   the \code{pstream__} argument in compiled Stan simulation functions.
#'   Defined here (not in RcppExports.R) so it survives
#'   \code{Rcpp::compileAttributes()} regeneration.
#' @noRd
rmdcev_get_stream <- function() {
	if (!requireNamespace("rstan", quietly = TRUE))
		stop("Package 'rstan' is required for simulation. Please install it.")
	rstan::get_stream()
}

#' @title CreatePsi
#' @param dat_vars_i psi covariate data for one individual (J rows, or npols*J rows
#'   in the policy case with a \code{policy} column)
#' @param est_pars_i psi parameter draws for that individual (nsims x n_psi matrix)
#' @description Compute the log-psi utility index (nsims x J matrix per individual),
#'   or a list of nsims matrices (each npols x J) in the multi-policy case.
#' @noRd
CreatePsi <- function(dat_vars_i, est_pars_i, J, NPsi_ij, psi_ascs, npols) {
	nsims <- nrow(est_pars_i)
	lpsi  <- matrix(0, nsims, J)

	# Add alternative-specific constants (ASCs occupy the first J-1 columns of est_pars_i)
	if (psi_ascs == 1) {
		asc_draws <- est_pars_i[, 1:(J - 1), drop = FALSE]
		lpsi <- lpsi + if (nsims == 1) c(0, asc_draws) else cbind(0, asc_draws)
	}

	# Columns of est_pars_i that correspond to covariate (non-ASC) parameters
	covar_cols <- if (psi_ascs == 1) seq.int(J, length.out = NPsi_ij) else seq_len(NPsi_ij)

	# Baseline case: one set of J alternatives
	if (NPsi_ij > 0 && nrow(dat_vars_i) == J) {
		lpsi <- lpsi + as.matrix(est_pars_i[, covar_cols, drop = FALSE]) %*%
		        t(as.matrix(dat_vars_i))
	}

	# Policy case: npols * J rows, with a `policy` grouping column
	if (nrow(dat_vars_i) > J) {
		if (NPsi_ij == 0) {
			# No covariates: replicate the baseline lpsi (with ASCs) across all policies
			lpsi <- lapply(CreateListsRow(lpsi), function(x)
				matrix(x, nrow = npols, ncol = J, byrow = TRUE))
		} else {
			# Covariates: compute lpsi separately for each policy, then rearrange
			covar_draws <- est_pars_i[, covar_cols, drop = FALSE]
			pol_data    <- dplyr::group_split(dat_vars_i, .data$policy, .keep = FALSE)
			lpsi_by_pol <- lapply(pol_data, function(d)
				lpsi + as.matrix(covar_draws) %*% t(as.matrix(d)))
			# lpsi_by_pol: list of npols matrices (nsims x J)
			# rearrange to a list of nsims matrices (each npols x J)
			lpsi_arr <- aperm(
				array(unlist(lpsi_by_pol), dim = c(nsims, J, npols)),
				perm = c(1, 3, 2)   # -> (nsims, npols, J)
			)
			lpsi <- lapply(seq_len(nsims), function(s) lpsi_arr[s, , ])
		}
	}
	lpsi
}

# ============================================================
# PrepareSimulationData helpers
# ============================================================

#' @title .extract_parameter_draws
#' @description Extract and sample parameter draws from a fitted mdcev object.
#'   Caps nsims to the number of available draws and drops diagnostic columns.
#' @noRd
.extract_parameter_draws <- function(object, nsims) {
	# Extract raw draws and determine the cap
	if (object$algorithm == "Bayes") {
		est_pars  <- extract_bayes_draws(object)
		max_draws <- nrow(est_pars)
	} else if (object$algorithm == "MLE" && object$std_errors == "mvn") {
		est_pars  <- tibble::as_tibble(object$stan_fit$theta_tilde)
		max_draws <- object$n_draws
	} else if (object$algorithm == "MLE" && object$std_errors == "deltamethod") {
		par      <- object$stan_fit$par
		par      <- par[!grepl("^log_like|^sum_log_lik|^tau_unif|^Sigma|^lp__|^theta",
		                       names(par))]
		est_pars  <- tibble::as_tibble(t(unlist(par)))
		max_draws <- 1L
	} else {
		stop("Unknown algorithm / std_errors combination in mdcev object.")
	}

	# Cap nsims and warn if needed
	if (nsims > max_draws) {
		warning("nsims capped to ", max_draws, " (available draws).")
		nsims <- max_draws
	}

	# Drop diagnostic columns and assign readable names for fixed-parameter models
	diag_pat <- "^log_like|^sum_log_lik|^tau_unif|^Sigma|^lp__"
	est_pars  <- est_pars[, !grepl(diag_pat, colnames(est_pars)), drop = FALSE]
	if (object$random_parameters == "fixed")
		names(est_pars)[seq_len(object$parms_info$n_vars$n_parms_total)] <-
			object$parms_info$parm_names$all_names

	est_pars[sample(nrow(est_pars), nsims), , drop = FALSE]
}

#' @title .pivot_draws_to_long
#' @description Pivot parameter draws from wide (cols = parameters) to long
#'   (cols: sim_id, parms, value). For Latent Class models, filters to the
#'   requested class.
#' @noRd
.pivot_draws_to_long <- function(est_pars, object, class) {
	est_sim <- est_pars %>%
		tibble::rowid_to_column("sim_id") %>%
		tidyr::pivot_longer(-"sim_id", names_to = "parms", values_to = "value")

	if (object$n_classes > 1) {
		if (object$stan_data$single_scale == 1)
			est_sim$parms[est_sim$parms == "scale"] <- paste0(class, ".scale")
		est_sim <- est_sim[grepl(class, est_sim$parms), ]
	}
	est_sim
}

#' @title .build_scale_sims
#' @description Return an nsims x 1 matrix of scale draws, or ones when scale is fixed.
#' @noRd
.build_scale_sims <- function(est_sim, fixed_scale1, nsims) {
	if (fixed_scale1 == 0) t(GrabParms(est_sim, "scale")) else matrix(1, nsims, 1)
}

#' @title .build_gamma_sims
#' @description Return an nsims x J matrix of non-random gamma draws, or NULL
#'   when gamma is individual-specific (random). model_num == 2 fixes gamma at 1.
#' @noRd
.build_gamma_sims <- function(est_sim, model_num, gamma_nonrandom, gamma_ascs, nsims, J) {
	if (model_num == 2) {
		matrix(1, nsims, J)
	} else if (gamma_nonrandom == 1) {
		g <- GrabParms(est_sim, "gamma")
		if (gamma_ascs == 0) g <- matrix(rep(g, each = J), ncol = J, byrow = TRUE)
		g
	} else {
		NULL   # random gamma handled per-individual in .build_individual_pars
	}
}

#' @title .build_alpha_sims
#' @description Return an nsims x (J+1) matrix of non-random alpha draws, or NULL
#'   when alpha is individual-specific. Column 1 = numeraire; 2:(J+1) = non-numeraire.
#'   model_num == 4 fixes all alpha at 0.
#' @noRd
.build_alpha_sims <- function(est_sim, model_num, alpha_nonrandom, nsims, J) {
	if (model_num == 4) {
		matrix(0, nsims, J + 1)
	} else if (alpha_nonrandom == 1) {
		a <- GrabParms(est_sim, "alpha")
		if (model_num == 1 || model_num == 5)
			a <- cbind(a, matrix(0, nsims, J))        # append zero alpha for non-numeraires
		else if (model_num == 3)
			a <- matrix(rep(a, each = J + 1), ncol = J + 1, byrow = TRUE)  # single alpha -> J+1
		a
	} else {
		NULL   # random alpha handled per-individual in .build_individual_pars
	}
}

#' @title .build_individual_pars
#' @description Dispatch to the fixed or random helper. Returns a list with:
#'   psi_sim_temp, phi_sim_temp, gamma_rand (NULL or list), alpha_rand (NULL or list).
#' @noRd
.build_individual_pars <- function(est_sim, est_pars, object, nsims) {
	if (object$random_parameters == "fixed")
		.indiv_pars_fixed(est_sim, object)
	else
		.indiv_pars_random(est_sim, est_pars, object, nsims)
}

#' @title .indiv_pars_fixed
#' @description For fixed-parameter models, replicate the pooled draw matrix for
#'   each individual (every individual shares the same parameter values per draw).
#' @noRd
.indiv_pars_fixed <- function(est_sim, object) {
	I <- object$stan_data$I
	list(
		psi_sim_temp = replicate(I, GrabParms(est_sim, "psi"), simplify = FALSE),
		phi_sim_temp = if (object$parms_info$n_vars$n_phi > 0)
			replicate(I, GrabParms(est_sim, "phi"), simplify = FALSE) else NULL,
		gamma_rand   = NULL,
		alpha_rand   = NULL
	)
}

#' @title .indiv_pars_random
#' @description For random-parameter models, compute individual betas via
#'   beta = mu + z * tau (with Cholesky scaling for correlated models),
#'   back-transform gamma (exp) and alpha (logit), and split into per-individual draws.
#' @noRd
.indiv_pars_random <- function(est_sim, est_pars, object, nsims) {
	J               <- object$stan_data$J
	model_num       <- object$stan_data$model_num
	gamma_nonrandom <- object$stan_data$gamma_nonrandom
	alpha_nonrandom <- object$stan_data$alpha_nonrandom

	# Step A: Build a (sim_id x parm_id) table of mu and tau draws
	est_sim_mu_tau <- est_sim %>%
		dplyr::filter(grepl("mu|tau", .data$parms)) %>%
		tidyr::separate("parms", into = c("parms", "parm_id"), sep = "\\.") %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id)) %>%
		tidyr::pivot_wider(names_from = "parms", values_from = "value") %>%
		dplyr::arrange(.data$sim_id)

	if (object$random_parameters == "corr")
		est_sim_mu_tau <- .compute_corr_tau(est_sim, est_sim_mu_tau, object)

	# Step B: Compute individual betas: beta_i = mu + z_i * tau, then name them
	est_sim <- est_sim %>%
		dplyr::filter(grepl("^z\\.", .data$parms)) %>%
		tidyr::separate("parms", into = c("parms", "id", "parm_id"), sep = "\\.") %>%
		tidyr::pivot_wider(names_from = "parms", values_from = "value") %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id), id = as.numeric(.data$id)) %>%
		dplyr::arrange(.data$id, .data$sim_id, .data$parm_id) %>%
		dplyr::left_join(est_sim_mu_tau, by = c("sim_id", "parm_id")) %>%
		dplyr::arrange(.data$id, .data$sim_id, .data$parm_id) %>%
		dplyr::mutate(
			beta  = .data$mu + .data$z * .data$tau,
			parms = rep(object$parms_info$parm_names$sd_names, nsims * object$n_individuals)
		) %>%
		dplyr::select(-"tau", -"mu", -"z")

	# Step C: Back-transform: gamma was estimated on log scale, alpha on logit scale
	est_sim <- dplyr::mutate(est_sim, beta = dplyr::case_when(
		gamma_nonrandom == 0 & grepl("gamma", .data$parms) ~ exp(.data$beta),
		alpha_nonrandom == 0 & grepl("alpha", .data$parms) ~ plogis(.data$beta),
		TRUE ~ .data$beta
	))

	# Step D: Extract random gamma and alpha draws, shaped for the C++ simulation functions
	gamma_rand <- if (gamma_nonrandom == 0) {
		g <- GrabIndividualParms(est_sim, "gamma")
		if (object$stan_data$gamma_ascs == 0)
			g <- lapply(g, function(x) matrix(rep(as.matrix(x), each = J), ncol = J, byrow = TRUE))
		list(gamma_sims = lapply(g, as.matrix))
	} else NULL

	alpha_rand <- if (alpha_nonrandom == 0) {
		a <- GrabIndividualParms(est_sim, "alpha")
		if (model_num == 1 || model_num == 5)
			a <- lapply(a, function(x) cbind(as.matrix(x), matrix(0, nrow(x), J)))
		else if (model_num == 3)
			a <- lapply(a, function(x) matrix(rep(as.matrix(x), each = J + 1), ncol = J + 1, byrow = TRUE))
		else
			a <- lapply(a, as.matrix)
		list(alpha_sims = a)
	} else NULL

	list(
		psi_sim_temp = lapply(GrabIndividualParms(est_sim, "psi"), as.matrix),
		phi_sim_temp = if (model_num == 5) GrabIndividualParms(est_sim, "phi") else NULL,
		gamma_rand   = gamma_rand,
		alpha_rand   = alpha_rand
	)
}

#' @title .compute_corr_tau
#' @description For correlated random-parameter models, scale each individual's tau
#'   by the Cholesky factor L_Omega (tau_scaled = tau %*% L_Omega) so that draws
#'   respect the estimated correlation structure.
#' @noRd
.compute_corr_tau <- function(est_sim, est_sim_mu_tau, object) {
	# Number of random parameters inferred from est_sim (avoids dependency on est_pars scope)
	num_rand <- length(grep("^mu\\.", unique(est_sim$parms)))

	# Split tau and L_Omega draws by sim_id for parallel mapply
	tau_by_sim    <- est_sim_mu_tau %>%
		dplyr::select("sim_id", "parm_id", "tau") %>%
		dplyr::group_split(.data$sim_id)
	lomega_by_sim <- est_sim %>%
		dplyr::filter(grepl("L_Omega", .data$parms)) %>%
		dplyr::arrange(.data$sim_id) %>%
		dplyr::group_split(.data$sim_id)
	sim_ids <- dplyr::distinct(dplyr::arrange(est_sim, .data$sim_id), .data$sim_id)

	# For each draw: scaled_tau (row vector) = tau (row vector) %*% L_Omega
	L <- mapply(
		function(tau_df, lomega_df) {
			l_omega <- matrix(lomega_df$value, nrow = num_rand, byrow = FALSE)
			as.vector(tau_df$tau %*% l_omega)
		},
		tau_by_sim, lomega_by_sim
	)
	L <- matrix(unlist(L), nrow = nrow(sim_ids), byrow = TRUE)
	colnames(L) <- paste0("parm_id", seq_len(num_rand))

	scaled_tau <- dplyr::bind_cols(sim_ids, tibble::as_tibble(L)) %>%
		tidyr::pivot_longer(-"sim_id", names_to = "parm_id", values_to = "tau") %>%
		dplyr::mutate(parm_id = as.numeric(gsub("[^0-9]", "", .data$parm_id))) %>%
		dplyr::arrange(.data$sim_id)

	# Return updated est_sim_mu_tau with scaled tau
	est_sim_mu_tau %>%
		dplyr::select("sim_id", "parm_id", "mu") %>%
		dplyr::left_join(scaled_tau, by = c("sim_id", "parm_id"))
}

#' @title .build_psi_sims
#' @description Build the baseline psi utility index for each individual.
#'   Returns a list of I matrices (nsims x J). When n_psi == 0, returns zero matrices.
#' @noRd
.build_psi_sims <- function(psi_sim_temp, object, npols, nsims) {
	sd <- object$stan_data
	if (object$parms_info$n_vars$n_psi > 0) {
		dat_vars <- tibble::tibble(id = rep(seq_len(sd$I), each = sd$J))
		if (sd$NPsi_ij > 0)
			dat_vars <- dplyr::bind_cols(dat_vars, tibble::as_tibble(sd$dat_psi))
		dat_vars <- dplyr::group_split(dat_vars, .data$id, .keep = FALSE)
		mapply(CreatePsi, dat_vars, psi_sim_temp,
		       J = sd$J, NPsi_ij = sd$NPsi_ij, psi_ascs = sd$psi_ascs, npols = npols,
		       SIMPLIFY = FALSE)
	} else {
		replicate(sd$I, matrix(0, nsims, sd$J), simplify = FALSE)
	}
}

#' @title .build_phi_sims
#' @description Build the baseline phi utility index for each individual.
#'   Returns a list of I matrices (nsims x J). When n_phi == 0, returns constant-one matrices.
#' @noRd
.build_phi_sims <- function(phi_sim_temp, object, nsims) {
	sd <- object$stan_data
	if (object$parms_info$n_vars$n_phi > 0) {
		dat_vars <- tibble::tibble(id = rep(seq_len(sd$I), each = sd$J))
		dat_vars <- dplyr::bind_cols(dat_vars, tibble::as_tibble(sd$dat_phi))
		dat_vars <- dplyr::group_split(dat_vars, .data$id, .keep = FALSE)
		# phi = exp(dat_phi %*% phi_draws'): J x nsims, then transposed to nsims x J
		mapply(function(x, y) t(exp(as.matrix(x) %*% t(as.matrix(y)))),
		       dat_vars, phi_sim_temp, SIMPLIFY = FALSE)
	} else {
		replicate(sd$I, array(1, dim = c(nsims, sd$J)), simplify = FALSE)
	}
}

#' @title .build_policy_psi_sims
#' @description Build psi utility indices for policy scenarios (model_num < 5).
#'   Returns a list of I elements; each is a list of npols matrices (nsims x J).
#' @noRd
.build_policy_psi_sims <- function(psi_sim_temp, object, policies, npols) {
	sd       <- object$stan_data
	dat_vars <- tibble::tibble(id = rep(seq_len(sd$I), each = sd$J))

	if (sd$NPsi_ij > 0) {
		# Stack policy-specific dat_psi with a policy index column, then bind individual ids
		pol_data <- mapply(cbind, policies[["dat_psi_p"]], "policy" = seq_len(npols),
		                   SIMPLIFY = FALSE)
		pol_data <- lapply(pol_data, function(x) dplyr::bind_cols(dat_vars, tibble::as_tibble(x)))
		dat_vars <- do.call(rbind, pol_data)
	}
	dat_vars <- dplyr::group_split(dat_vars, .data$id, .keep = FALSE)
	mapply(CreatePsi, dat_vars, psi_sim_temp,
	       J = sd$J, NPsi_ij = sd$NPsi_ij, psi_ascs = sd$psi_ascs, npols = npols,
	       SIMPLIFY = FALSE)
}

#' @title .build_policy_phi_sims
#' @description Build phi utility indices for policy scenarios (model_num == 5).
#'   Returns a list of I lists, each with npols matrices (nsims x J).
#' @noRd
.build_policy_phi_sims <- function(phi_sim_temp, object, policies, npols) {
	sd       <- object$stan_data
	dat_vars <- tibble::tibble(id = rep(seq_len(sd$I), each = sd$J))

	pol_data <- mapply(cbind, policies[["dat_phi_p"]], "policy" = seq_len(npols),
	                   SIMPLIFY = FALSE)
	pol_data <- lapply(pol_data, function(x) dplyr::bind_cols(dat_vars, tibble::as_tibble(x)))
	dat_vars <- dplyr::group_split(do.call(rbind, pol_data), .data$id, .keep = FALSE)

	mapply(
		function(x, y) {
			lapply(dplyr::group_split(x, .data$policy, .keep = FALSE), function(xx)
				t(exp(as.matrix(xx) %*% t(as.matrix(y)))))
		},
		dat_vars, phi_sim_temp,
		SIMPLIFY = FALSE
	)
}
