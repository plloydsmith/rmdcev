#' @title PrepareSimulationData
#' @description Prepare data for WTP/demand simulation from a fitted mdcev object.
#' @param object An object of class \code{mdcev}.
#' @param policies A list produced by \code{\link{CreateBlankPolicies}} containing
#'   \code{price_p} (additive price changes) and optionally \code{dat_psi_p} /
#'   \code{dat_phi_p} (alternative-attribute changes).
#' @param nsims Number of parameter draws to use for uncertainty quantification.
#' @param class Class label for Latent Class models (e.g. \code{"class1"}).
#' @return A list with \code{df_indiv} (individual-level data), \code{df_common}
#'   (shared simulation inputs), and \code{sim_options} (model metadata).
#' @export
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#'
#' data_rec <- mdcev.data(data_rec, subset = id <= 500, id.var = "id",
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 0, data = data_rec,
#'                model = "hybrid0", algorithm = "MLE",
#'                std_errors = "mvn")
#'
#' policies <- CreateBlankPolicies(npols = 2, mdcev_est,
#'              price_change_only = TRUE)
#'
#' df_sim <- PrepareSimulationData(mdcev_est, policies)
#'
#'}
PrepareSimulationData <- function(object, policies, nsims = 30, class = "class1") {

	sd    <- object$stan_data
	npols <- length(policies$price_p)

	# Step 1: Extract and sample parameter draws.
	# nsims may be silently capped to the number of available draws.
	est_pars <- .extract_parameter_draws(object, nsims)
	nsims    <- nrow(est_pars)

	# Step 2: Pivot draws wide-to-long; filter to one latent class if needed.
	est_sim <- .pivot_draws_to_long(est_pars, object, class)

	# Step 3: Common simulation matrices shared across all individuals.
	scale_sims          <- .build_scale_sims(est_sim, sd$fixed_scale1, nsims)
	gamma_sim_nonrandom <- .build_gamma_sims(est_sim, sd$model_num, sd$gamma_nonrandom,
	                                          sd$gamma_ascs, nsims, sd$J)
	alpha_sim_nonrandom <- .build_alpha_sims(est_sim, sd$model_num, sd$alpha_nonrandom,
	                                          nsims, sd$J)

	# Step 4: Individual-level parameter draws (psi, phi, and random gamma/alpha).
	indiv <- .build_individual_pars(est_sim, est_pars, object, nsims)

	# Step 5: Baseline psi and phi utility indices (nsims x J per individual).
	psi_sims <- .build_psi_sims(indiv$psi_sim_temp, object, npols, nsims)
	phi_sims <- .build_phi_sims(indiv$phi_sim_temp, object, nsims)

	# Step 6: Policy-scenario psi/phi (zero-dimension matrices when price_change_only).
	psi_p_sims <- replicate(sd$I, matrix(0, 0, 0), simplify = FALSE)
	phi_p_sims <- replicate(sd$I, matrix(0, 0, 0), simplify = FALSE)
	if (!policies$price_change_only) {
		if (sd$model_num < 5)
			psi_p_sims <- .build_policy_psi_sims(indiv$psi_sim_temp, object, policies, npols)
		else
			phi_p_sims <- .build_policy_phi_sims(indiv$phi_sim_temp, object, policies, npols)
	}

	# Step 7: Assemble individual-level data (one element per individual, used by purrr::pmap).
	df_indiv <- c(
		list(income     = as.list(sd$income)),
		list(quant_j    = CreateListsRow(sd$quant_j)),
		list(price      = CreateListsRow(cbind(1, sd$price_j))),
		list(psi_sims   = psi_sims),
		list(phi_sims   = phi_sims),
		list(psi_p_sims = psi_p_sims),
		list(phi_p_sims = phi_p_sims),
		indiv$gamma_rand,   # NULL or list(gamma_sims = ...)
		indiv$alpha_rand    # NULL or list(alpha_sims = ...)
	)

	list(
		df_indiv   = df_indiv,
		df_common  = list(
			price_p_list        = policies$price_p,
			gamma_sim_nonrandom = gamma_sim_nonrandom,
			alpha_sim_nonrandom = alpha_sim_nonrandom,
			scale_sims          = scale_sims
		),
		sim_options = list(
			n_classes         = object$n_classes,
			model_num         = sd$model_num,
			price_change_only = policies$price_change_only
		)
	)
}
