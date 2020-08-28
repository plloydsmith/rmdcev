#' @title PrepareSimulationData
#' @description Prepare Data for WTP simulation
#' @param object an object of class `mdcev`
#' @param policies list containing
#' price_p with additive price increases, and
#' dat_psi_p with new psi data
#' @param nsims Number of simulation draws to use for parameter uncertainty
#' @param class The class number for Latent Class models.
#' @return A list with individual-specific data (df_indiv) and common data (df_common)
#' and n_classes for number of classes and model_num for model type
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
PrepareSimulationData <- function(object,
								  policies,
								  nsims = 30,
								  class = "class1"){

	# Get parameter estimates in matrix form
	if (object$algorithm == "Bayes") {
		est_pars <- as.data.frame(rstan::extract(object$stan_fit, permuted = TRUE, inc_warmup = FALSE))
	} else if (object$algorithm == "MLE") {
		if(object$std_errors == "mvn") {
			est_pars <- as_tibble(object[["stan_fit"]][["theta_tilde"]])
			if (nsims > object$n_draws){
				nsims <- object$n_draws
				warning("Number of simulations > Number of mvn draws from mdcev object. nsims has been set to: ", nsims)
			}
		} else if (object$std_errors == "deltamethod"){
			est_pars <- object[["stan_fit"]][["par"]]
			est_pars <- est_pars[!grepl(c("^log_like|^sum_log_lik|^tau_unif|^Sigma|^lp__|^theta"),names(est_pars))]
			est_pars <- as_tibble(t(unlist(est_pars)))
			if (nsims > 1){
				nsims <- 1
				warning("Number of simulations > 1. Use mvn draws to incorporate parameter uncertainty nsims has been set to: ", nsims)
			}
		}
	}

	# drop extra variables
	est_pars <- est_pars[,!grepl(c("^log_like|^sum_log_lik|^tau_unif|^Sigma|^lp__"), colnames(est_pars))]

	# Rename variables
	if (object$random_parameters == "fixed"){
		names(est_pars)[1:object$parms_info$n_vars$n_parms_total] <- object$parms_info$parm_names$all_names
	}

	# Checks on simulation options
	if (object$algorithm == "Bayes" && nsims > nrow(est_pars)) {
		nsims <- nrow(est_pars)
		warning("Number of simulations > Number of posterior draws from mdcev object. nsims has been set to: ", nsims)
	}

	# Sample from parameter estimate draws
	est_pars <- est_pars[sample(nrow(est_pars), nsims), ]

	est_sim <- est_pars  %>%
		tibble::rowid_to_column("sim_id") %>%
		tidyr::pivot_longer(-sim_id, names_to = "parms", values_to = "value")

#object <- mdcev_corr
#	if (object$random_parameters == "corr")
#		stop("Demand and welfare simulation not set up for correlated RP-MDCEV models yet.", "\n")
	if(object$n_classes > 1){
		est_sim <- est_sim[grepl(class, est_sim$parms),]
	}

	sim_data <- ProcessSimulationData(est_sim, object, policies, nsims)
	df_common <- sim_data
	df_common$df_indiv <- NULL

	df_indiv <- sim_data$df_indiv

	sim_options <- list(n_classes = object$n_classes,
					model_num = object$stan_data$model_num,
					price_change_only = policies$price_change_only)

	output <- list(df_indiv = df_indiv,
				   df_common = df_common,
				   sim_options = sim_options)

return(output)
}

#' @title ProcessSimulationData
#' @description Internal function for WTP simulation
#' @inheritParams PrepareSimulationData
#' @param est_sim Cleaned up parameter simulations from PrepareSimulationData
#' @keywords internal
ProcessSimulationData <- function(est_sim, object, policies, nsims){

	J <- object$stan_data$J
	I <- object$stan_data$I
	NPsi_ij = object$stan_data$NPsi_ij
	model_num <- object$stan_data$model_num
	random_parameters <- object$random_parameters
	gamma_nonrandom <- object$stan_data$gamma_nonrandom
	alpha_nonrandom <- object$stan_data$alpha_nonrandom
	npols <- length(policies$price_p)
	psi_ascs = object$stan_data$psi_ascs

	# scales (always fixed parameter)
	if (object$stan_data$fixed_scale1 == 0){
#		scale_sim = est_sim$scale
		scale_sims <- t(GrabParms(est_sim, "scale"))
	} else if (object$stan_data$fixed_scale1 == 1){
		scale_sims = matrix(1, nsims, 1)
	}

	# Get fixed parameters for gammas
	if (model_num == 2){
		gamma_sim_nonrandom <- matrix(1, nsims, J)
	} else if (model_num != 2 & gamma_nonrandom == 1){
#		gamma_sim_nonrandom = est_sim$gamma
		gamma_sim_nonrandom <- GrabParms(est_sim, "gamma")
		if (object$stan_data$gamma_ascs == 0){
			gamma_sim_nonrandom <- matrix(rep(gamma_sim_nonrandom, each=J), ncol=J, byrow=TRUE)
		}
	} else if (gamma_nonrandom == 0)
		gamma_sim_nonrandom <- NULL

	# alphas
	if (model_num == 4){
		alpha_sim_nonrandom <- matrix(1e-3, nsims, J+1)
	} else if (model_num != 4 & alpha_nonrandom == 1){
#		alpha_sim_nonrandom = est_sim$alpha
		alpha_sim_nonrandom <- GrabParms(est_sim, "alpha")

		if (model_num == 1 || model_num == 5) {
			alpha_sim_nonrandom <- cbind(alpha_sim_nonrandom, matrix(0, nsims, J) )
		} else if (model_num == 3){
			alpha_sim_nonrandom <- matrix(rep(alpha_sim_nonrandom, each=J+1), ncol=J+1, byrow=TRUE)
		}
	} else if (alpha_nonrandom == 0){
		alpha_sim_nonrandom <- NULL
	}

	if (!is.null(alpha_sim_nonrandom)){
		alpha_sim_rand <- NULL
	}

	if (!is.null(gamma_sim_nonrandom)){
		gamma_sim_rand <- NULL
	}
# Get individual parameters
# ensure psi variables in same form as psi_sim_rand
if (random_parameters == "fixed"){
	psi_sim_temp <- GrabParms(est_sim, "psi") # change back to est_sim
	psi_sim_temp = replicate(I, psi_sim_temp, simplify=FALSE)

	if (object[["parms_info"]][["n_vars"]][["n_phi"]] > 0){
		phi_sim_temp <- GrabParms(est_sim, "phi") # change back to est_sim
		phi_sim_temp = replicate(I, phi_sim_temp, simplify=FALSE)
	}
} else if (random_parameters != "fixed"){

	est_sim_mu_tau <- est_sim %>%
		dplyr::filter(grepl(c("mu|tau"), parms)) %>%
		tidyr::separate(parms, into = c("parms", "parm_id"), sep = "\\.") %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id)) %>%
		tidyr::spread(parms, value)  %>%
		dplyr::arrange(sim_id)

	if(random_parameters == "corr"){

		num_rand <- object[["stan_fit"]]@par_dims[["mu"]]

		est_sim_tau <- est_sim_mu_tau %>%
			dplyr::select(sim_id, .data$parm_id, .data$tau) %>%
			dplyr::group_split(sim_id)

		est_sim_lomega <- est_sim %>%
			dplyr::filter(grepl(c("L_Omega"), parms))   %>%
			dplyr::arrange(sim_id) %>%
			dplyr::group_split(sim_id)

		sim_id <- est_sim %>%
			dplyr::arrange(sim_id) %>%
			dplyr::distinct(sim_id)

		L <- mapply(function(x, y){
						l_omega <- matrix(y$value, nrow = num_rand, byrow=F)
						L <- as.vector(x$tau %*% l_omega)
						return(L)
					}, est_sim_tau, est_sim_lomega)

		L <- matrix(unlist(L), nrow = nrow(sim_id), byrow = T )
		colnames(L) <- c(paste0(rep("parm_id", num_rand), 1:num_rand))

		est_sim_tau <- bind_cols(sim_id, as_tibble(L)) %>%
			tidyr::gather(parm_id, tau, -sim_id) %>%
			dplyr::mutate(parm_id = as.numeric(gsub("[^0-9]", "", .data$parm_id))) %>%
			dplyr::arrange(sim_id)

		est_sim_mu_tau <- est_sim_mu_tau %>%
			dplyr::select(sim_id, .data$parm_id, .data$mu) %>%
			dplyr::left_join(est_sim_tau, by = c("sim_id", "parm_id"))
	}

	est_sim <- est_sim %>%
		dplyr::filter(grepl(c("^z|sim_id"), parms)) %>%
		tidyr::separate(parms, into = c("parms", "id", "parm_id"), sep = "\\.") %>%
		tidyr::spread(parms, value) %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id),
			   id= as.numeric(id)) %>%
		dplyr::arrange(id, sim_id, .data$parm_id)  %>%
		dplyr::left_join(est_sim_mu_tau, by = c("sim_id", "parm_id")) %>%
		dplyr::arrange(id, sim_id, .data$parm_id)	%>%
		dplyr::mutate(beta = .data$mu + .data$z *.data$tau,
			   parms = rep(object[["parms_info"]][["parm_names"]][["sd_names"]],
			   			nsims*object[["n_individuals"]] )) %>%
		dplyr::select(-.data$tau, -.data$mu, -.data$z)

	# Transform gamma and alpha estimates
	if(gamma_nonrandom == 0){
		est_sim <- est_sim %>%
			dplyr::mutate(beta = ifelse(grepl(c("gamma"), parms), exp(beta), beta))
	}
	if(alpha_nonrandom == 0){
		est_sim <- est_sim %>%
			dplyr::mutate(beta = ifelse(grepl(c("alpha"), parms), 1 / (1 + exp(-beta)), beta))
	}

	if (gamma_nonrandom == 0){
		gamma_sim_rand <- GrabIndividualParms(est_sim, "gamma")
		if (object$stan_data$gamma_ascs == 0){
			gamma_sim_rand <- lapply(gamma_sim_rand, function(x){
				matrix(rep(x, each=J), ncol=J, byrow=TRUE)})
		}
		gamma_sim_rand <- list(lapply(gamma_sim_rand, as.matrix))
		names(gamma_sim_rand) <- "gamma_sims"
	}

	if (alpha_nonrandom == 0){
		alpha_sim_rand <- GrabIndividualParms(est_sim, "alpha")
		if (model_num == 1 || model_num == 5) {
			alpha_sim_rand <- lapply(alpha_sim_rand, function(x){
				cbind(x, matrix(0, nsims, J) )})
		} else if (model_num == 3){
			alpha_sim_rand <- lapply(alpha_sim_rand, function(x){
				matrix(rep(x,each=J+1), ncol=J+1, byrow=TRUE)})
		}
		alpha_sim_rand <- list(lapply(alpha_sim_rand, as.matrix))
		names(alpha_sim_rand) <- "alpha_sims"
	}

	psi_sim_temp <- GrabIndividualParms(est_sim, "psi")
	psi_sim_temp <- lapply(psi_sim_temp, as.matrix)

	# Get parameters for phis
	if (model_num == 5){
		phi_sim_temp <- list(GrabIndividualParms(est_sim, "phi"))
	} else if (model_num != 5){
		phi_sim_temp <- list(array(0, dim = c(0,0)))
	}
}

# Create Psi variables
if (object[["parms_info"]][["n_vars"]][["n_psi"]] > 0){
	dat_vars <- tibble(id = rep(1:I, each = J))

	if (NPsi_ij > 0)
		dat_vars <- bind_cols(dat_vars, as_tibble(object$stan_data$dat_psi))

	dat_vars <- dat_vars %>%
		group_split(id, .keep = F)

		psi_sims <- mapply(CreatePsi, dat_vars, psi_sim_temp,
					  J = J, NPsi_ij=NPsi_ij, psi_ascs=psi_ascs, npols = npols,
					  SIMPLIFY = FALSE)
} else {
	psi_sims = replicate(I, matrix(0, nsims, J), simplify=FALSE)
}

psi_sims <- list(psi_sims)
names(psi_sims) <- "psi_sims"


# phi
# Get parameters for phis
if (object[["parms_info"]][["n_vars"]][["n_phi"]] > 0){

	dat_id <- tibble(id = rep(1:I, each = J))
	dat_vars <- bind_cols(dat_id, as_tibble(object$stan_data$dat_phi))

	dat_vars <- dat_vars %>%
		group_split(id, .keep = F)

	phi_sims <- mapply(function(x, y){
		phi_temp <- exp(as.matrix(x) %*% t(as.matrix(y))) # exponentiate to get back to proper transformation
		phi_temp <- t(phi_temp)
	}, dat_vars, phi_sim_temp, SIMPLIFY = FALSE)

	phi_sims <- list(phi_sims)
	names(phi_sims) <- "phi_sims"

} else if (object[["parms_info"]][["n_vars"]][["n_phi"]] == 0){
	phi_sims <- array(1, dim = c(nsims,J))
	phi_sims = replicate(I, phi_sims, simplify=FALSE)
	phi_sims <- list(phi_sims)
	names(phi_sims) <- "phi_sims"
}

# Set up psi_p and phi_p

phi_p_sims <- replicate(I, matrix(0, 0, 0), simplify=FALSE)
psi_p_sims <- replicate(I, matrix(0, 0, 0), simplify=FALSE)

if (policies$price_change_only == FALSE) {

	if (model_num < 5){
		dat_vars <- tibble(id = rep(1:I, each = J))

		if (NPsi_ij > 0){
			dat_vars_temp <- mapply(cbind, policies[["dat_psi_p"]], "policy"=1:npols, SIMPLIFY=F)
			dat_vars_temp <- lapply(dat_vars_temp, function(x){
				bind_cols(dat_vars, as_tibble(x))})

		dat_vars <- do.call(rbind, dat_vars_temp)
		}

		dat_vars <- dat_vars %>%
			group_split(id, .keep = F)

		psi_p_sims <- mapply(CreatePsi, dat_vars, psi_sim_temp,
							J = J, NPsi_ij=NPsi_ij, psi_ascs=psi_ascs,npols = npols,
							SIMPLIFY = FALSE)



	} else if (model_num == 5){
		dat_vars <- tibble(id = rep(1:I, each = J))

			dat_vars_temp <- mapply(cbind, policies[["dat_phi_p"]], "policy"=1:npols, SIMPLIFY=F)
			dat_vars_temp <- lapply(dat_vars_temp, function(x){
				bind_cols(dat_vars, as_tibble(x))})

			dat_vars <- do.call(rbind, dat_vars_temp)


		dat_vars <- dat_vars %>%
			group_split(id, .keep = F)

		phi_p_sims <- mapply(function(x, y){
			x_i = x %>%
				group_split(policy, .keep = F)
			phi_p_temp <- lapply(x_i, function(xx){
			phi_p_temp <- exp(as.matrix(xx) %*% t(as.matrix(y))) # exponentiate to get back to proper transformation
			phi_p_temp <- t(phi_p_temp)
			return(phi_p_temp)
			})
		}, dat_vars, phi_sim_temp, SIMPLIFY = FALSE)


	}
}
psi_p_sims <- list(psi_p_sims)
names(psi_p_sims) <- "psi_p_sims"
phi_p_sims <- list(phi_p_sims)
names(phi_p_sims) <- "phi_p_sims"

# Set baseline individual data into lists
income <- list(as.list(object$stan_data$income))
names(income) <- "income"

quant_j <- list(CreateListsRow(object$stan_data$quant_j))
names(quant_j) <- "quant_j"

price <- cbind(1, object$stan_data$price_j) #add numeraire price to price matrix (<-1)
price <- list(CreateListsRow(price))
names(price) <- "price"

# Pull individual level data into one list
df_indiv <- c(income, quant_j, price, psi_sims, phi_sims, psi_p_sims, phi_p_sims,
			  gamma_sim_rand, alpha_sim_rand)

out <- list(df_indiv = df_indiv,
			price_p_list = policies$price_p,
			gamma_sim_nonrandom = gamma_sim_nonrandom,
			alpha_sim_nonrandom = alpha_sim_nonrandom,
			scale_sims = scale_sims)
return(out)
}
