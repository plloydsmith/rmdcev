#' @title PrepareSimulationData
#' @description Prepare Data for WTP simulation
#' @param stan_est Stan fit model from FitMDCEV
#' @param policies list containing
#' price_p with additive price increases, and
#' dat_psi_p with new psi data
#' @param nsims Number of simulation draws to use for parameter uncertainty
#' @return A list with individual-specific data (df_indiv) and common data (df_common)
#' and n_classes for number of classes and model_num for model type
#' @export
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#' mdcev_est <- FitMDCEV(psi_formula = ~ 1,
#' data = subset(data_rec, id < 500),
#' model = "hybrid0",
#' algorithm = "MLE")
#'
#' policies <- CreateBlankPolicies(npols = 2,
#' ngoods = mdcev_est[["stan_data"]][["J"]],
#' dat_psi = mdcev_est[["stan_data"]][["dat_psi"]],
#' price_change_only = TRUE)
#'
#' df_sim <- PrepareSimulationData(mdcev_est, policies)
#'}
PrepareSimulationData <- function(stan_est,
								  policies,
								  nsims = 30){
#stan_est <- result
#	if (stan_est$random_parameters == "corr")
#		stop("Demand and welfare simulation not set up for correlated RP-MDCEV models yet.", "\n")

	model_num <- stan_est$stan_data$model_num

	# Checks on simulation options
	if (stan_est$algorithm == "MLE" && nsims > stan_est$n_draws){
		nsims <- stan_est$n_draws
		warning("Number of simulations > Number of mvn draws from stan_est. nsims has been set to: ", nsims)
	} else if (stan_est$algorithm == "Bayes" && nsims > max(stan_est$est_pars$sim_id)) {
		nsims <- max(stan_est$est_pars$sim_id)
		warning("Number of simulations > Number of posterior draws from stan_est. nsims has been set to: ", nsims)
	}

	# Sample from parameter estimate draws
	est_sim <- stan_est$est_pars %>%
		dplyr::distinct(sim_id) %>%
		dplyr::sample_n(., nsims) %>%
		dplyr::left_join(stan_est$est_pars, by = "sim_id")

	if(stan_est$n_classes == 1){
		sim_welfare <- ProcessSimulationData(est_sim, stan_est, policies, nsims)
		df_common <- sim_welfare
		df_common$df_indiv <- NULL

		df_indiv <- sim_welfare$df_indiv

	} else if(stan_est$n_classes > 1){

		est_sim_lc <- suppressWarnings(est_sim %>% # suppress warnings about scale not having a class parameter
									   	dplyr::filter(!stringr::str_detect(.data$parms, "delta")) %>%
									   	tidyr::separate(.data$parms, into = c("class","parms"), sep = 7) %>%
									   	dplyr::mutate(class = gsub("[^[:alnum:]]", "", class)))

		est_sim_lc <- split( est_sim_lc , f = est_sim_lc$class )
		names(est_sim_lc) <- rep("est_sim", stan_est$n_classes)

		est_sim_lc <- purrr::map(est_sim_lc, function(x){ x %>%
				select(-class)})

		sim_welfare <- purrr::map(est_sim_lc, ProcessSimulationData, stan_est, policies, nsims)

		df_common <- purrr::map(sim_welfare, `[`, c("price_p_list", "gamma_sim_fixed", "alpha_sim_fixed", "scale_sim"))
		names(df_common) <- rep("df_common", stan_est$n_classes)

		df_indiv <- purrr::flatten(purrr::map(sim_welfare, `[`, c("df_indiv")))
	}

	sim_options <- list(n_classes = stan_est$n_classes,
					model_num = model_num,
					price_change_only = policies$price_change_only)

	df_wtp <- list(df_indiv = df_indiv,
			   df_common = df_common,
			   sim_options = sim_options)

return(df_wtp)
}

#' @title ProcessSimulationData
#' @description Internal function for WTP simulation
#' @inheritParams PrepareSimulationData
#' @param est_sim Cleaned up parameter simulations from PrepareSimulationData
#' @keywords internal
ProcessSimulationData <- function(est_sim, stan_est, policies, nsims){

#	stan_est <- mdcev_corr2

	J <- stan_est$stan_data$J
	I <- stan_est$stan_data$I
	model_num <- stan_est$stan_data$model_num
	random_parameters <- stan_est$random_parameters
	gamma_fixed <- stan_est$stan_data$gamma_fixed
	alpha_fixed <- stan_est$stan_data$alpha_fixed
	npols <- length(policies$price_p)


	# scales (always fixed parameter)
	if (stan_est$stan_data$fixed_scale1 == 0){
		scale_sim <- t(GrabParms(est_sim, "scale"))
	} else if (stan_est$stan_data$fixed_scale1 == 1){
		scale_sim = matrix(1, nsims, 1)
	}


	# Get fixed parameters for gammas
	if (model_num == 2){
		gamma_sim_fixed <- matrix(1, nsims, J)
	} else if (model_num != 2 & gamma_fixed == 1){
		gamma_sim_fixed <- t(GrabParms(est_sim, "gamma"))
	} else if (gamma_fixed == 0){
		gamma_sim_fixed <- NULL
	}


	# alphas
	if (model_num == 4){
		alpha_sim_fixed <- matrix(1e-3, nsims, J+1)
	} else if (model_num != 4 & alpha_fixed == 1){
		alpha_sim_fixed <- t(GrabParms(est_sim, "alpha"))

		if (model_num == 1) {
			alpha_sim_fixed <- cbind(alpha_sim_fixed, matrix(0, nsims, J) )
		} else if (model_num == 3){
			alpha_sim_fixed <- matrix(rep(alpha_sim_fixed, each=J+1), ncol=J+1, byrow=TRUE)
		}
	} else if (alpha_fixed == 0){
		alpha_sim_fixed <- NULL
	}

	# Get individual parameters
if(random_parameters != "fixed"){

	est_sim_mu_tau <- est_sim %>%
		dplyr::filter(grepl(c("mu|tau"), parms)) %>%
		tidyr::separate(parms, into = c("parms", "parm_id"), sep = "\\.") %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id)) %>%
		tidyr::spread(parms, value)  %>%
		dplyr::arrange(sim_id)

	if(random_parameters == "corr"){

		num_rand <- stan_est[["stan_fit"]]@par_dims[["mu"]]

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

		L <- purrr::map2(est_sim_tau, est_sim_lomega, function(x, y){
			l_omega <- matrix(y$value, nrow = num_rand, byrow=F)
			L <- as.vector(x$tau %*% l_omega)
			return(L)
		})

	L <- matrix(unlist(L), nrow = nrow(sim_id), byrow = T )
	colnames(L) <- c(paste0(rep("parm_id", num_rand), 1:num_rand))

	est_sim_tau <- bind_cols(sim_id, tbl_df(L)) %>%
		tidyr::gather(.data$parm_id, .data$tau, -sim_id) %>%
		dplyr::mutate(parm_id = as.numeric(gsub("[^0-9]", "", .data$parm_id))) %>%
		dplyr::arrange(sim_id)

	est_sim_mu_tau <- est_sim_mu_tau %>%
		dplyr::select(sim_id, .data$parm_id, .data$mu) %>%
		dplyr::left_join(est_sim_tau, by = c("sim_id", "parm_id"))

	} else if(random_parameters == "uncorr"){

	est_sim_mu_tau <- est_sim %>%
		dplyr::filter(grepl(c("mu|tau"), parms)) %>%
		tidyr::separate(parms, into = c("parms", "parm_id"), sep = "\\.") %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id)) %>%
		tidyr::spread(parms, value)  %>%
		dplyr::arrange(sim_id)
	}

	est_sim <- est_sim %>%
		dplyr::filter(stringr::str_detect(parms, "^z")) %>%
		tidyr::separate(parms, into = c("parms", "id", "parm_id"), sep = "\\.") %>%
		tidyr::spread(parms, value) %>%
		dplyr::mutate(parm_id = as.numeric(.data$parm_id),
			   id= as.numeric(id)) %>%
		dplyr::arrange(id, sim_id, .data$parm_id)  %>%
		dplyr::left_join(est_sim_mu_tau, by = c("sim_id", "parm_id")) %>%
		dplyr::arrange(id, sim_id, .data$parm_id)	%>%
		dplyr::mutate(beta = .data$mu + .data$z *.data$tau,
			   parms = rep(stan_est[["parms_info"]][["parm_names"]][["sd_names"]],
			   			nsims*stan_est[["n_individuals"]] )) %>%
		dplyr::select(-.data$tau, -.data$mu, -.data$z)

	# Transform gamma and alpha estimates
	if(stan_est[["stan_data"]][["gamma_fixed"]]==0){
		est_sim <- est_sim %>%
			dplyr::mutate(beta = ifelse(grepl(c("gamma"), parms), exp(beta), beta))
	}
	if(stan_est[["stan_data"]][["alpha_fixed"]]==0){
		est_sim <- est_sim %>%
			dplyr::mutate(beta = ifelse(grepl(c("alpha"), parms), 1 / (1 + exp(-beta)), beta))
	}

	if (gamma_fixed == 0){
		gamma_sim_rand <- est_sim %>%
			dplyr::filter(grepl(c("gamma"), parms)) %>%
			dplyr::select(id, sim_id, .data$parm_id, beta) %>%
			tidyr::spread(.data$parm_id, beta) %>%
			dplyr::select(-sim_id) %>%
			dplyr::group_split(id, keep = F)

		gamma_sim_rand <- purrr::map(gamma_sim_rand, ~as.matrix(.))

		gamma_sim_rand <- list(gamma_sim_rand)
		names(gamma_sim_rand) <- "gamma_sim"
	}

	if (alpha_fixed == 0){
		alpha_sim_rand <- est_sim %>%
			dplyr::filter(grepl(c("alpha"), parms)) %>%
			dplyr::select(id, sim_id, .data$parm_id, beta) %>%
			tidyr::spread(.data$parm_id, beta) %>%
			dplyr::select(-sim_id) %>%
			dplyr::group_split(id, keep = F)

		alpha_sim_rand <- purrr::map(alpha_sim_rand, ~as.matrix(.))

		alpha_sim_rand <- list(alpha_sim_rand)
		names(alpha_sim_rand) <- "alpha_sim"
	}

	psi_sim_rand <- est_sim %>%
		dplyr::filter(grepl(c("psi"), parms)) %>%
		dplyr::select(id, sim_id, .data$parm_id, beta) %>%
		tidyr::spread(.data$parm_id, beta) %>%
		dplyr::select(-sim_id) %>%
		dplyr::group_split(id, keep = F)
}

if (random_parameters == "fixed"){
	# psi
	psi_temp <- GrabParms(est_sim, "psi") # change back to est_sim


	psi_temp <- CreateListsCol(psi_temp)
	psi_sim <- purrr::map(psi_temp, MultiplyMatrix, mat_temp = stan_est$stan_data$dat_psi, n_rows = I)

	psi_sim <- DoCbind(psi_sim)
	psi_sim <- CreateListsRow(psi_sim)
	psi_sim <- purrr::map(psi_sim, function(x){matrix(x , nrow = nsims, byrow = TRUE)})

	psi_sim <- list(psi_sim)
	names(psi_sim) <- "psi_sim"

	if (policies$price_change_only == FALSE) {
	# psi_p
	psi_p_sim <- purrr::map(psi_temp, function(psi){ purrr::map(policies[["dat_psi_p"]], MultiplyMatrix, x = psi, n_rows = I)})
	psi_p_sim <- purrr::map(psi_p_sim, DoCbind)
	psi_p_sim <- DoCbind(psi_p_sim)
	psi_p_sim <- CreateListsRow(psi_p_sim)
	psi_p_sim <- purrr::map(psi_p_sim, function(x){aperm(array(x, dim = c(J, npols, nsims)), perm=c(2,1,3))})

	# Ensure psi_p_sim is a list of J lists each with nsims lists of npol X ngood matrices
	psi_p_sim <- purrr::map(psi_p_sim, function(x){lapply(seq_len(nsims), function(i) x[,,i])})

	psi_p_sim <- list(psi_p_sim)
	names(psi_p_sim) <- "psi_p_sim"

	} else if (policies$price_change_only == TRUE){
		psi_p_sim <- NULL
	}

} else if (random_parameters != "fixed"){

	dat_id <- tibble(id = rep(1:stan_est$n_individuals, each = stan_est$stan_data$J))

	dat_psi <- bind_cols(dat_id, tbl_df(stan_est$stan_data$dat_psi)) %>%
		group_split(id, keep = F)

	psi_sim <- purrr::map2(psi_sim_rand, dat_psi, function(x, y){

		psi_sim <- CreateListsRow(x)
		dat_psi_1 <- as.matrix(y)

		out <- purrr::map(psi_sim, function(xx){
			psi <- dat_psi_1 %*% t(as.matrix(xx))} )

		out <- matrix(unlist(out), byrow=TRUE, nrow=length(out) )
		return(out)
	})

psi_sim <- list(psi_sim)
names(psi_sim) <- "psi_sim"

# Can't change psi_p for random parameters yet
psi_p_sim <- NULL
}

if (!is.null(alpha_sim_fixed)){
	alpha_sim_rand <- NULL
}

if (!is.null(gamma_sim_fixed)){
	gamma_sim_rand <- NULL
}


# Set baseline individual data into lists
income <- list(as.list(stan_est$stan_data$income))
names(income) <- "income"

quant_j <- list(CreateListsRow(stan_est$stan_data$quant_j))
names(quant_j) <- "quant_j"

price <- cbind(1, stan_est$stan_data$price_j) #add numeraire price to price matrix (<-1)
price <- list(CreateListsRow(price))
names(price) <- "price"

# Pull individual level data into one list
df_indiv <- c(income, quant_j, price, psi_sim, psi_p_sim,
			  gamma_sim_rand, alpha_sim_rand)

out <- list(df_indiv = df_indiv,
			price_p_list = policies$price_p,
			gamma_sim_fixed = gamma_sim_fixed,
			alpha_sim_fixed = alpha_sim_fixed,
			scale_sim = scale_sim)
return(out)
}
