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

PrepareSimulationData <- function(stan_est, policies, nsims = 30){

	# Checks on simulation options
	model_num <- stan_est$stan_data$model_num

	if (nsims > stan_est$n_draws) {
		nsims <- stan_est$n_draws
		warning("Number of simulations > Number of Draws from stan_est. nsims has been set to: ", nsims)
	}

	# Sample from parameter estimate draws
	est_sim <- stan_est$est_pars %>%
		distinct(.data$sim_id) %>%
		sample_n(., nsims ) %>%
		left_join(stan_est$est_pars, by = "sim_id")

	if(stan_est$n_classes == 1){
		sim_welfare <- ProcessSimulationData(est_sim, stan_est, policies, nsims)
		df_common <- sim_welfare
		df_common$df_indiv <- NULL

		df_indiv <- sim_welfare$df_indiv

	} else if(stan_est$n_classes > 1){

		est_sim_lc <- suppressWarnings(est_sim %>% # suppress warnings about scale not having a class parameter
									   	filter(!str_detect(.data$parms, "beta")) %>%
									   	separate_(.data$parms, into = c("parms", "class", "good")) %>%
									   	mutate(good = ifelse(is.na(as.numeric(.data$good)), "0", .data$good )) %>%
									   	unite_(parms, parms, good))

		est_sim_lc <- split( est_sim_lc , f = est_sim_lc$class )
		names(est_sim_lc) <- rep("est_sim", stan_est$n_classes)

		est_sim_lc <- map(est_sim_lc, function(x){ x %>%
				select(-class)})

		sim_welfare <- map(est_sim_lc, ProcessSimulationData, stan_est, policies, nsims)

		df_common <- map(sim_welfare, `[`, c("price_p_list", "gamma_sim_list", "alpha_sim_list", "scale_sim"))
		names(df_common) <- rep("df_common", stan_est$n_classes)

		df_indiv <- flatten(map(sim_welfare, `[`, c("df_indiv")))


	}
	df_wtp <- list(df_indiv = df_indiv,
				   df_common = df_common,
				   n_classes = stan_est$n_classes,
				   model_num = model_num)
return(df_wtp)
}



#'@importFrom rlang .data
#'
#'
#'
ProcessSimulationData <- function(est_sim, stan_est, policies, nsims){

	J <- stan_est$stan_data$J
	I <- stan_est$stan_data$I

	# gammas
	if (stan_est$stan_data$model_num == 2)
		gamma_sim <- matrix(1, nsims, J)
	else if (stan_est$stan_data$model_num != 2)
		gamma_sim <- t(GrabParms(est_sim, "gamma"))

	gamma_sim_list <- CreateListsRow(gamma_sim)	# Put in a list for each simulation

	# alphas
	if (stan_est$stan_data$model_num != 4){
		alpha_sim <- t(GrabParms(est_sim, "alpha"))

		if (stan_est$stan_data$model_num == 1)
			alpha_sim <- cbind(alpha_sim, matrix(0, nsims, J) )
		else if (stan_est$stan_data$model_num == 3)
			alpha_sim <- matrix(rep(alpha_sim,each=J+1), ncol=J+1, byrow=TRUE)

	} else if (stan_est$stan_data$model_num ==4)
		alpha_sim <- matrix(1e-6, nsims, J+1)

	alpha_sim_list <- CreateListsRow(alpha_sim)

	# scales
	if (stan_est$stan_data$fixed_scale == 0)
		scale_sim <- t(GrabParms(est_sim, "scale"))
	else if (stan_est$stan_data$fixed_scale == 1)
		scale_sim = matrix(1, nsims, 1)

	# psi
	psi_temp <- GrabParms(est_sim, "psi")

	npols <- length(policies$price_p)

	psi_temp <- CreateListsCol(psi_temp)
	psi_sim <- map(psi_temp, MultiplyMatrix, mat_temp = stan_est$stan_data$dat_psi, n_rows = I)

	psi_sim <- DoCbind(psi_sim)
	psi_sim <- CreateListsRow(psi_sim)
	psi_sim <- map(psi_sim, function(x){matrix(x , nrow = nsims, byrow = TRUE)})

	psi_sim <- list(psi_sim)
	names(psi_sim) <- "psi_sim"

	# psi_p
	psi_p_sim <- map(psi_temp, function(psi){ map(policies[["dat_psi_p"]], MultiplyMatrix, x = psi, n_rows = I)})
	psi_p_sim <- map(psi_p_sim, DoCbind)
	psi_p_sim <- DoCbind(psi_p_sim)
	psi_p_sim <- CreateListsRow(psi_p_sim)
	psi_p_sim <- map(psi_p_sim, function(x){aperm(array(x, dim = c(J, npols, nsims)), perm=c(2,1,3))})

	# Ensure psi_p_sim is a list of J lists each with nsims lists of npol X ngood matrices
	psi_p_sim <- map(psi_p_sim, function(x){lapply(seq_len(nsims), function(i) x[,,i])})

	psi_p_sim <- list(psi_p_sim)
	names(psi_p_sim) <- "psi_p_sim"

	# Set baseline individual data into lists
	inc <- list(as.list(stan_est$stan_data$inc))
	names(inc) <- "inc"

	quant_j <- list(CreateListsRow(stan_est$stan_data$j_quant))
	names(quant_j) <- "quant_j"

	price <- cbind(1, stan_est$stan_data$j_price) #add numeraire price to price matrix (<-1)
	price <- list(CreateListsRow(price))
	names(price) <- "price"

	# Pull individual level data into one list
	df_indiv <- c(inc, quant_j, price, psi_sim, psi_p_sim)

	out <- list(df_indiv = df_indiv,
				price_p_list = policies$price_p,
				gamma_sim_list = gamma_sim_list,
				alpha_sim_list = alpha_sim_list,
				scale_sim = scale_sim)
	return(out)
}
