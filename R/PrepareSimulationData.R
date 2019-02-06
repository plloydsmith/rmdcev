PrepareSimulationData <- function(est_sim, stan_est, policies, nsims){

	gamma_sim <- est_sim %>%
		filter(str_detect(parms, "gamma")) %>%
		spread(sim_id, value) %>%
		select(-parms) %>%
		as.matrix(.) %>%
		t(.)

	if (stan_est$stan_data$model_num != 4){
		alpha_sim <- est_sim %>%
			filter(str_detect(parms, "alpha")) %>%
			spread(sim_id, value) %>%
			select(-parms) %>%
			as.matrix(.) %>%
			t(.)
		if (stan_est$stan_data$model_num == 1){
			alpha_sim <- cbind(alpha_sim, matrix(0, nsims, stan_est$stan_data$J) )
		} else if (stan_est$stan_data$model_num == 3)
			alpha_sim <- matrix(rep(alpha_sim,each=stan_est$stan_data$J+1), ncol=stan_est$stan_data$J+1, byrow=TRUE)

	} else if (stan_est$stan_data$model_num ==4)
		alpha_sim = matrix(1e-6, nsims, stan_est$stan_data$J+1)

	if (stan_est$stan_data$fixed_scale == 0){
		scale_sim <- est_sim %>%
			filter(str_detect(parms, "scale")) %>%
			spread(sim_id, value) %>%
			select(-parms) %>%
			as.matrix(.) %>%
			t(.)
	} else if (stan_est$stan_data$fixed_scale == 1){
		scale_sim = matrix(1, nsims, 1)
	}

	psi_temp <- est_sim %>%
		filter(str_detect(parms, "psi")) %>%
		spread(sim_id, value) %>%
		select(-parms) %>%
		as.matrix(.)

	npols <- length(policies$price_p)

	psi_temp <- CreateListsCol(psi_temp)
	psi_sim <- map(psi_temp, MultiplyMatrix, mat_temp = stan_est$stan_data$dat_psi, nrows = stan_est$stan_data$I)

	psi_sim <- DoCbind(psi_sim)
	psi_sim <- CreateListsRow(psi_sim)
	psi_sim <- map(psi_sim, function(x){matrix(x , nrow = nsims, byrow = TRUE)})

	psi_p_sim <- map(psi_temp, function(psi){ map(policies[["dat_psi_p"]], MultiplyMatrix, x = psi, nrows = stan_est$stan_data$I)})
	psi_p_sim <- map(psi_p_sim, DoCbind)
	psi_p_sim <- DoCbind(psi_p_sim)
	psi_p_sim <- CreateListsRow(psi_p_sim)
	psi_p_sim <- map(psi_p_sim, function(x){aperm(array(x, dim = c(stan_est$stan_data$J, npols, nsims)), perm=c(2,1,3))})

	# Ensure psi_p_sim is a list of J lists each with nsims lists of npol X ngood matrices
	psi_p_sim <- map(psi_p_sim, function(x){lapply(seq_len(nsims), function(i) x[,,i])})


	# Put in a list for each simulation
	gamma_sim_list <- CreateListsRow(gamma_sim)
	alpha_sim_list <- CreateListsRow(alpha_sim)
	scale_sim <- scale_sim

	psi_sim <- list(psi_sim)
	names(psi_sim) <- "psi_sim"

	psi_p_sim <- list(psi_p_sim)
	names(psi_p_sim) <- "psi_p_sim"

	###########################################################################
	# Set baseline data into lists
	inc <- list(as.list(stan_est$stan_data$inc))
	names(inc) <- "inc"

	quant_j <- list(CreateListsRow(stan_est$stan_data$j_quant))
	names(quant_j) <- "quant_j"

	price <- cbind(1, stan_est$stan_data$j_price) #add numeraire price to price matrix (<-1)
	price <- list(CreateListsRow(price))
	names(price) <- "price"

	###########################################################################
	# Pull individual level data into one list
	df_indiv <- c(inc, quant_j, price, psi_sim, psi_p_sim)

	out <- list(df_indiv = df_indiv,
				price_p_list = policies$price_p,
				gamma_sim_list = gamma_sim_list,
				alpha_sim_list = alpha_sim_list,
				scale_sim = scale_sim)
	return(out)
}

