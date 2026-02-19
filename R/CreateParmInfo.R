#' @title CreateParmInfo
#' @description Create parameter names and number of parameters
#' @param stan_data data for model
#' @param alt_names name of alternatives.
#' @inheritParams mdcev
#' @return A list of parameter names and numbers
#' @noRd
CreateParmInfo <- function(stan_data, alt_names, algorithm, random_parameters){

	J <- stan_data$J

	if (stan_data$psi_ascs == 1){
		n_psi_asc <- J - 1
		psi_asc_names <- alt_names[-1]
	} else if (stan_data$psi_ascs == 0){
		n_psi_asc <- 0
		psi_asc_names <- NULL
	}

	if (stan_data$NPsi_ij > 0)
		psi_non_asc_names <- colnames(stan_data[["dat_psi"]])
	else
		psi_non_asc_names <- NULL

	# add in alternative-specific attributes
	n_psi <- n_psi_asc + stan_data$NPsi_ij
	psi_names <- c(psi_asc_names, psi_non_asc_names)

	if (n_psi > 0)
		psi_names <- paste0(rep('psi', n_psi), sep="_", psi_names)

	# alpha
	if (stan_data$model_num == 1 || stan_data$model_num == 5){
		n_alpha <- 1
		alpha_names <- 'alpha_num'
	} else if (stan_data$model_num == 2){
		n_alpha <- J + 1
		alpha_names <- c('alpha_num', paste0(rep('alpha', n_alpha-1), sep = "_", alt_names))
	} else if (stan_data$model_num == 3){
		n_alpha <- 1
		alpha_names <- 'alpha'
	} else if (stan_data$model_num == 4){
		n_alpha <- 0
		alpha_names <- NULL
	}

	# gamma
	if (stan_data$model_num == 2){
		n_gamma <- 0
		gamma_names <- NULL
	} else {
		if (stan_data$gamma_ascs == 0){
			n_gamma <- 1
			gamma_names <- 'gamma'
		} else if (stan_data$gamma_ascs == 1){
			n_gamma <- J
			gamma_names <- paste0(rep('gamma', n_gamma), sep = "_", alt_names)
		}
	}

	# Phi parameters
	if (stan_data$model_num != 5 || stan_data$NPhi == 0){
		phi_names <- NULL
		n_phi <- 0
	} else {
		n_phi <- stan_data$NPhi
		phi_names <- colnames(stan_data[["dat_phi"]])
		phi_names <- paste0(rep('phi', n_phi), sep="_", phi_names)
	}

	if (stan_data$fixed_scale1 == 1){
		n_scale <- 0
		scale_names <- NULL
	} else {
		n_scale <- 1
		scale_names <- "scale"
	}

	parm_names <- list(psi_names   = psi_names,
					   phi_names   = phi_names,
					   gamma_names = gamma_names,
					   alpha_names = alpha_names,
					   scale_names = scale_names,
					   all_names   = c(psi_names, phi_names, gamma_names, alpha_names, scale_names))

	n_parameters <- n_psi + n_phi + n_alpha + n_gamma + n_scale
	n_vars <- list(n_psi = n_psi, n_phi = n_phi, n_alpha = n_alpha, n_gamma = n_gamma, n_scale = n_scale)

	if (stan_data$K > 1){
		n_vars <- lapply(n_vars, function(x){x* stan_data$K})
		all_names <- GenClassNames(parm_names$all_names, stan_data$K)

		if (stan_data$single_scale == 1){
			n_vars$n_scale = 1
			all_names <- c(all_names[1:(length(all_names)-stan_data$K)], scale_names)
		}
		# membership names
		n_vars$n_beta <- stan_data$L * (stan_data$K - 1)
		delta.names <- GenClassNames(colnames(stan_data[["data_class"]]), stan_data$K)
		parm_names$delta.names <- grep("class1", delta.names, invert=TRUE, value = TRUE)

		parm_names$all_names <- c(all_names, parm_names$delta.names)
	}

	#standard deviations
	if(random_parameters != "fixed"){
		if(stan_data$gamma_nonrandom == 1){
			n_gamma_rp <- 0
			gamma_sd_names <- NULL
		} else if(stan_data$gamma_nonrandom == 0){
			n_gamma_rp <- n_gamma
			gamma_sd_names <- gamma_names
		}

		if(stan_data$alpha_nonrandom == 1){
			n_alpha_rp <- 0
			alpha_sd_names <- NULL
		} else if(stan_data$alpha_nonrandom == 0){
			n_alpha_rp <- n_alpha
			alpha_sd_names <- alpha_names
		}

		n_vars$n_std_dev <- n_psi + n_phi + n_gamma_rp + n_alpha_rp
		parm_names$sd_names <- paste0("sd.", c(psi_names, phi_names, gamma_sd_names, alpha_sd_names))

		if (random_parameters == "corr"){
			n_vars$n_std_corr <- n_vars$n_std_dev * (n_vars$n_std_dev-1) / 2
		}
	}

	n_vars$n_parms_total <- Reduce("+",n_vars)

	parms_info <- list(n_vars = n_vars,
					   parm_names = parm_names,
					   alt_names = alt_names)

	return(parms_info)
}


#' @title GenClassNames
#' @description Create names for LC
#' @param parms_names list of parameter names
#' @inheritParams mdcev
#' @return A vector of LC names
#' @noRd
GenClassNames <- function(parms_names, n_classes){
	if(length(parms_names) > 0){
		classes <- paste0(rep('class', n_classes), sep = "",c(1:n_classes))
		parms_names <- paste0(rep(classes, length(parms_names)), sep=".", rep(parms_names,  each = n_classes))
	} else {
		parms_names <- NULL
	}
	return(parms_names)
}
