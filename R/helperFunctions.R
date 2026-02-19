#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @return A list
#' @noRd
CreateListsRow <- function(x){
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CleanInit
#' @description Adds appropriate dimensions for initial values
#' @param init_input cleaned version of inits
#' @noRd
CleanInit <- function(init_input){
	if (!is.list(init_input)){
		temp = init_input
	} else{

	temp = list(scale = NULL)

	# Add dimension to starting values
	temp <- lapply(init_input, function(x){
		x <- matrix(x, nrow = 1, length(x))
	})
	if(!is.null(init_input$scale))
		temp$scale <- array(init_input$scale, dim = 1)
	}
	return(temp)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @noRd
CreateListsCol <- function(x){
	if (is.vector(x))
		out <- as.list(x)
	else
		out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @param x matrix B to be multiplied
#' @param mat_temp matrix A to be multiplied
#' @param n_rows number of rows for final matrix
#' @noRd
MultiplyMatrix <- function(mat_temp, x, n_rows){
	out <- matrix(mat_temp %*% x , nrow = n_rows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @param x list to be converted
#' @noRd
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
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
CreateBlankPolicies <- function(npols, model, price_change_only = TRUE){

	price_p <- CreateListsRow(matrix(0, nrow = npols, ncol = model[["stan_data"]][["J"]] + 1))
	model_num <- model[["stan_data"]][["model_num"]]

	dat_psi_p <- NULL
	dat_phi_p <- NULL

	if (price_change_only == FALSE){
		if (model_num < 5 && model[["parms_info"]][["n_vars"]][["n_psi"]] == 0)
			stop("No psi variables to vary! Use price_change_only == TRUE option.")

		if (model_num == 5 && model[["parms_info"]][["n_vars"]][["n_phi"]] == 0)
			stop("No phi variables to vary! Use price_change_only == TRUE option.")

		if (model_num < 5)
			dat_psi_p <- lapply(seq_len(npols), function(X) model[["stan_data"]][["dat_psi"]])
		else if (model_num == 5)
			dat_phi_p <- lapply(seq_len(npols), function(X) model[["stan_data"]][["dat_phi"]])
	}

	out <- list(price_p = price_p, dat_psi_p = dat_psi_p, dat_phi_p = dat_phi_p,
				price_change_only = price_change_only)
	return(out)
}

#' @title CreatePsiMatrix
#' @param psi_j matrix (JXn_psi_j) of alt-specific attributes
#' @param psi_i matrix (IXn_psi_i) of people-specific attributes
#' @description Creates the Psi data matrix for use in mdcev model
#' @noRd
CreatePsiMatrix <- function(psi_j = NULL, psi_i = NULL){
	if(!is.na(psi_i))
		psi_i <- lapply(psi_i, function(x) {rep(x, each= nrow(psi_j))})
	if(!is.na(psi_j))
		psi_j <- lapply(psi_j, function(x) {rep(x, times=nrow(psi_i))})

	dat_psi <- c(psi_j, psi_i)
return(dat_psi)
}

#' @title GrabParms
#' @param data est_par object from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @noRd
GrabParms <- function(data, parm_name){
	out <- data %>%
		dplyr::filter(grepl(parm_name, .data$parms)) %>%
		tidyr::pivot_wider(id_cols = sim_id, names_from = "parms", values_from = "value") %>%
		dplyr::select(-sim_id) %>%
		as.matrix(.)
	return(out)
}


#' @title GrabIndividualParms
#' @param est_sim est_sim from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @noRd
GrabIndividualParms <- function(est_sim, parm_name){
	out <- est_sim %>%
		dplyr::filter(grepl(c(parm_name), parms)) %>%
		dplyr::select(id, sim_id, "parm_id", beta) %>%
		tidyr::spread("parm_id", beta) %>%
		dplyr::select(-sim_id) %>%
		dplyr::group_split(id, .keep = F)
	return(out)
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
            dplyr::rename(parms = variable, n_eff = ess_bulk, Rhat = rhat)
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

#' @title CreatePsi
#' @param dat_vars_i psi data for each person
#' @param est_pars_i psi parameter estimates for each person
#' @description Works at individual level and creates the psi
#' variables for each simulation, policy, alternative
#' @noRd
CreatePsi = function(dat_vars_i, est_pars_i, J, NPsi_ij, psi_ascs, npols){
#	dat_vars_i = dat_vars[[3]]
#	est_pars_i = psi_sim_temp[[3]]
	lpsi = matrix(0, nrow(est_pars_i), J)
	if (psi_ascs == 1){
		psi_non_ascs_start = J
		psi_non_ascs_end = J+NPsi_ij-1
		if (nrow(est_pars_i) == 1)
			lpsi = lpsi + c(0, est_pars_i[,1:(J-1)])
		else
			lpsi = lpsi + cbind(0, est_pars_i[,1:(J-1)]) ##  alternative specific constants
	} else if (psi_ascs == 0){
		psi_non_ascs_start = 1
		psi_non_ascs_end = NPsi_ij
	}

	if ((NPsi_ij > 0) && (nrow(dat_vars_i) == J)){
		psi_non_ascs = est_pars_i[,psi_non_ascs_start:psi_non_ascs_end, drop=FALSE]
		lpsi = lpsi + as.matrix(psi_non_ascs) %*% t(as.matrix(dat_vars_i))
	}

	if (nrow(dat_vars_i) > J){
		if (NPsi_ij == 0){
			lpsi = CreateListsRow(lpsi)
			lpsi = lapply(lpsi, function(x){
				lpsi= matrix(x,nrow=npols,ncol=length(x),byrow=TRUE)
			})
		} else if (NPsi_ij > 0){
			dat_vars_i = dat_vars_i %>%
				group_split(policy, .keep = F)
			psi_non_ascs <- est_pars_i[,psi_non_ascs_start:psi_non_ascs_end, drop=FALSE]
			lpsi <- lapply(dat_vars_i, function(xx){
				lpsi = lpsi + as.matrix(psi_non_ascs) %*% t(as.matrix(xx))
				return(lpsi)
			})
			lpsi = aperm(array(unlist(lpsi),
							   dim = c(nrow(lpsi[[1]]), ncol(lpsi[[1]]), length(lpsi))),
						 perm=c(1,3,2))
			lpsi = lapply(seq(dim(lpsi)[1]), function(x) lpsi[ x, , ])
		}
	}
	return(lpsi)
}
