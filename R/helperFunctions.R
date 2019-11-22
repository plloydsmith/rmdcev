#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @return A list
#' @export
#' @examples
#'
#' tmp <- matrix(0, nrow = 10, ncol = 5)
#' tmp_list <- CreateListsRow(tmp)
#'
CreateListsRow <- function(x){
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @export
#' @examples
#'
#' tmp <- matrix(0, nrow = 10, ncol = 5)
#' tmp_list <- CreateListsCol(tmp)
#'
CreateListsCol <- function(x){
	out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @param x matrix B to be multiplied
#' @param mat_temp matrix A to be multipled
#' @param n_rows number of rows for final matrix
#' @keywords internal
MultiplyMatrix <- function(mat_temp, x, n_rows){
	out <- matrix(mat_temp %*% x , nrow = n_rows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @param x list to be converted
#' @keywords internal
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
#' @param price_change_only Logical value for whether to include policy changes to dat_psi.
#' TRUE implies that only price changes are used in simulation.
#' @param npols Number of policies to simulate
#' @param nalts Number of non-numeraire alts
#' @param dat_psi Psi data matrix used in estimation
#' @export
#' @examples
#' \donttest{
#' CreateBlankPolicies(npols = 2, nalts = 10, dat_psi = NULL, price_change_only = TRUE)
#'}
CreateBlankPolicies <- function(npols, nalts, dat_psi, price_change_only){

	price_p <- CreateListsRow(matrix(0, nrow = npols, ncol = nalts + 1))

	if (price_change_only == FALSE)
		dat_psi_p <- lapply(seq_len(npols), function(X) dat_psi)
	else if (price_change_only == TRUE)
		dat_psi_p <- NULL

	out <- list(price_p = price_p, dat_psi_p = dat_psi_p,
				price_change_only = price_change_only)
	return(out)
}

#' @title CreatePsiMatrix
#' @param psi_j matrix (JXn_psi_j) of alt-specific attributes
#' @param psi_i matrix (IXn_psi_i) of people-specific attributes
#' @description Creates the Psi data matrix for use in mdcev model
#' @keywords internal
CreatePsiMatrix <- function(psi_j = NULL, psi_i = NULL){
	if(!is.na(psi_i))
		psi_i <- purrr::map(psi_i, function(x) {rep(x, each= nrow(psi_j))})
	if(!is.na(psi_j))
		psi_j <- purrr::map(psi_j, function(x) {rep(x, times=nrow(psi_i))})

	dat_psi <- c(psi_j, psi_i)
return(dat_psi)
}

#' @title GrabParms
#' @param data est_par object from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @keywords internal
GrabParms <- function(data, parm_name){
	out <- data %>%
		dplyr::filter(stringr::str_detect(.data$parms, parm_name)) %>%
		tidyr::spread_(key_col = 'sim_id',
				value_col = 'value') %>%
		dplyr::select(-parms) %>%
		as.matrix(.)
	return(out)
}
