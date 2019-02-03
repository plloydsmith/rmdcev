#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @return A list
#' @export
#'
CreateListsRow <- function(x){
	#
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @export
CreateListsCol <- function(x){
	#
	out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @export
MultiplyMatrix <- function(mat_temp, x, nrows){
	# m - matrix
	# x - vector
	# nrows - number of rows in m
	out <- matrix(mat_temp %*% x , nrow = nrows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @export
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
#' @export
CreateBlankPolicies <- function(npols, ngoods, dat_psi){
	price_p <- CreateListsRow(matrix(0, nrow = npols, ncol = ngoods + 1))
	dat_psi_p <- lapply(seq_len(npols), function(X) dat_psi)
	out <- list(price_p = price_p, dat_psi_p = dat_psi_p)
	return(out)
}
