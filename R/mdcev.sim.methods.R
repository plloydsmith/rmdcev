############################
# S3 method for mdcev.sim
#############################

#' @rdname mdcev.sim
#' @method print mdcev.sim
#' @param x,object an object of class `mdcev.sim`
#' @param digits the number of digits,
#' @param width the width of the printing,
#' @export
print.mdcev.sim <- function(x, digits = max(3, getOption("digits") - 3),
							width = getOption("width"), ...){

	cat("\nSize of simulations:\n")
	print.default(dim(x))
	cat("\n")
	invisible(x)
}

#' @rdname mdcev.sim
#' @method summary mdcev.sim
#' @param ci choose confidence interval for simulations. Default is 95 percent.
#' @export
summary.mdcev.sim <- function(object, ci = 0.95, ...){

	# check if list then demand
	if(is.list(object[[1]])){
		nobs <- length(object)
		nsims <- length(object[[1]])
		nalts <- ncol(object[[1]][[1]])
		npols <- nrow(object[[1]][[1]])

		out <- tibble::tibble(demand = unlist(object),
							  id = rep(1:nobs, each = nsims*nalts*npols),
							  alt = rep(0:(nalts-1), each = npols, times = nsims*nobs),
							  policy = rep(paste0(rep("policy",npols), 1:npols), times =nobs*nsims*nalts ),
							  sim_id = rep(1:nsims, each = nalts*npols, times = nobs)) %>%
			dplyr::group_by(.data$alt, .data$policy, .data$sim_id) %>%
			dplyr::summarise(demand = mean(demand)) %>%
			dplyr::group_by(.data$policy, .data$alt) %>%
			dplyr::summarise(mean = round(mean(demand),2),
							 std_dev = round(stats::sd(demand),2),
							 ci_lo = round(stats::quantile(demand, (1-ci)/2),2),
							 ci_hi = round(stats::quantile(demand, ci+(1-ci)/2),2))

		colnames(out) <- c("policy", "alt", "mean", "std.dev",
						   paste0("ci_lo",(1-ci)/2*100, "%"),
						   paste0("ci_hi",(ci+(1-ci)/2)*100, "%"))
	}
	else{
		out <- apply(simplify2array(object),1:2, mean)
		colnames(out)<- paste0(rep("policy",ncol(out)), 1:ncol(out))

		out <- tbl_df(out) %>%
			tidyr::gather(policy, wtp) %>%
			dplyr::group_by(.data$policy) %>%
			dplyr::summarise(mean = mean(wtp),
							 std_dev = stats::sd(wtp),
							 ci_lo = stats::quantile(wtp, (1-ci)/2),
							 ci_hi = stats::quantile(wtp, ci+(1-ci)/2))

		colnames(out) <- c("policy", "mean", "std.dev",
						   paste0("ci_lo",(1-ci)/2*100, "%"),
						   paste0("ci_hi",(ci+(1-ci)/2)*100, "%"))

	}

	object$CoefTable    <- out
	class(object)       <- c("summary.mdcev.sim", "mdcev.sim")
	return(object)
}

#' @rdname mdcev.sim
#' @method print summary.mdcev.sim
#' @export
print.summary.mdcev.sim <- function(x, digits = max(3, getOption("digits") - 2),
									width = getOption("width"),
									...){
	print(x$CoefTable, digits = digits)

	invisible(x)
}
