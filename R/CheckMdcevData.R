#' @title CheckMdcevData
#' @description Check MDCEV data
#' @inheritParams FitMDCEV
#' @export
CheckMdcevData <- function(data){

	message("Checking data...")

	if(!"id" %in% colnames(data))
		stop("Data must have id column for individual")

	if(!"good" %in% colnames(data))
		stop("Data must have good column for non-numeraire alternatives")

	if(!"quant" %in% colnames(data))
		stop("Data must have quant column for consumption")

	if(!"price" %in% colnames(data))
		stop("Data must have price column for non-numeraire alternatives")

	if(!"income" %in% colnames(data))
		stop("Data must have income column for individual's income")

	check <- tbl_df(data) %>%
		mutate(expend = .data$price * .data$quant) %>%
		group_by(.data$id) %>%
		summarise(numeraire = mean(.data$income) - sum(.data$expend)) %>%
		select(.data$numeraire)

	if (sum(check$numeraire < 0) > 0){
		print(paste0("Numeraire is less than 0 for individuals in rows: ",
					 toString(as.character(which(check$numeraire < 0)))))
		stop()
	}

	check <- tbl_df(data) %>%
		arrange(id, good)

	if (identical(identical(check$id, data$id), FALSE) || identical(identical(check$good, data$good), FALSE) ){
		stop("Data must be arranged by id then good. Can use dplyr command arrange(id, good).")
	}

	if (sum(data$price <= 0) > 0){
		print(paste0("Price is less than or equal to 0 for at least one individual good in rows: ",
					 toString(as.character(which(data$price <= 0)))))
		stop()
	}
	message("Data is good")
}
