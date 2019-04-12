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

	if(!"inc" %in% colnames(data))
		stop("Data must have inc column for individual's income")

	check <- tbl_df(data) %>%
		mutate(expend = price * quant) %>%
		group_by(id) %>%
		summarise(numeraire = mean(inc) - sum(expend)) %>%
		select(numeraire)

	if (sum(check$numeraire < 0) > 0){
		print(paste0("Numeraire is less than 0 for individuals in rows: ",
					 toString(as.character(which(check$numeraire < 0)))))
		stop()
}
	if (sum(data$price <= 0) > 0){
		print(paste0("Price is less than or equal to 0 for at least one individual good in rows: ",
					 toString(as.character(which(data$price <= 0)))))
}
	message("Data is good")
}
