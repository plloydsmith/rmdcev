#' @title mdcev.datacheck
#' @description Check mdcev data
#' @inheritParams mdcev
#' @noRd
mdcev.datacheck <- function(data_input){

	message("Checking data...")

	id.name <- attr(data_input, "id")
	quant.name <- attr(data_input, "choice")
	price.name <- attr(data_input, "price")
	income.name <- attr(data_input, "income")

	if(!id.name %in% colnames(data_input))
		stop("Data must have id column for individual")

	if(!quant.name %in% colnames(data_input))
		stop("Data must have choice column for consumption")

	if(!price.name %in% colnames(data_input))
		stop("Data must have price column for non-numeraire alternatives")

	if(!income.name %in% colnames(data_input))
		stop("Data must have income column for individual's income")

	data_input$expend_alt <- data_input[[price.name]] * data_input[[quant.name]]
	id.expend <- stats::aggregate(expend_alt ~ get(id.name), data = data_input, FUN = sum )
	id.income <- data_input[!duplicated(data_input[,c(id.name,income.name)]),]
	id.expend$expend_numeraire <- id.income[[income.name]] - id.expend$expend_alt
	id.expend$expend_alt <- NULL

	if (sum(id.expend$expend_numeraire < 0) > 0){
		print(paste0("The outside good is not chosen for individuals: ",
					 toString(id.expend[which(id.expend$expend_numeraire < 0), ]),
			  "\n Ensure that the total expenditure on the alternatives is less than income."))
		stop()
	}

	if (sum(data_input$price <= 0) > 0){
		print(paste0("Price is less than or equal to 0 for at least one individual alt in rows: ",
					 toString(as.character(which(data_input$price <= 0)))))
		stop()
	}
	message("Data is good")
}
