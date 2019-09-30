#' data.frame for mdcev model
#'
#' shape a `data.frame` in a suitable form for the use of the
#' `mdcev` function and complete some data checks
#'
#' @name mdcev.data
#' @param data a `data.frame`,
#' @param id.var the name of the variable that contains the individual
#'     index.
#' @param alt.var the name of the variable that contains the
#'     alternative index or the name under which the alternative index
#'     will be stored (the default  name is `alt`),
#' @param choice the variable indicating the consumption of non-numeraire
#'     alternatives that is made: it has to be a numerical vector
#' @param price the variable indicating the price of the non-numeraire
#'     alternatives. Default is "price"
#' @param income the variable indicating the income of the individual.
#'     Default is "income".
#' @param alt.levels the name of the alternatives: if null,
#'     they are guessed from the `alt.var` argument,
#' @param drop.index should the index variables be dropped from the
#'     `data.frame`,
#' @param subset a logical expression which defines the subset of
#'     observations to be selected,
#' @param ... further arguments.
#' @return A `mdcev.data` object, which is a `data.frame` in long
#'     format, *i.e.* one line for each alternative.  It has a `index`
#'     attribute, which is a `data.frame` that contains the
#'     index of the individual (`id`) and the index of the alternative
#'     (`alt`).
#' @export
mdcev.data <- function(data,
					   id.var = NULL,
					   alt.var = NULL,
					   choice = NULL,
					   price = "price",
					   income = "income",
					   alt.levels = NULL,
					   drop.index = FALSE,
					   subset = NULL, ...){
	# id, alt : the index variables

	# if a subset argument is provided, subset the original data frame
	cldata <- match.call(expand.dots = TRUE)
	if (match("subset", names(cldata), 0)){
		m <- match(c("data", "subset"), names(cldata), 0)
		cldata <- cldata[c(1, m)]
		names(cldata)[2] <- "x"
		cldata[[1]] <- as.name("subset")
		data <- eval(cldata, parent.frame())
	}

	if (is.null(alt.var) && is.null(alt.levels))
		stop("at least one of alt.var and alt.levels should be filled")

	if (! is.null(alt.levels)){
		J <- length(alt.levels)
		n <- nrow(data) / J
		alt <- factor(rep(alt.levels, n), levels = alt.levels)
		if (!is.null(alt.var) && !is.null(data[[alt.var]])){
			warning(paste("variable", alt.var, "exists and will be replaced"))
			alt.is.variable <- TRUE
		}
		else alt.is.variable <- FALSE
		alt.name <- ifelse(is.null(alt.var), "alt", alt.var)
	}
	else{
		alt.name <- alt.var
		alt.is.variable <- TRUE
		if (!is.factor(data[[alt.name]])) data[[alt.name]] <- factor(data[[alt.name]])
		alt.levels <- levels(data[[alt.name]])
		J <- length(alt.levels)
		alt <- data[[alt.name]]
	}
	n <- nrow(data) / J
	if (! is.null(id.var)){
		idpos <- which(names(data) == id.var)
		id.var <- as.factor(data[[id.var]])
	}

	if (! is.null(choice)){
		choice.name <- choice
		choice <- data[[choice]]
		data[[choice.name]] <- choice
	}
	id <- as.factor(id.var)
	alt <- as.factor(alt)
	row.names(data) <- paste(id, alt, sep = ".")

	altpos <- which(names(data) == alt.name)

	if (drop.index){
		if (! is.null(id.var)) data <- data[, -c(altpos, idpos)]
		else data <- data[, -c(altpos)]
	}

	index <- data.frame(id = id, alt = alt)
	rownames(index) <- rownames(data)
	attr(data, "index") <- index
	attr(data, "class") <- c("mdcev.data", "data.frame")
	if (! is.null(choice)) attr(data, "choice") <- choice.name
	attr(data, "price") <- price
	attr(data, "income") <- income

	mdcev.datacheck(data)
	return(data)
}
