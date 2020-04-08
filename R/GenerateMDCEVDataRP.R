#' @title GenerateMDCEVDataRP
#' @description Simulate random parameter data for MDCEV model
#' @inheritParams GenerateMDCEVData
#' @param corr Whether to draw correlated random parameters (=1) or uncorrelated (=0)
#' @return A `mdcev.data` object, which is a `data.frame` in long
#'     format. Also includes parms_true with parameter values
#' @export
#' @examples
#' \donttest{
#' data <- GenerateMDCEVDataRP(model = "hybrid0")
#'}
GenerateMDCEVDataRP <- function(model,
								nobs = 1000,
								nalts = 10,
								inc_lo = 100000,
								inc_hi = 150000,
								price_lo = 100,
								price_hi = 500,
								alpha_parms = 0.5,
								scale_parms = 1,
								gamma_parms = stats::runif(nalts, 1, 10),
								psi_j_parms = c(-5, 0.5, 2),
								nerrs = 1,
								tol = 1e-20,
								max_loop = 999,
								corr = 0){

	income <-  stats::runif(nobs, inc_lo, inc_hi)
	price <- matrix(stats::runif(nobs*nalts, price_lo, price_hi), nobs, nalts)

	num_psi <- length(psi_j_parms)
	RP <- num_psi + nalts + 1

	if (model == "gamma"){
		model_num <- 1
		algo_gen <- 1
		a <- c(rep(-Inf, num_psi), rep(0.01, RP-num_psi))
		b <- c(rep(Inf, num_psi), rep(Inf, RP-num_psi-1), .99)
	} else if (model == "alpha"){
		model_num <- 2
		alpha_parms <- 0 + stats::runif(nalts+1, 0.01, .98)
		algo_gen <- 1
		a <- c(rep(-Inf, num_psi), rep(0.01, RP-num_psi))
		b <- c(rep(Inf, num_psi), rep(.99, RP-num_psi))
	} else if (model == "hybrid"){
		model_num <- 3
		algo_gen <- 0
		a <- c(rep(-Inf, num_psi), rep(0.01, RP-num_psi))
		b <- c(rep(Inf, num_psi), rep(Inf, RP-num_psi-1), .99)
	} else if (model == "hybrid0"){
		RP <- RP - 1 # subtract one for fixed alpha
		model_num <- 4
		algo_gen <- 0
		alpha_parms <- NULL
		a <- c(rep(-Inf, num_psi), rep(0.01, RP-num_psi))
		b <- c(rep(Inf, num_psi), rep(Inf, RP-num_psi))
	} else
		stop("No model specificied. Choose a model")

	beta <- c(psi_j_parms, gamma_parms, alpha_parms)

	# Random correlation matrix to construct a covariance matrix of the individual betas
	tau <- stats::rgamma(RP, 4, 2)
	if (corr == 0){
		Sigma <- diag(tau)
	} else if (corr == 1){
		Sigma <- diag(tau) %*% stats::cor(matrix(stats::rnorm(RP*(RP+1)), RP+1, RP)) %*% diag(tau)
	}


	beta_individual <- tmvtnorm::rtmvnorm(n=nobs, mean=beta, sigma=Sigma, lower=a, upper=b)

	indexes <- tibble(individual = rep(1:nobs, each = nalts),
					  task = rep(1:nobs, each = nalts),
					  row = 1:(nobs*nalts)) %>%
		group_by(task) %>%
		summarise(task_individual = first(individual),
				  start = first(row),
				  end = last(row))

	# Create psi variables that vary over alternatives
	psi_j <- cbind(rep(1,nalts), # add constant term
				   matrix(stats::runif(nalts*(num_psi-1), 0 , 1), nrow = nalts))
	psi_j <-  rep(1, nobs) %x% psi_j

	psi_beta <- beta_individual[,1:num_psi] %x% rep(1, nalts)

	psi_sims <- matrix(rowSums(psi_j * psi_beta), ncol = nalts, byrow = TRUE)
	psi_sims <- CreateListsRow(psi_sims)
	psi_sims <- list(psi_sims )
	names(psi_sims) <- "psi_sims"


	# Create gamma variables that vary over alternatives
	#	dat_gamma <- rep(1, nobs) %x% diag(nalts)
	#	gamma_sims <- matrix(rowSums(dat_gamma * gamma_beta), ncol = nalts, byrow = TRUE)

	gamma_sims <- beta_individual[,(num_psi+1):(num_psi+nalts)] #%x% rep(1, nalts)

	if (model == "gamma"){
		alpha_sims <- cbind(beta_individual[,RP], matrix(0, nobs, nalts))
	} else if (model == "alpha"){
		alpha_sims <- beta_individual[,(num_psi+nalts+1):RP]
		gamma_sims <- matrix(1, nobs, nalts)
	} else if (model == "hybrid"){
		alpha_sims <- matrix(rep(beta_individual[,(num_psi+nalts+1):RP],nalts+1), nrow = nobs, ncol = nalts+1, byrow = F)
	} else if (model == "hybrid0"){
		alpha_sims <- matrix(1e-6, nobs, nalts+1)
	} else
		stop("No model specificied. Choose a model")


	gamma_sims <- CreateListsRow(gamma_sims)
	gamma_sims <- list(gamma_sims )
	names(gamma_sims) <- "gamma_sims"

	alpha_sims <- CreateListsRow(alpha_sims)
	alpha_sims <- list(alpha_sims )
	names(alpha_sims) <- "alpha_sims"

	income_list <- list(as.list(income))
	names(income_list) <- "income" # price normalized MU at zero

	price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
	price_list <- list(CreateListsRow(price_list))
	names(price_list) <- "price" # price normalized MU at zero

	df_indiv <- c(income_list, price_list, psi_sims, gamma_sims, alpha_sims)


	PRNG <-rstan::get_rng(seed = 3)
	o <- rstan::get_stream()

	quant <- purrr::pmap(df_indiv, CalcmdemandOne_rng,
				  scale_sim=scale_parms,
				  nerrs=nerrs, algo_gen = algo_gen,
				  tol = tol, max_loop = max_loop,
				  PRNG, o)

	# Convert simulated data into estimation data
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)
	quant <- quant[,2:(ncol(quant))]
	quant <- as.vector(t(quant))
	price <- as.vector(t(price))

	id <- rep(1:nobs, each = nalts)
	alt <- rep(1:nalts, times = nobs)
	income <- rep(income, each = nalts)

	data <- as.data.frame(cbind(id, alt, quant, price, psi_j, income))

	data <- mdcev.data(data,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

	parms_true <- df_indiv
	parms_true$income <- parms_true$price <-  NULL

	out <- list(data = data,
				parms_true = beta_individual,
				Sigma = Sigma,
				beta = beta)
	return(out)
}

