#' The 'rmdcev' package.
#'
#' @description The rmdcev R package estimates and simulates multiple discrete-continuous extreme value (MDCEV)
#' demand models with observed and unobserved individual heterogneity. Fixed parameter, latent class,
#' and random parameter models can be estimated. These models are estimated using maximum likelihood or
#' Bayesian estimation techniques and are implemented in Stan, which is a C++ package for performing full
#' Bayesian inference (see http://mc-stan.org/). The package also supports Phaneuf and von Haefen's Kuhn-Tucker
#' model specification.
#'
#' @docType package
#' @name rmdcev
#' @aliases rmdcev
#' @useDynLib rmdcev, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")
	utils::globalVariables(c(".", "quant", "lp__", "good", "data",
							"task", "individual", "CalcWTP_rng",
							"CalcMarshallianDemand_rng", "CalcWTPPriceOnly_rng",
							"CalcmdemandOne_rng", "parms", "sim_id", "value", "policy", "std_dev",
							"names_b", "parm_num", "n_eff", "Rhat"))
