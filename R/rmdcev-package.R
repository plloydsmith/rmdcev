#' rmdcev: Estimating and simulating Kuhn-Tucker and multiple discrete-continuous extreme value (MDCEV) demand models
#'
#' @description The rmdcev R package estimate and simulates Kuhn-Tucker demand models with individual
#' heterogeneity. The models supported by rmdcev are the multiple-discrete continuous extreme value (MDCEV)
#' model and Kuhn-Tucker specification common in the environmental economics literature on recreation demand.
#' Latent class and random parameters specifications can be implemented and the models are fit using maximum
#' likelihood estimation or Bayesian estimation. All models are implemented in Stan, which is a C++ package
#' for performing full Bayesian inference (see Stan Development Team, 2019) <https://mc-stan.org/>. The rmdcev package also implements
#' demand forecasting (Pinjari and Bhat (2011) <https://repositories.lib.utexas.edu/handle/2152/23880>) and
#' welfare calculation (Lloyd-Smith (2018) <doi:10.1016/j.jocm.2017.12.002>) for policy simulation.
#'
#' @name rmdcev
#' @aliases rmdcev
#' @useDynLib rmdcev, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#'
#' @keywords internal
"_PACKAGE"
#' @author Patrick Lloyd-Smith \email{patrick.lloydsmith@usask.ca}
#'
#' @references
#' Bhat, CR (2008). The multiple discrete-continuous extreme value (MDCEV) model: Role of utility function parameters, identification considerations, and model extensions. Transportation Research Part B: Methodological, 42(3): 274-303.\href{https://doi.org/10.1016/j.trb.2007.06.002}{(link)}

#' Lloyd-Smith, P (2018). A new approach to calculating welfare measures in Kuhn-Tucker demand models. Journal of Choice Modeling 26:19–27. \href{https://doi.org/10.1016/j.jocm.2017.12.002}{(link)}

#' Pinjari, AR, Bhat, CR (2011). Computationally Efficient Forecasting Procedures for Kuhn-Tucker Consumer Demand Model Systems: Application to Residential Energy Consumption Analysis. Department of Civil and Environmental Engineering, University of South Florida. \href{https://repositories.lib.utexas.edu/handle/2152/23880}{(link)}

#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. \href{https://mc-stan.org/}{(link)}
NULL

