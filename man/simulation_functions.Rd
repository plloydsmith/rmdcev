% Generated manually — do not overwrite with roxygen2
\name{simulation_functions}
\alias{simulation_functions}
\alias{DrawError_rng}
\alias{MarshallianDemand}
\alias{HicksianDemand}
\alias{ComputeUtilJ}
\alias{CalcWTP_rng}
\alias{CalcMarshallianDemand_rng}
\title{Low-level MDCEV simulation functions}
\description{
Compiled C++ routines that underlie \code{\link{mdcev.sim}} and
\code{\link{GenerateMDCEVData}}.  These functions are exported for
advanced users and test code; most users should call the higher-level
wrappers instead.

\describe{
  \item{\code{DrawError_rng}}{Draw preference error terms for one
    individual using the MLHS or independent draw approach.}
  \item{\code{MarshallianDemand}}{Compute Marshallian (uncompensated)
    demand for one individual given marginal utilities at zero.}
  \item{\code{HicksianDemand}}{Compute Hicksian (compensated) demand
    for one individual given a reference utility level.}
  \item{\code{ComputeUtilJ}}{Evaluate the KT utility function for one
    individual.}
  \item{\code{CalcWTP_rng}}{Simulate welfare (WTP) across
    \code{nerrs} error draws for one individual under a price policy.}
  \item{\code{CalcMarshallianDemand_rng}}{Simulate Marshallian demand
    across \code{nerrs} error draws for one individual under a price
    policy.}
}
}
\arguments{
  \item{base_rng__}{External pointer to the Boost RNG created by
    \code{\link{rmdcev_get_rng}}.}
  \item{pstream__}{External pointer to the output stream created by
    \code{\link{rmdcev_get_stream}}.}
  \item{model_num}{Integer model type: 1 = gamma, 2 = alpha,
    3 = hybrid, 4 = hybrid0, 5 = kt_ee, 6 = gamma1.}
  \item{income}{Individual income (scalar).}
  \item{quant_j}{Vector of non-numeraire quantities.}
  \item{quant_num}{Numeraire quantity consumed.}
  \item{price_j}{Vector of non-numeraire prices.}
  \item{price}{Full price vector (numeraire first).}
  \item{psi_j, psi}{Psi utility parameters.}
  \item{phi_j, phi}{Phi parameters (kt_ee model).}
  \item{gamma_j, gamma}{Gamma satiation parameters.}
  \item{alpha}{Alpha parameters.}
  \item{scale}{Scale parameter.}
  \item{nalts}{Number of non-numeraire alternatives.}
  \item{nerrs}{Number of error draws.}
  \item{cond_error}{1 = conditional error draws; 0 = unconditional.}
  \item{draw_mlhs}{1 = Modified Latin Hypercube Sampling; 0 = iid.}
  \item{algo_gen}{Algorithm: 0 = hybrid (fast, requires alpha = 0);
    1 = general bisection.}
  \item{MUzero}{Vector of marginal utilities at zero consumption.}
  \item{util}{Reference utility level (Hicksian demand).}
  \item{tol}{Convergence tolerance for bisection.}
  \item{max_loop}{Maximum bisection iterations.}
}
\seealso{
  \code{\link{mdcev.sim}}, \code{\link{PrepareSimulationData}},
  \code{\link{rmdcev_get_rng}}, \code{\link{rmdcev_get_stream}}
}
