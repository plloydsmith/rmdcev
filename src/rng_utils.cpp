// rng_utils.cpp - Provides boost::ecuyer1988 RNG and output stream
// for use by Stan simulation functions compiled with StanHeaders.
// This replaces rstan::get_rng() and rstan::get_stream() so that
// the package does not require rstan at runtime.

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]

#include <stan/math/prim/fun/Eigen.hpp>
#include <Rcpp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/ecuyer1988.hpp>

// [[Rcpp::export]]
SEXP rmdcev_get_rng(int seed = 0) {
    boost::ecuyer1988* rng = new boost::ecuyer1988(static_cast<uint32_t>(seed));
    return Rcpp::XPtr<boost::ecuyer1988>(rng, true);
}

// [[Rcpp::export]]
SEXP rmdcev_get_stream() {
    std::ostream* o = &Rcpp::Rcout;
    return Rcpp::XPtr<std::ostream>(o, false);
}
