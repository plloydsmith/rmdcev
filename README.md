---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# rmdcev

The rmdcev R package estimates and simulates multiple discrete-continuous extreme value (MDCEV) demand models with observed and unobserved individual heterogneity. Fixed parameter, latent class, and random parameter models can be estimated. These models are estimated using maximum likelihood or Bayesian estimation techniques and are implemented in Stan, which is a C++ package for performing full Bayesian inference (see http://mc-stan.org/). The package also supports Phaneuf and von Haefen's Kuhn-Tucker model specification.

## Installation

You can install the released version of rmdcev from GitHub

``` r
install_github("plloydsmith/rmdcev")
```

Development is in progress. Currently users can estimate the following models:

1. Fixed parameter models
2. Latent class models

Under development
1. Random parameters models
2. Phaneuf/von Haefen Kuhn Tucker model

For more details on the model specification and estimation:

Bhat, C.R. (2008) ["The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of Utility Function Parameters, Identification Considerations, and Model Extensions."](https://www.sciencedirect.com/science/article/pii/S0191261507000677) Transportation Research Part B, 42(3): 274-303. 

For more details on the welfare simulation:

Lloyd-Smith, P (2018). ["A New Approach to Calculating Welfare Measures in Kuhn-Tucker Demand Models."](https://www.sciencedirect.com/science/article/pii/S1755534517300994) Journal of Choice Modeling, 26: 19-27
