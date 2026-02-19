# Benchmark MLE estimation speed across model types and data sizes
# Usage: source("inst/benchmarks/bench_estimation.R")

library(bench)
library(rmdcev)

data(data_rec, package = "rmdcev")

prep_data <- function(n) {
  mdcev.data(
    data_rec,
    subset = id <= n,
    id.var = "id",
    alt.var = "alt",
    choice = "quant"
  )
}

# ---- Standard models (hybrid0, gamma, alpha) ----
results_mle <- bench::press(
  nobs = c(500, 1000),
  model = c("hybrid0", "gamma", "alpha"),
  {
    d <- prep_data(nobs)
    init_list <- switch(
      model,
      "hybrid0" = list(scale = 1),
      "gamma" = list(),
      "alpha" = list(),
      stop("Unknown model")
    )
    bench::mark(
      mdcev(
        ~0,
        data = d,
        model = model,
        algorithm = "MLE",
        print_iterations = FALSE,
        init = init_list
      ),
      iterations = 3,
      check = FALSE
    )
  }
)

print(results_mle[, c("nobs", "model", "median", "mem_alloc")])

# ---- kt_ee model (requires phi covariates) ----
prep_data_ktee <- function(n) {
  d <- data_rec
  d$beach <- ifelse(d$alt == "beach", 1, 0)
  mdcev.data(
    d,
    subset = id <= n,
    id.var = "id",
    alt.var = "alt",
    choice = "quant"
  )
}

results_ktee <- bench::press(
  nobs = c(500, 1000),
  {
    d <- prep_data_ktee(nobs)
    bench::mark(
      mdcev(
        formula = ~ ageindex | 0 | beach,
        data = d,
        model = "kt_ee",
        gamma_ascs = 0,
        algorithm = "MLE",
        print_iterations = FALSE,
        prior_phi_sd = 1
      ),
      iterations = 3,
      check = FALSE
    )
  }
)

print(results_ktee[, c("nobs", "median", "mem_alloc")])
