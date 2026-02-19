# Benchmark mdcev.sim() speed across error draws, sim type, and nsims
# Usage: source("inst/benchmarks/bench_simulation.R")

library(bench)
library(rmdcev)

data(data_rec, package = "rmdcev")

# Fit one model outside the benchmark loop
d <- mdcev.data(
  data_rec,
  subset = id <= 100,
  id.var = "id",
  alt.var = "alt",
  choice = "quant"
)

result <- mdcev(
  ~0,
  data = d,
  model = "hybrid0",
  algorithm = "MLE",
  print_iterations = FALSE
)

npols <- 2
policies <- CreateBlankPolicies(npols, result, price_change_only = TRUE)
df_sim <- PrepareSimulationData(result, policies, nsims = 5)

# ---- nerrs x sim_type grid ----
results_sim <- bench::press(
  nerrs = c(100, 500),
  sim_type = c("welfare", "demand"),
  {
    bench::mark(
      mdcev.sim(
        df_sim$df_indiv,
        df_common = df_sim$df_common,
        sim_options = df_sim$sim_options,
        cond_err = 1,
        nerrs = nerrs,
        sim_type = sim_type
      ),
      iterations = 5,
      check = FALSE
    )
  }
)

print(results_sim[, c("nerrs", "sim_type", "median", "mem_alloc")])

# ---- nsims (number of parameter draws) grid ----
results_nsims <- bench::press(
  nsims = c(20, 50),
  {
    df_s <- PrepareSimulationData(result, policies, nsims = nsims)
    bench::mark(
      mdcev.sim(
        df_s$df_indiv,
        df_common = df_s$df_common,
        sim_options = df_s$sim_options,
        cond_err = 1,
        nerrs = 10,
        sim_type = "welfare"
      ),
      iterations = 5,
      check = FALSE
    )
  }
)

print(results_nsims[, c("nsims", "median", "mem_alloc")])

# ---- conditional vs unconditional errors ----
results_cond <- bench::press(
  cond_err = c(0L, 1L),
  {
    bench::mark(
      mdcev.sim(
        df_sim$df_indiv,
        df_common = df_sim$df_common,
        sim_options = df_sim$sim_options,
        cond_err = cond_err,
        nerrs = 50,
        sim_type = "welfare"
      ),
      iterations = 5,
      check = FALSE
    )
  }
)

print(results_cond[, c("cond_err", "median", "mem_alloc")])
