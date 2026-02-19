# Benchmark core exported C++ functions directly
# Usage: source("inst/benchmarks/bench_cpp.R")

library(bench)
library(rmdcev)

data(data_rec, package = "rmdcev")

# ---- Setup: fit one model and extract parameters ----
d <- mdcev.data(data_rec, subset = id <= 100,
                id.var = "id", alt.var = "alt", choice = "quant")

result <- mdcev(~ 0, data = d, model = "hybrid0",
                algorithm = "MLE", print_iterations = FALSE)

nalts     <- result$stan_data[["J"]]
model_num <- result$stan_data[["model_num"]]
npols     <- 2
policies  <- CreateBlankPolicies(npols, result, price_change_only = TRUE)
df_sim    <- PrepareSimulationData(result, policies, nsims = 1)

income  <- df_sim[["df_indiv"]][["income"]][[1]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[1]]
price   <- df_sim[["df_indiv"]][["price"]][[1]]

quant_num <- as.numeric(income - quant_j %*% price[-1])

psi_j   <- result[["stan_data"]][["dat_psi"]][1:nalts, , drop = FALSE] %*%
            t(result[["stan_fit"]][["par"]][["psi"]])
phi_j   <- rep(0, nalts)
gamma_j <- result[["stan_fit"]][["par"]][["gamma"]]
gamma   <- c(1, gamma_j)
alpha   <- rep(0, nalts + 1)
scale   <- result[["stan_fit"]][["par"]][["scale"]]

tol_e    <- 1e-20
tol_l    <- 1e-20
max_loop <- 999

rng <- rmdcev_get_rng(seed = 42)
o   <- rmdcev_get_stream()

# Draw errors once for a fixed input to Marshallian/Hicskian benchmarks
error <- DrawError_rng(quant_num, quant_j, price[-1],
                       psi_j, phi_j, gamma_j, alpha, scale,
                       model_num = model_num, nalts = nalts,
                       nerrs = 1, cond_error = 1, draw_mlhs = 1,
                       rng, o)

psi_b_err <- exp(c(0, psi_j) + error[[1]])
MUzero_b  <- psi_b_err / price

util <- ComputeUtilJ(income, quant_j, price[-1],
                     psi_b_err, phi_j, gamma[-1], alpha,
                     nalts, model_num, o)
price_p  <- price + c(0.001, rep(1, nalts))
MUzero_p <- psi_b_err / price_p

# ---- MarshallianDemand: algo_gen 0 vs 1 ----
cat("\n--- MarshallianDemand: algo_gen 0 vs 1 ---\n")
results_md <- bench::mark(
  algo0 = MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
                             nalts, algo_gen = 0, model_num, tol_e, max_loop, o),
  algo1 = MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
                             nalts, algo_gen = 1, model_num, tol_e, max_loop, o),
  check = FALSE
)
print(results_md[, c("expression", "median", "mem_alloc", "n_itr")])

# ---- HicksianDemand: algo_gen 0 vs 1 ----
cat("\n--- HicksianDemand: algo_gen 0 vs 1 ---\n")
results_hd <- bench::mark(
  algo0 = HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
                          nalts, algo_gen = 0, model_num, tol_l, max_loop, o),
  algo1 = HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
                          nalts, algo_gen = 1, model_num, tol_l, max_loop, o),
  check = FALSE
)
print(results_hd[, c("expression", "median", "mem_alloc", "n_itr")])

# ---- DrawError_rng: varying nerrs ----
cat("\n--- DrawError_rng: nerrs grid (conditional) ---\n")
results_draw <- bench::press(
  nerrs = c(1, 10, 50),
  {
    rng2 <- rmdcev_get_rng(seed = 1)
    bench::mark(
      DrawError_rng(quant_num, quant_j, price[-1],
                    psi_j, phi_j, gamma_j, alpha, scale,
                    model_num = model_num, nalts = nalts,
                    nerrs = nerrs, cond_error = 1, draw_mlhs = 1,
                    rng2, o),
      iterations = 20, check = FALSE
    )
  }
)
print(results_draw[, c("nerrs", "median", "mem_alloc")])

# ---- CalcWTP_rng and CalcMarshallianDemand_rng: varying nerrs ----
cat("\n--- CalcWTP_rng: nerrs grid ---\n")
results_wtp <- bench::press(
  nerrs = c(1, 10, 50),
  {
    rng3 <- rmdcev_get_rng(seed = 2)
    bench::mark(
      CalcWTP_rng(income, price, price_p, psi_j, phi_j, gamma_j, alpha, scale,
                  model_num = model_num, nalts = nalts, nerrs = nerrs,
                  cond_error = 1, draw_mlhs = 1, algo_gen = 0,
                  tol_e = tol_e, tol_l = tol_l, max_loop = max_loop,
                  rng3, o),
      iterations = 10, check = FALSE
    )
  }
)
print(results_wtp[, c("nerrs", "median", "mem_alloc")])

cat("\n--- CalcMarshallianDemand_rng: nerrs grid ---\n")
results_cmd <- bench::press(
  nerrs = c(1, 10, 50),
  {
    rng4 <- rmdcev_get_rng(seed = 3)
    bench::mark(
      CalcMarshallianDemand_rng(income, price, psi_j, phi_j, gamma_j, alpha, scale,
                                model_num = model_num, nalts = nalts, nerrs = nerrs,
                                cond_error = 0, draw_mlhs = 1, algo_gen = 0,
                                tol_e = tol_e, max_loop = max_loop,
                                rng4, o),
      iterations = 10, check = FALSE
    )
  }
)
print(results_cmd[, c("nerrs", "median", "mem_alloc")])
