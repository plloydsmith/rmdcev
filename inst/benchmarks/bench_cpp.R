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

# Aliases for non-exported internal C++ functions
rmdcev_get_rng        <- rmdcev:::rmdcev_get_rng
rmdcev_get_stream     <- rmdcev:::rmdcev_get_stream
DrawError_rng         <- rmdcev:::DrawError_rng
MarshallianDemand     <- rmdcev:::MarshallianDemand
HicksianDemand        <- rmdcev:::HicksianDemand
ComputeUtilJ          <- rmdcev:::ComputeUtilJ
CalcWTP_rng           <- rmdcev:::CalcWTP_rng
CalcmdemandOne_rng    <- rmdcev:::CalcmdemandOne_rng

income  <- df_sim[["df_indiv"]][["income"]][[1]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[1]]
price   <- df_sim[["df_indiv"]][["price"]][[1]]

quant_num <- as.numeric(income - quant_j %*% price[-1])

# Compute psi_j for one individual (handles ASC-only, attribute-only, or mixed)
psi_ascs <- result[["stan_data"]][["psi_ascs"]]
NPsi_ij  <- result[["stan_data"]][["NPsi_ij"]]
par_psi  <- result[["stan_fit"]][["par"]][["psi"]]

if (psi_ascs == 1 && NPsi_ij == 0) {
    # ASC-only: par_psi has J-1 values; reference alternative gets 0
    psi_j <- c(0, par_psi)
} else if (NPsi_ij > 0) {
    # Attribute-based psi (may include leading J-1 ASC values when psi_ascs == 1)
    attr_psi <- if (psi_ascs == 1) par_psi[nalts:length(par_psi)] else par_psi
    asc_part <- if (psi_ascs == 1) c(0, par_psi[1:(nalts - 1)]) else rep(0, nalts)
    psi_j <- asc_part + result[["stan_data"]][["dat_psi"]][1:nalts, , drop = FALSE] %*% attr_psi
} else {
    psi_j <- rep(0, nalts)
}

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

# Draw errors once for a fixed input to Marshallian/Hicksian benchmarks
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

# ---- CalcmdemandOne_rng: varying nerrs ----
# Benchmarks the full demand pipeline for a single individual
cat("\n--- CalcmdemandOne_rng: nerrs grid ---\n")
results_cmd <- bench::press(
  nerrs = c(1, 10, 50),
  {
    rng3 <- rmdcev_get_rng(seed = 2)
    bench::mark(
      CalcmdemandOne_rng(income, price, psi_j, phi_j, gamma_j, alpha, scale,
                         nerrs = nerrs, model_num = model_num,
                         algo_gen = 0, tol = tol_e, max_loop = max_loop,
                         rng3, o),
      iterations = 10, check = FALSE
    )
  }
)
print(results_cmd[, c("nerrs", "median", "mem_alloc")])

# ---- CalcWTP_rng: varying nerrs using simulation-prepared data ----
cat("\n--- CalcWTP_rng: nerrs grid ---\n")
indiv1_psi_sims   <- df_sim[["df_indiv"]][["psi_sims"]][[1]]
indiv1_phi_sims   <- df_sim[["df_indiv"]][["phi_sims"]][[1]]
indiv1_psi_p_sims <- df_sim[["df_indiv"]][["psi_p_sims"]][[1]]
indiv1_phi_p_sims <- df_sim[["df_indiv"]][["phi_p_sims"]][[1]]

results_wtp <- bench::press(
  nerrs = c(1, 10, 50),
  {
    rng4 <- rmdcev_get_rng(seed = 3)
    bench::mark(
      CalcWTP_rng(
        income            = income,
        quant_j           = quant_j,
        price             = price,
        price_p_policy    = df_sim[["df_common"]][["price_p_list"]],
        psi_p_sims        = indiv1_psi_p_sims,
        phi_p_sims        = indiv1_phi_p_sims,
        psi_sims          = indiv1_psi_sims,
        phi_sims          = indiv1_phi_sims,
        gamma_sims        = df_sim[["df_common"]][["gamma_sim_nonrandom"]],
        alpha_sims        = df_sim[["df_common"]][["alpha_sim_nonrandom"]],
        scale_sims        = df_sim[["df_common"]][["scale_sims"]],
        nerrs             = nerrs,
        cond_error        = 1L,
        draw_mlhs         = 1L,
        algo_gen          = 0L,
        model_num         = model_num,
        price_change_only = 1L,
        tol               = tol_e,
        max_loop          = max_loop,
        rng4, o
      ),
      iterations = 10, check = FALSE
    )
  }
)
print(results_wtp[, c("nerrs", "median", "mem_alloc")])
