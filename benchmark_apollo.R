suppressPackageStartupMessages({
  library(rmdcev)
  library(apollo)
  library(dplyr)
  library(tidyr)
})

# ── helpers ──────────────────────────────────────────────────────────────────

to_wide <- function(data_long, nalts) {
  data_long <- as.data.frame(data_long)
  all_b <- grep("^b[0-9]", names(data_long), value = TRUE)
  ind_b <- all_b[sapply(all_b, function(v) {
    all(tapply(data_long[[v]], data_long$id, function(x) length(unique(x)) == 1))
  })]
  ind_vars  <- data_long %>% distinct(id, income, across(all_of(ind_b)))
  alt_vars  <- data_long %>%
    select(id, alt, quant, price, all_of(setdiff(all_b, ind_b))) %>%
    pivot_wider(names_from = alt, values_from = c(quant, price, setdiff(all_b, ind_b)))
  db <- left_join(ind_vars, alt_vars, by = "id")
  pj_xj <- Reduce("+", lapply(1:nalts, function(j)
    db[[paste0("price_", j)]] * db[[paste0("quant_", j)]]))
  db$quant_num <- db$income - pj_xj
  db
}

make_apollo_beta <- function(nalts, npsi_j = 1) {
  betas <- c(psi_ind = 0)
  for (k in seq_len(npsi_j)) betas[paste0("psi_j", k)] <- 0
  for (j in seq_len(nalts))  betas[paste0("ln_g", j)]  <- 0
  betas["ln_sigma"] <- 0
  betas
}

make_apollo_prob_concrete <- function(nalts, ind_col = "b2") {
  all_str      <- paste0('c("outside", ', paste0('"alt', seq_len(nalts), '"', collapse = ", "), ")")
  V_items      <- paste(paste0("alt", seq_len(nalts), " = psi_ind * ", ind_col,
                               " + psi_j1 * b1_", seq_len(nalts)), collapse = ",\n    ")
  alpha_items  <- paste(c("outside = 0", paste0("alt", seq_len(nalts), " = 0")), collapse = ", ")
  gamma_items  <- paste(c("outside = 1",
                           paste0("alt", seq_len(nalts), " = exp(ln_g", seq_len(nalts), ")")),
                        collapse = ", ")
  choice_items <- paste(c("outside = quant_num",
                           paste0("alt", seq_len(nalts), " = quant_", seq_len(nalts))),
                        collapse = ", ")
  cost_items   <- paste(c("outside = 1",
                           paste0("alt", seq_len(nalts), " = price_", seq_len(nalts))),
                        collapse = ", ")
  avail_items  <- paste(c("outside = 1", paste0("alt", seq_len(nalts), " = 1")), collapse = ", ")
  eval(parse(text = sprintf(
'function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  V <- list(outside = 0, %s)
  mdcev_settings <- list(
    alternatives = %s, budget = income, V = V,
    alpha = list(%s), gamma = list(%s), sigma = exp(ln_sigma),
    continuousChoice = list(%s), cost = list(%s), avail = list(%s)
  )
  P[["model"]] <- apollo_mdcev(mdcev_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}', V_items, all_str, alpha_items, gamma_items, choice_items, cost_items, avail_items)))
}

# ── Setup: n=1000, J=50 ──────────────────────────────────────────────────────
n <- 1000
J <- 50
set.seed(42)

cat(sprintf("Generating data (n=%d, J=%d)...\n", n, J))
sim <- GenerateMDCEVData(
  model = "hybrid0", nobs = n, nalts = J,
  psi_j_parms = 0.5, psi_i_parms = -1.5,
  gamma_parms = rep(5, J), scale_parms = 1.0
)

cat("Fitting rmdcev...\n")
fit_r <- mdcev(
  ~ b1 + b2, data = sim$data, psi_ascs = 0, model = "hybrid0",
  algorithm = "MLE", print_iterations = FALSE, backend = "rstan"
)

cat("Fitting Apollo...\n")
database <- to_wide(sim$data, J)
ac  <- list(modelName = "apl_sim", modelDescr = "", indivID = "id",
            outputDirectory = tempdir(), nCores = 1)
ab  <- make_apollo_beta(J)
pfn <- make_apollo_prob_concrete(J)
suppressMessages(suppressWarnings({
  ai <- apollo_validateInputs(apollo_beta = ab, apollo_fixed = character(0),
                              database = database, apollo_control = ac)
  fit_a <- apollo_estimate(ab, character(0), pfn, ai,
             estimate_settings = list(writeIter = FALSE, hessianRoutine = "none",
                                      printLevel = 0L))
}))

# 1% price-change policy
bp  <- matrix(as.numeric(fit_r$stan_data$price_j), nrow = n)
pc  <- 0.01 * colMeans(bp)
pol <- CreateBlankPolicies(npols = 2, fit_r, price_change_only = TRUE)
pol$price_p[[2]] <- c(0, pc)
df_s <- PrepareSimulationData(fit_r, pol, nsims = 1)

db_pol <- database
for (j in seq_len(J))
  db_pol[[paste0("price_", j)]] <- database[[paste0("price_", j)]] + pc[j]
pj_xj_pol <- Reduce("+", lapply(seq_len(J), function(j)
  db_pol[[paste0("price_", j)]] * db_pol[[paste0("quant_", j)]]))
db_pol$quant_num <- db_pol$income - pj_xj_pol
suppressMessages(suppressWarnings({
  ai_pol <- apollo_validateInputs(apollo_beta = fit_a$estimate, apollo_fixed = character(0),
                                  database = db_pol, apollo_control = ac)
}))

# ── Timing: 3 replicates each ────────────────────────────────────────────────
cat("\nTiming rmdcev (3 reps)...\n")
rmdcev_times <- sapply(1:3, function(r) {
  t <- system.time(
    mdcev.sim(df_s$df_indiv, df_common = df_s$df_common,
              sim_options = df_s$sim_options,
              cond_err = 0, nerrs = 100, sim_type = "demand")
  )["elapsed"]
  cat(sprintf("  rep %d: %.2f s\n", r, t)); t
})

cat("\nTiming Apollo (3 reps, 2 predictions each)...\n")
apollo_times <- sapply(1:3, function(r) {
  suppressMessages(suppressWarnings({
    t1 <- system.time(apollo_prediction(fit_a, pfn, ai,
            prediction_settings = list(nRep = 100, silent = TRUE)))["elapsed"]
    t2 <- system.time(apollo_prediction(fit_a, pfn, ai_pol,
            prediction_settings = list(nRep = 100, silent = TRUE)))["elapsed"]
  }))
  t <- t1 + t2
  cat(sprintf("  rep %d: %.2f s (base=%.2f, policy=%.2f)\n", r, t, t1, t2)); t
})

# ── Results ───────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Demand simulation timing: n=%d, J=%d, nerrs=100 ===\n", n, J))
cat(sprintf("rmdcev : mean=%.2f s  sd=%.2f s\n", mean(rmdcev_times), sd(rmdcev_times)))
cat(sprintf("Apollo : mean=%.2f s  sd=%.2f s\n", mean(apollo_times),  sd(apollo_times)))
cat(sprintf("Speedup: %.1fx\n", mean(apollo_times) / mean(rmdcev_times)))
