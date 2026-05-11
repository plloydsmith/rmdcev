# Parameter recovery test for the correlated random parameters model.
#
# Generates data from a gamma model where psi ASCs for alts 2 and 3 are drawn
# from a bivariate normal distribution with known means, standard deviations,
# and correlation.
#
# WHY scale = 0.1 (small Gumbel noise):
#   In a standard MDCEV, the Gumbel error variance (pi^2/6 â‰ˆ 1.64 for scale=1)
#   is ~4.5x larger than the psi signal variance (tau^2 = 0.36).  Two effects then
#   destroy the correlation signal in quantity space:
#     1. Noise dominates signal (SNR â‰ˆ 0.18): individual quantities reflect random
#        errors more than correlated psi values.
#     2. Budget substitution: when psi_2 and psi_3 are both high (positive rho),
#        both goods compete for the same income, creating a mechanical negative
#        quantity correlation that cancels the positive psi correlation.
#   With scale = 0.1, the noise variance drops to ~0.016, giving SNR â‰ˆ 0.96.
#   The individual betas are then tightly identified from the data, the MCMC
#   funnel problem is suppressed, and rho converges with 300 sampling draws.
#   The scale is NOT fixed (fixed_scale1 = FALSE) so the model estimates it freely;
#   it should converge to â‰ˆ 0.1.

test_that("Correlated RP gamma model recovers population mu, tau, and rho", {
    skip_on_cran()

    # -----------------------------------------------------------------------
    # 1. True parameter values
    # -----------------------------------------------------------------------
    n_indiv       <- 500
    n_alts        <- 3
    mu_psi_true   <- c(-2.0, 0.5)   # psi ASC means for alts 2 and 3
    tau_psi_true  <- c(0.6, 0.6)    # psi ASC SDs
    rho_true      <- 0.7             # correlation between the two psi ASCs
    gamma_true    <- 3.0
    alpha_true    <- 0.5
    scale_true    <- 0.1             # small scale â†’ high SNR â†’ rho identifiable

    # -----------------------------------------------------------------------
    # 2. Generate synthetic data from bivariate normal psi ASCs
    # -----------------------------------------------------------------------
    set.seed(456)
    income    <- rep(1000, n_indiv)
    price_mat <- matrix(1.0, n_indiv, n_alts)

    # Full covariance matrix for the two random ASCs
    Sigma_true <- diag(tau_psi_true) %*%
        matrix(c(1, rho_true, rho_true, 1), 2, 2) %*%
        diag(tau_psi_true)
    L_chol <- t(chol(Sigma_true))   # lower-triangular Cholesky factor

    # Draw correlated psi ASCs: n_indiv x 2 matrix
    z_draws  <- matrix(rnorm(n_indiv * 2), nrow = 2, ncol = n_indiv)
    psi_corr <- t(L_chol %*% z_draws + mu_psi_true)   # n_indiv x 2

    # Individual-level psi ASCs: alt 1 = 0 (reference), alts 2-3 are random
    psi_ind       <- matrix(0.0, n_indiv, n_alts)
    psi_ind[, 2]  <- psi_corr[, 1]
    psi_ind[, 3]  <- psi_corr[, 2]

    alpha_parms <- c(alpha_true, rep(0, n_alts))

    PRNG <- rmdcev_get_rng(seed = 99L)
    o    <- rmdcev_get_stream()

    df_gen <- list(
        income = as.list(income),
        price  = CreateListsRow(cbind(1, price_mat)),
        psi_j  = lapply(seq_len(n_indiv), function(i) psi_ind[i, ]),
        phi_j  = replicate(n_indiv, rep(0, n_alts), simplify = FALSE)
    )

    quant_list <- purrr::pmap(df_gen, CalcmdemandOne_rng,
        gamma_j   = rep(gamma_true, n_alts),
        alpha     = alpha_parms,
        scale     = scale_true,
        nerrs     = 1,
        model_num = 1,
        algo_gen  = 1,
        tol       = 1e-6,
        max_loop  = 999,
        PRNG, o)

    quant_mat <- matrix(unlist(quant_list), nrow = n_indiv, byrow = TRUE)

    data_rp <- data.frame(
        id     = rep(seq_len(n_indiv), each = n_alts),
        alt    = rep(seq_len(n_alts),  times = n_indiv),
        quant  = as.vector(t(quant_mat[, -1])),
        price  = as.vector(t(price_mat)),
        income = rep(income, each = n_alts)
    )
    data_rp <- mdcev.data(data_rp, id.var = "id", alt.var = "alt", choice = "quant")

    # -----------------------------------------------------------------------
    # 3. Estimate correlated RP gamma model
    #    fixed_scale1 = FALSE: scale is estimated (will converge to â‰ˆ scale_true)
    # -----------------------------------------------------------------------
    output_rp <- mdcev(
        formula            = ~ 0,
        data               = data_rp,
        model              = "gamma",
        algorithm          = "Bayes",
        random_parameters  = "corr",
        psi_ascs           = 1,
        gamma_ascs         = 0,
        gamma_nonrandom    = 1,
        alpha_nonrandom    = 1,
        fixed_scale1       = 0,
        n_chains           = 2,
        n_iterations       = 600,
        print_iterations   = FALSE,
        show_stan_warnings = FALSE,
        backend            = "rstan"
    )

    # -----------------------------------------------------------------------
    # 4. Check parameter recovery
    # -----------------------------------------------------------------------
    draws   <- extract_bayes_draws(output_rp)
    mu_est  <- colMeans(draws[, grep("^mu\\.",  names(draws)), drop = FALSE])
    tau_est <- colMeans(draws[, grep("^tau\\.", names(draws)), drop = FALSE])

    # Population means should be close to the true values
    expect_equal(as.numeric(mu_est), mu_psi_true, tolerance = 0.4,
                 label = "RP posterior mean of mu matches true psi ASC means")

    # Population SDs should be positive and close to the true values
    expect_true(all(tau_est > 0),
                label = "All RP tau estimates are positive")
    expect_equal(as.numeric(tau_est), tau_psi_true, tolerance = 0.4,
                 label = "RP posterior mean of tau matches true psi ASC SDs")

    # L_Omega[2,1] equals rho for a 2-parameter lower-triangular Cholesky factor.
    # The low-noise data design ensures rho is identifiable with 300 sampling draws.
    l_omega_21 <- mean(draws[, "L_Omega.2.1"])
    expect_equal(l_omega_21, rho_true, tolerance = 0.3,
                 label = "L_Omega[2,1] recovers the true ASC correlation")

    # -----------------------------------------------------------------------
    # 5. Verify that welfare simulation runs and returns non-empty results
    # -----------------------------------------------------------------------
    policies <- CreateBlankPolicies(2, output_rp, price_change_only = TRUE)
    df_sim   <- PrepareSimulationData(output_rp, policies, nsims = 10)

    wtp <- mdcev.sim(
        df_sim$df_indiv,
        df_common    = df_sim$df_common,
        sim_options  = df_sim$sim_options,
        cond_err     = 1,
        nerrs        = 1,
        sim_type     = "welfare",
        suppressTime = TRUE
    )

    # Zero price change â†’ zero WTP
    expect_equal(sum(abs(summary(wtp)$CoefTable$mean)), 0, tolerance = 0.01,
                 label = "WTP is zero for zero price change")

    # Simulation must have actually run (nsims rows, not 0)
    expect_true(nrow(wtp[[1]]) > 0,
                label = "Simulation returned non-empty results")
})
