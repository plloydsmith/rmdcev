# Parameter recovery test for the uncorrelated random parameters model.
#
# Generates data from a gamma model where psi ASCs vary across individuals
# according to a known normal distribution, estimates the uncorrelated RP model,
# and checks that the population means (mu) and standard deviations (tau) are
# recovered within tolerance.
#
# Also verifies that welfare simulation produces non-empty results, which
# requires the scale_sims fix in StanSimulate.R to be active.

test_that("Uncorrelated RP gamma model recovers population mu and tau", {
    skip_on_cran()

    # -----------------------------------------------------------------------
    # 1. True parameter values
    # -----------------------------------------------------------------------
    n_indiv       <- 500
    n_alts        <- 3
    mu_psi_true   <- c(-2.0, 0.5)   # psi ASC means for alts 2 and 3
    tau_psi_true  <- c(0.6, 0.6)    # psi ASC SDs
    gamma_true    <- 3.0
    alpha_true    <- 0.5

    # -----------------------------------------------------------------------
    # 2. Generate synthetic data
    # -----------------------------------------------------------------------
    set.seed(123)
    income    <- rep(1000, n_indiv)
    price_mat <- matrix(1.0, n_indiv, n_alts)

    # Individual-level psi ASCs: alt 1 = 0 (reference), alts 2-3 are random
    psi_ind       <- matrix(0.0, n_indiv, n_alts)
    for (j in seq(2, n_alts)) {
        psi_ind[, j] <- rnorm(n_indiv, mu_psi_true[j - 1], tau_psi_true[j - 1])
    }

    alpha_parms <- c(alpha_true, rep(0, n_alts))

    PRNG <- rmdcev_get_rng(seed = 42L)
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
        scale     = 1.0,
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
    # 3. Estimate uncorrelated RP gamma model
    # -----------------------------------------------------------------------
    output_rp <- mdcev(
        formula            = ~ 0,
        data               = data_rp,
        model              = "gamma",
        algorithm          = "Bayes",
        random_parameters  = "uncorr",
        psi_ascs           = 1,
        gamma_ascs         = 0,
        gamma_nonrandom    = 1,
        alpha_nonrandom    = 1,
        fixed_scale1       = 1,
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

    # -----------------------------------------------------------------------
    # 5. Verify that welfare simulation runs and returns non-empty results.
    #    Before the StanSimulate.R scale_sims fix, df_common$scale_sim was NULL,
    #    causing the C++ loop to execute 0 times (nsims = 0 in Stan).
    #    After the fix, nrow(wtp[[1]]) == nsims > 0.
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
                label = "Simulation returned non-empty results (scale_sims fix)")
})
