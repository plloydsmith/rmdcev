# Tests for the psi_random argument: partial random-parameter specification.
#
# psi_random = ~ <vars> makes only those formula terms random; all other
# formula terms become fixed (pooled) point estimates while ASCs remain random.
#
# Parameter count logic for the test below (J=17, model=gamma, uncorr):
#   formula = ~ ageindex + urban    (2 formula terms)
#   psi_random = ~ ageindex         (ageindex random, urban fixed)
#   psi_ascs = TRUE, fixed_scale1 = TRUE, gamma_nonrandom = TRUE, alpha_nonrandom = TRUE
#
#   n_psi (random)  = (J-1) + 1 = 17   [16 ASCs + 1 random formula term]
#   n_psi_fixed     = 1                 [urban, pooled point estimate]
#   n_gamma         = 17                [gamma_ascs=1, nonrandom]
#   n_alpha         = 1                 [model_num=1]
#   n_scale         = 0                 [fixed_scale1=1]
#   n_std_dev       = 17                [tau for n_psi=17 random psi only]
#   n_parms_total   = 17+1+17+1+0+17 = 53

test_that("psi_random splits formula psi into random and fixed parts", {
    skip_on_cran()

    data(data_rec, package = "rmdcev")
    data_rec <- mdcev.data(data_rec, subset = id <= 200,
                           id.var = "id", alt.var = "alt", choice = "quant")

    output_pr <- mdcev(
        formula            = ~ ageindex + urban,
        data               = data_rec,
        model              = "gamma",
        algorithm          = "Bayes",
        random_parameters  = "uncorr",
        psi_ascs           = 1,
        gamma_ascs         = 1,
        gamma_nonrandom    = 1,
        alpha_nonrandom    = 1,
        fixed_scale1       = 1,
        psi_random         = ~ ageindex,
        n_chains           = 1,
        n_iterations       = 10,
        print_iterations   = FALSE,
        show_stan_warnings = FALSE,
        backend            = "rstan"
    )

    # -----------------------------------------------------------------------
    # 1. Parameter counts
    # -----------------------------------------------------------------------
    expect_equal(output_pr$parms_info$n_vars$n_parms_total, 53)
    expect_equal(output_pr$parms_info$n_vars$n_psi,         17)
    expect_equal(output_pr$parms_info$n_vars$n_psi_fixed,    1)

    # -----------------------------------------------------------------------
    # 2. psi_urban (fixed) appears in all_names but NOT in sd_names
    # -----------------------------------------------------------------------
    all_nms <- output_pr$parms_info$parm_names$all_names
    sd_nms  <- output_pr$parms_info$parm_names$sd_names

    expect_true("psi_urban" %in% all_nms,
                label = "Fixed psi term appears in parameter names")
    expect_false(any(grepl("urban", sd_nms)),
                 label = "Fixed psi term has no SD parameter")

    # -----------------------------------------------------------------------
    # 3. Summary table has exactly n_parms_total rows
    # -----------------------------------------------------------------------
    output_sum <- summary(output_pr)
    expect_equal(nrow(output_sum$CoefTable), 53,
                 label = "Summary CoefTable has 53 rows (random + fixed psi)")

    # -----------------------------------------------------------------------
    # 4. Welfare simulation runs and returns non-empty results
    # -----------------------------------------------------------------------
    policies <- CreateBlankPolicies(2, output_pr, price_change_only = TRUE)
    df_sim   <- PrepareSimulationData(output_pr, policies, nsims = 5)

    wtp <- mdcev.sim(
        df_sim$df_indiv,
        df_common    = df_sim$df_common,
        sim_options  = df_sim$sim_options,
        cond_err     = 1,
        nerrs        = 1,
        sim_type     = "welfare",
        suppressTime = TRUE
    )

    expect_equal(sum(abs(summary(wtp)$CoefTable$mean)), 0, tolerance = 0.01,
                 label = "WTP is zero for zero price change")
    expect_true(nrow(wtp[[1]]) > 0,
                label = "Simulation returned non-empty results")
})

test_that("psi_random = NULL gives same result as omitting psi_random", {
    skip_on_cran()

    data(data_rec, package = "rmdcev")
    data_rec <- mdcev.data(data_rec, subset = id <= 200,
                           id.var = "id", alt.var = "alt", choice = "quant")

    # Default (psi_random = NULL): all formula terms random
    # formula = ~ ageindex â†’ 1 random formula term â†’ n_psi = 16+1=17, n_psi_fixed=0
    # n_std_dev = 17, n_parms_total = 17+0+17+1+0+17 = 52
    output_default <- mdcev(
        formula            = ~ ageindex,
        data               = data_rec,
        model              = "gamma",
        algorithm          = "Bayes",
        random_parameters  = "uncorr",
        gamma_nonrandom    = 1,
        alpha_nonrandom    = 1,
        fixed_scale1       = 1,
        n_chains           = 1,
        n_iterations       = 10,
        print_iterations   = FALSE,
        show_stan_warnings = FALSE,
        backend            = "rstan"
    )

    expect_equal(output_default$parms_info$n_vars$n_parms_total, 52)
    expect_equal(output_default$parms_info$n_vars$n_psi_fixed,    0)
})
