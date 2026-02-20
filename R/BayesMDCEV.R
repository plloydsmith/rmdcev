#' @title BayesMDCEV
#' @description Fit a MDCEV model using Bayesian estimation and Stan
#' @param bayes_options list of Bayes options
#' @param stan_data data for model
#' @inheritParams mdcev
#' @param keep.samples default is FALSE
#' @param include.stanfit default is TRUE
#' @param backend character, either "cmdstanr" (default) or "rstan"
#' @import dplyr
#' @noRd
BayesMDCEV <- function(stan_data, bayes_options,
                       keep.samples = FALSE,
                       include.stanfit = TRUE,
                       backend = "cmdstanr", ...) {

    if (bayes_options$n_iterations <= 0)
        stop("The specified number of iterations must be greater than 0.")

    # Create indices for individual level psi parameters
    indexes <- tibble::tibble(
        individual = rep(1:stan_data$I, each = stan_data$J),
        task = rep(1:stan_data$I, each = stan_data$J),
        row = 1:(stan_data$I * stan_data$J)
    ) %>%
        dplyr::group_by(.data$task) %>%
        dplyr::summarise(
            task_individual = dplyr::first(.data$individual),
            start = dplyr::first(.data$row),
            end = dplyr::last(.data$row)
        )

    stan_data$start         <- indexes$start
    stan_data$end           <- indexes$end
    stan_data$task_individual <- indexes$task_individual
    stan_data$task          <- indexes$task
    stan_data$IJ            <- stan_data$I * stan_data$J
    stan_data$lkj_shape     <- bayes_options$lkj_shape_prior

    stan_data$K            <- 1
    stan_data$L            <- 0
    stan_data$data_class   <- matrix(0, stan_data$I, 0)

    if (bayes_options$random_parameters == "fixed") {
        message("Using Bayes to estimate KT model")
        model_name <- "mdcev"
    } else if (bayes_options$random_parameters == "uncorr") {
        message("Using Bayes to estimate uncorrelated RP-KT")
        model_name <- "mdcev_rp"
        stan_data$corr <- 0
    } else if (bayes_options$random_parameters == "corr") {
        message("Using Bayes to estimate correlated RP-KT")
        model_name <- "mdcev_rp"
        stan_data$corr <- 1
    }

    if (backend == "cmdstanr") {
        if (!requireNamespace("cmdstanr", quietly = TRUE))
            stop("Package 'cmdstanr' is required for backend = 'cmdstanr'. ",
                 "Install it with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))")

        stan_fit <- RunCmdStanSampling(stan_data, model_name, bayes_options, ...)

    } else if (backend == "rstan") {
        if (!requireNamespace("rstan", quietly = TRUE))
            stop("Package 'rstan' is required for backend = 'rstan'. ",
                 "Install it with: install.packages('rstan')")

        stan.model <- get_rstan_model(model_name)

        if (bayes_options$show_stan_warnings == FALSE)
            suppressWarnings(stan_fit <- RunStanSampling(stan_data, stan.model, bayes_options, ...))
        else
            stan_fit <- RunStanSampling(stan_data, stan.model, bayes_options, ...)
    } else {
        stop("backend must be 'cmdstanr' or 'rstan'")
    }

    if (bayes_options$n_chains == 1)
        chain_index <- 1
    else
        chain_index <- bayes_options$n_chains + 1

    result <- list()
    result$stan_fit <- stan_fit
    result$backend  <- backend

    if (backend == "cmdstanr") {
        draws_df <- stan_fit$draws(variables = "sum_log_lik", format = "draws_matrix")
        result$log.likelihood <- mean(as.numeric(draws_df))
    } else {
        result$log.likelihood <- rstan::get_posterior_mean(
            result$stan_fit, pars = "sum_log_lik"
        )[, chain_index]
    }

    return(result)
}

#' @title RunCmdStanSampling
#' @description Wrapper function for cmdstanr sampling.
#' @inheritParams BayesMDCEV
#' @param model_name Name of Stan model ("mdcev" or "mdcev_rp")
#' @return A CmdStanMCMC object.
#' @noRd
RunCmdStanSampling <- function(stan_data, model_name, bayes_options, ...) {
    stan_file <- system.file("stan", paste0(model_name, ".stan"), package = "rmdcev")
    if (!nzchar(stan_file))
        stop("Stan model file not found: ", model_name, ".stan")

    # Workaround: cmdstanr wraps include-paths with spaces in single quotes,
    # which Windows does not interpret as quoting. Fall back to installed path.
    if (.Platform$OS.type == "windows" && grepl(" ", stan_file)) {
        lib_paths <- .libPaths()
        installed_pkg <- file.path(lib_paths, "rmdcev")
        installed_pkg <- installed_pkg[file.exists(installed_pkg)][1]
        if (!is.na(installed_pkg)) {
            candidate <- file.path(installed_pkg, "stan",
                                   paste0(model_name, ".stan"))
            if (file.exists(candidate))
                stan_file <- candidate
        }
    }

    mod <- cmdstanr::cmdstan_model(stan_file, include_paths = dirname(stan_file))

    init_fn <- bayes_options$initial.parameters
    # cmdstanr init handling:
    # - "random" / NULL  → use 0.5 (small radius; avoids large random starts
    #   for the RP z-matrix while still giving non-zero gradients for HMC)
    # - flat named list  → wrap as list of n_chains lists for cmdstanr sample
    if (!is.list(init_fn) && (is.null(init_fn) || identical(init_fn, "random"))) {
        init_fn <- 0.5
    } else if (is.list(init_fn) && length(init_fn) > 0 &&
               !is.list(init_fn[[1L]])) {
        init_fn <- rep(list(init_fn), bayes_options$n_chains)
    }

    n_warmup  <- floor(bayes_options$n_iterations / 2)
    n_samples <- bayes_options$n_iterations - n_warmup

    mod$sample(
        data            = stan_data,
        chains          = bayes_options$n_chains,
        parallel_chains = bayes_options$n_cores,
        init            = init_fn,
        iter_warmup     = n_warmup,
        iter_sampling   = n_samples,
        seed            = as.integer(bayes_options$seed),
        max_treedepth   = bayes_options$max_tree_depth,
        adapt_delta     = bayes_options$adapt_delta,
        show_messages   = bayes_options$show_stan_warnings,
        ...
    )
}

#' @title RunStanSampling
#' @description Wrapper function for \code{rstan::sampling} (rstan backend).
#' @inheritParams BayesMDCEV
#' @param stan.model Compiled rstan stanmodel object
#' @param ... Additional parameters to pass on to \code{rstan::sampling}.
#' @return A stanfit object.
#' @noRd
RunStanSampling <- function(stan_data, stan.model, bayes_options, ...) {
    if (stan_data$fixed_scale1 == 0 && !is.list(bayes_options$initial.parameters))
        bayes_options$initial.parameters <- rep(
            list(list(scale = array(1, dim = 1))),
            bayes_options$n_chains
        )

    rstan::sampling(
        stan.model,
        data    = stan_data,
        chains  = bayes_options$n_chains,
        cores   = bayes_options$n_cores,
        init    = bayes_options$initial.parameters,
        iter    = bayes_options$n_iterations,
        seed    = bayes_options$seed,
        control = list(
            max_treedepth = bayes_options$max_tree_depth,
            adapt_delta   = bayes_options$adapt_delta
        ),
        ...
    )
}
