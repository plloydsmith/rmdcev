#' @title maxlikeMDCEV
#' @description Fit a MDCEV model with MLE
#' @param stan_data data for model formatted from processMDCEVdata
#' @inheritParams mdcev
#' @param mle_options modeling options for MLE
#' @param parms_info information on parameters
#' @param backend character, either "cmdstanr" (default) or "rstan"
#' @noRd
maxlikeMDCEV <- function(stan_data,
                         mle_options,
                         parms_info,
                         backend = "cmdstanr", ...) {

    result <- list()

    if (!is.list(mle_options$initial.parameters) || mle_options$n_classes == 1) {

        message("Using MLE to estimate KT model")

        # Ensure single class used for base model
        stan_data_temp              <- stan_data
        stan_data_temp$K            <- 1
        stan_data_temp$L            <- 0
        stan_data_temp$data_class   <- matrix(0, stan_data$I, 0)

        mle_options$init <- CleanInit(mle_options$initial.parameters)

        # Set starting values for scale to 1
        if (stan_data$fixed_scale1 == 0 && is.null(mle_options$init[["scale"]]))
            mle_options$init[["scale"]] <- array(1, dim = 1)

        if (backend == "cmdstanr") {
            stan_fit <- RunCmdStanOptimizing(
                stan_data_temp, "mdcev", mle_options, ...
            )
        } else {
            stan_fit <- rstan::optimizing(
                get_rstan_model("mdcev"),
                data        = stan_data_temp,
                as_vector   = FALSE,
                seed        = mle_options$seed,
                verbose     = mle_options$print_iterations,
                iter        = mle_options$max_iterations,
                init        = mle_options$init,
                draws       = mle_options$n_draws,
                hessian     = mle_options$hessian,
                ...
            )
        }

        if (mle_options$keep_loglik == 0)
            stan_fit <- ReduceStanFitSize(stan_fit, parms_info)

        result$stan_fit            <- stan_fit
        result$stan_fit$par[["theta"]] <- NULL
        result$stan_fit$par[["delta"]] <- NULL
        result$log.likelihood      <- stan_fit[["par"]][["sum_log_lik"]]
        result$effective.sample.size <- sum(stan_data$weights)
    }

    if (mle_options$n_classes > 1) {
        if (!is.list(mle_options$initial.parameters)) {
            result$mdcev_fit           <- result$stan_fit
            result$mdcev_log.likelihood <- result$log.likelihood

            # Build LC initial values from single-class estimates
            init.par  <- stan_fit$par
            init.psi  <- init.par$psi

            if (length(init.psi) > 0) {
                init.shift <- seq(-0.02, 0.02, length.out = length(init.psi))
                for (i in seq_along(init.psi))
                    init.psi[i] <- init.psi[i] + init.shift[i]
            }
            init.psi <- matrix(init.psi,
                               nrow  = stan_data$K,
                               ncol  = length(init.psi),
                               byrow = TRUE)

            init <- list(psi = init.psi)

            if (stan_data$fixed_scale1 == 0 && stan_data$single_scale == 0)
                init$scale <- rep(stan_fit$par[["scale"]], stan_data$K)
            else if (stan_data$fixed_scale1 == 0 && stan_data$single_scale == 1)
                init$scale <- stan_fit$par[["scale"]]

            if (stan_data$model_num %in% c(1, 3, 5))
                init$alpha <- matrix(rep(init.par$alpha, stan_data$K),
                                     nrow = stan_data$K, ncol = 1)
            else if (stan_data$model_num == 2)
                init$alpha <- matrix(rep(init.par$alpha, stan_data$K),
                                     nrow = stan_data$K, ncol = stan_data$J + 1)

            if (stan_data$model_num != 2) {
                init.gamma <- init.par$gamma
                init.shift <- seq(-0.02, 0.02, length.out = length(init.gamma))
                for (i in seq_along(init.gamma))
                    init.gamma[i] <- init.gamma[i] + init.shift[i]

                if (stan_data$gamma_ascs == 1)
                    init$gamma <- matrix(rep(init.gamma, stan_data$K),
                                         nrow = stan_data$K, ncol = stan_data$J)
                else if (stan_data$gamma_ascs == 0)
                    init$gamma <- matrix(rep(init.gamma, stan_data$K),
                                         nrow = stan_data$K, ncol = 1)
            }

            if (stan_data$model_num == 5)
                init$phi <- matrix(rep(init.par$phi, stan_data$K),
                                   nrow = stan_data$K, ncol = stan_data$NPhi)

        } else {
            init <- mle_options$initial.parameters
        }

        message("Using MLE to estimate LC-KT")

        if (backend == "cmdstanr") {
            mle_options$init <- init
            stan_fit <- RunCmdStanOptimizing(stan_data, "mdcev", mle_options, ...)
        } else {
            stan_fit <- rstan::optimizing(
                get_rstan_model("mdcev"),
                data      = stan_data,
                as_vector = FALSE,
                seed      = mle_options$seed,
                init      = init,
                verbose   = mle_options$print_iterations,
                iter      = mle_options$max_iterations,
                draws     = mle_options$n_draws,
                hessian   = mle_options$hessian,
                ...
            )
        }

        if (mle_options$keep_loglik == 0)
            stan_fit <- ReduceStanFitSize(stan_fit, parms_info)

        result$stan_fit         <- stan_fit
        result$log.likelihood   <- stan_fit[["par"]][["sum_log_lik"]]
        class_probabilities     <- stan_fit[["par"]][["theta"]]
        colnames(class_probabilities) <- paste0("class", seq_len(mle_options$n_classes))
        result$class_probabilities <- class_probabilities
    }

    result$backend <- backend
    return(result)
}

#' @title RunCmdStanOptimizing
#' @description Run CmdStan optimization (MLE) for MDCEV models.
#' @param stan_data List of data passed to Stan.
#' @param model_name Name of the Stan model file (without extension).
#' @param mle_options List of MLE options.
#' @param ... Additional arguments passed to cmdstanr \code{$optimize()}.
#' @return A list with \code{par} (named parameter list), \code{hessian},
#'   \code{theta_tilde}, and \code{return_code} mimicking rstan's optimizing output.
#' @noRd
RunCmdStanOptimizing <- function(stan_data, model_name, mle_options, ...) {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
        stop("Package 'cmdstanr' is required for backend = 'cmdstanr'.")

    stan_file <- system.file("stan", paste0(model_name, ".stan"), package = "rmdcev")
    if (!nzchar(stan_file))
        stop("Stan model file not found: ", model_name, ".stan")

    mod <- cmdstanr::cmdstan_model(stan_file)

    init_val <- if (!is.null(mle_options$init)) mle_options$init else "random"

    fit <- mod$optimize(
        data         = stan_data,
        init         = init_val,
        seed         = as.integer(mle_options$seed),
        iter         = mle_options$max_iterations,
        jacobian     = FALSE,
        ...
    )

    # Convert CmdStan output to rstan-like list for downstream compatibility
    draws <- fit$draws(format = "draws_df")
    par_names <- setdiff(colnames(draws), c(".chain", ".iteration", ".draw", "lp__"))

    par_vec <- as.numeric(draws[1, par_names])
    names(par_vec) <- par_names

    # Reconstruct named list with correct array dimensions (matching rstan::optimizing output)
    par_list <- cmdstan_pars_to_list(par_vec)

    # Hessian: cmdstanr optimize does not expose the Hessian natively;
    # use rstan fallback if hessian = TRUE, otherwise set to NULL.
    hessian_mat <- NULL
    if (isTRUE(mle_options$hessian)) {
        if (requireNamespace("rstan", quietly = TRUE)) {
            rstan_fit <- rstan::optimizing(
                get_rstan_model(model_name),
                data      = stan_data,
                as_vector = FALSE,
                seed      = mle_options$seed,
                verbose   = FALSE,
                iter      = mle_options$max_iterations,
                init      = if (!is.null(mle_options$init)) mle_options$init else "random",
                draws     = mle_options$n_draws,
                hessian   = TRUE
            )
            hessian_mat    <- rstan_fit$hessian
            theta_tilde    <- rstan_fit$theta_tilde
            return_code    <- rstan_fit$return_code
        } else {
            warning("Hessian computation requires rstan. Install rstan or set hessian = FALSE. Returning NULL hessian.")
            theta_tilde <- matrix(par_vec, nrow = 1)
            return_code <- fit$return_codes()
        }
    } else {
        theta_tilde <- matrix(par_vec, nrow = 1)
        return_code <- fit$return_codes()
    }

    list(
        par         = par_list,
        hessian     = hessian_mat,
        theta_tilde = theta_tilde,
        return_code = if (length(return_code) > 0) return_code[1] else 0L
    )
}

#' @title cmdstan_pars_to_list
#' @description Convert a flat named parameter vector from CmdStan to a named
#'   list with properly-dimensioned arrays, matching the output of
#'   \code{rstan::optimizing(as_vector = FALSE)}.
#' @param par_vec Named numeric vector of parameters.
#' @return A named list.
#' @noRd
cmdstan_pars_to_list <- function(par_vec) {
    # Parameter names from CmdStan look like "psi[1]", "psi[2]", "gamma[1,1]", etc.
    # Group them into arrays.
    par_list  <- list()
    base_names <- sub("\\[.*\\]$", "", names(par_vec))
    unique_bases <- unique(base_names)

    for (bn in unique_bases) {
        idx <- which(base_names == bn)
        vals <- par_vec[idx]
        if (length(vals) == 1 && !grepl("\\[", names(vals))) {
            par_list[[bn]] <- unname(vals)
        } else {
            # Parse indices to rebuild correct array dimensions
            raw_idx <- sub(paste0("^", bn, "\\[(.*)\\]$"), "\\1", names(vals))
            dim_mat <- do.call(rbind, lapply(strsplit(raw_idx, ","), as.integer))
            if (ncol(dim_mat) == 1) {
                arr <- array(unname(vals), dim = max(dim_mat))
            } else {
                dims <- apply(dim_mat, 2, max)
                arr  <- array(unname(vals), dim = dims)
            }
            par_list[[bn]] <- arr
        }
    }
    par_list
}

#' @title ReduceStanFitSize
#' @description Reduces the size of the stan fit object by removing log-likelihood
#'   draws and trimming \code{theta_tilde} to model parameters only.
#' @param stan_fit A stanfit-like list (from rstan or cmdstanr wrapper).
#' @param parms_info Information on parameters.
#' @return A reduced-size fit list.
#' @noRd
ReduceStanFitSize <- function(stan_fit, parms_info) {
    stan_fit[["par"]][["log_like"]] <- NULL
    stan_fit[["theta_tilde"]] <- stan_fit[["theta_tilde"]][
        , seq_len(parms_info$n_vars$n_parms_total), drop = FALSE
    ]
    return(stan_fit)
}
