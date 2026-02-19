
#' @title StanSimulate
#' @description Use Stan functions to simulate Welfare and Demand
#' @inheritParams mdcev.sim
#' @return out list
#' @noRd
StanSimulate <- function(df_indiv, df_common, sim_options, stan_seed){

	PRNG <- rmdcev_get_rng(seed = stan_seed)
	o <- rmdcev_get_stream()

	sim_options[["price_change_only"]] <- as.integer(sim_options[["price_change_only"]])


	if (sim_options$sim_type == "welfare")
		stan_function <- CalcWTP_rng
	else
		stan_function <- CalcMarshallianDemand_rng


		if(!is.null(df_common$gamma_sim_nonrandom) &
		   !is.null(df_common$alpha_sim_nonrandom)){

		out <- purrr::pmap(df_indiv, stan_function,
					price_p_policy=df_common$price_p_list,
					gamma_sims=df_common$gamma_sim_nonrandom,
					alpha_sims=df_common$alpha_sim_nonrandom,
					scale_sims=df_common$scale_sim,
					nerrs=sim_options$nerrs,
					cond_error=sim_options$cond_error,
					draw_mlhs=sim_options$draw_mlhs,
					algo_gen=sim_options$algo_gen,
					model_num=sim_options$model_num,
					price_change_only = sim_options$price_change_only,
					tol = sim_options$tol,
					max_loop = sim_options$max_loop,
					PRNG, o)

		} else if(!is.null(df_common$gamma_sim_nonrandom) &
				  is.null(df_common$alpha_sim_nonrandom)){

			out <- purrr::pmap(df_indiv, stan_function,
							   price_p_policy=df_common$price_p_list,
							   gamma_sims=df_common$gamma_sim_nonrandom,
							   scale_sims=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   price_change_only = sim_options$price_change_only,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)

		} else if(is.null(df_common$gamma_sim_nonrandom) &
				  !is.null(df_common$alpha_sim_nonrandom)){

			out <- purrr::pmap(df_indiv, stan_function,
							   price_p_policy=df_common$price_p_list,
							   alpha_sims=df_common$alpha_sim_nonrandom,
							   scale_sims=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   price_change_only = sim_options$price_change_only,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)

		} else if(is.null(df_common$gamma_sim_nonrandom) &
				  is.null(df_common$alpha_sim_nonrandom)){

			out <- purrr::pmap(df_indiv, stan_function,
							   price_p_policy=df_common$price_p_list,
							   scale_sims=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
						       price_change_only = sim_options$price_change_only,
						       tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)
		}

	return(out)
}
