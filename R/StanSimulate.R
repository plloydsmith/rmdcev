
#' @title StanWelfare
#' @description Use Stan functions to simulate Welfare
#' @inheritParams mdcev.sim
#' @return wtp list
#' @keywords internal
StanWelfare <- function(df_indiv, df_common, sim_options, stan_seed){

	PRNG <-rstan::get_rng(seed = stan_seed)
	o <- rstan::get_stream()

#	df_indiv <- df_wtp$df_indiv
#	df_common <- df_wtp$df_common

	message("Simulating welfare...")

	if (sim_options$price_change_only == FALSE){
		wtp <- purrr::pmap(df_indiv, CalcWTP_rng,
				price_p_policy=df_common$price_p_fixed,
				gamma_sim=df_common$gamma_sim_fixed,
				alpha_sim=df_common$alpha_sim_fixed,
				scale_sim=df_common$scale_sim,
				nerrs=sim_options$nerrs,
				cond_error=sim_options$cond_error,
				draw_mlhs=sim_options$draw_mlhs,
				algo_gen=sim_options$algo_gen,
				model_num=sim_options$model_num,
				tol = sim_options$tol,
				max_loop = sim_options$max_loop,
				PRNG, o)
	} else if (sim_options$price_change_only == TRUE){

		if(!is.null(df_common$gamma_sim_fixed) &
		   !is.null(df_common$alpha_sim_fixed)){

		wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
					price_p_policy=df_common$price_p_list,
					gamma_sim=df_common$gamma_sim_fixed,
					alpha_sim=df_common$alpha_sim_fixed,
					scale_sim=df_common$scale_sim,
					nerrs=sim_options$nerrs,
					cond_error=sim_options$cond_error,
					draw_mlhs=sim_options$draw_mlhs,
					algo_gen=sim_options$algo_gen,
					model_num=sim_options$model_num,
					tol = sim_options$tol,
					max_loop = sim_options$max_loop,
					PRNG, o)

		} else if(!is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p_policy=df_common$price_p_list,
							   gamma_sim=df_common$gamma_sim_fixed,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  !is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p_policy=df_common$price_p_list,
							   alpha_sim=df_common$alpha_sim_fixed,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p_policy=df_common$price_p_list,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop,
							   PRNG, o)
		}
	}
	return(wtp)
}

#' @title StanDemand
#' @description Use Stan functions to simulate Marshallian demand
#' @inheritParams mdcev.sim
#' @return demand with nsim lists of npolsXnalts+1 matrices
#' @keywords internal
StanDemand <- function(df_indiv, df_common, sim_options, stan_seed){

	PRNG <-rstan::get_rng(seed = stan_seed)
	o <- rstan::get_stream()

	message("Simulating demand...")

	if (sim_options$price_change_only == FALSE){
		demand <- purrr::pmap(df_indiv, CalcMarshallianDemand_rng,
						   price_p_policy=df_common$price_p_list,
						   gamma_sim=df_common$gamma_sim_fixed,
						   alpha_sim=df_common$alpha_sim_fixed,
						   scale_sim=df_common$scale_sim,
						   nerrs=sim_options$nerrs,
						   cond_error=sim_options$cond_error,
						   draw_mlhs=sim_options$draw_mlhs,
						   algo_gen=sim_options$algo_gen,
						   model_num=sim_options$model_num,
						   tol = sim_options$tol,
						   max_loop = sim_options$max_loop,
							  PRNG, o)
	} else if (sim_options$price_change_only == TRUE){

		if(!is.null(df_common$gamma_sim_fixed) &
		   !is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p_policy=df_common$price_p_list,
								  gamma_sim=df_common$gamma_sim_fixed,
								  alpha_sim=df_common$alpha_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop,
								  PRNG, o)

		} else if(!is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p_policy=df_common$price_p_list,
								  gamma_sim=df_common$gamma_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop,
								  PRNG, o)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  !is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p_policy=df_common$price_p_list,
								  alpha_sim=df_common$alpha_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop,
								  PRNG, o)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p_policy=df_common$price_p_list,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop,
								  PRNG, o)
		}
	}
	return(demand)
}
