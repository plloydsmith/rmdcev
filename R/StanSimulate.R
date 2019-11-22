
#' @title StanWelfare
#' @description Use Stan functions to simulate Welfare
#' @param df_indiv list of income, quant_j, price, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return wtp list
#' @keywords internal
StanWelfare <- function(df_indiv, df_common, sim_options){

#	df_indiv <- df_wtp$df_indiv
#	df_common <- df_wtp$df_common

	message("Simulating welfare...")

	if (sim_options$price_change_only == FALSE){
		wtp <- purrr::pmap(df_indiv, CalcWTP_rng,
				price_p=df_common$price_p_fixed,
				gamma_sim=df_common$gamma_sim_fixed,
				alpha_sim=df_common$alpha_sim_fixed,
				scale_sim=df_common$scale_sim,
				nerrs=sim_options$nerrs,
				cond_error=sim_options$cond_error,
				draw_mlhs=sim_options$draw_mlhs,
				algo_gen=sim_options$algo_gen,
				model_num=sim_options$model_num,
				tol = sim_options$tol,
				max_loop = sim_options$max_loop)
	} else if (sim_options$price_change_only == TRUE){

		if(!is.null(df_common$gamma_sim_fixed) &
		   !is.null(df_common$alpha_sim_fixed)){

		wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
					price_p=df_common$price_p_list,
					gamma_sim=df_common$gamma_sim_fixed,
					alpha_sim=df_common$alpha_sim_fixed,
					scale_sim=df_common$scale_sim,
					nerrs=sim_options$nerrs,
					cond_error=sim_options$cond_error,
					draw_mlhs=sim_options$draw_mlhs,
					algo_gen=sim_options$algo_gen,
					model_num=sim_options$model_num,
					tol = sim_options$tol,
					max_loop = sim_options$max_loop)

		} else if(!is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p=df_common$price_p_list,
							   gamma_sim=df_common$gamma_sim_fixed,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  !is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p=df_common$price_p_list,
							   alpha_sim=df_common$alpha_sim_fixed,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
							   price_p=df_common$price_p_list,
							   scale_sim=df_common$scale_sim,
							   nerrs=sim_options$nerrs,
							   cond_error=sim_options$cond_error,
							   draw_mlhs=sim_options$draw_mlhs,
							   algo_gen=sim_options$algo_gen,
							   model_num=sim_options$model_num,
							   tol = sim_options$tol,
							   max_loop = sim_options$max_loop)
		}
	}
	return(wtp)
}

#' @title StanDemand
#' @description Use Stan functions to simulate Marshallian demand
#' @param df_indiv list of income, quant_j, price, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return demand with nsim lists of npolsXnalts+1 matrices
#' @keywords internal
StanDemand <- function(df_indiv, df_common, sim_options){

	message("Simulating demand...")

	if (sim_options$price_change_only == FALSE){
		demand <- purrr::pmap(df_indiv, CalcMarshallianDemand_rng,
						   price_p=df_common$price_p_list,
						   gamma_sim=df_common$gamma_sim_fixed,
						   alpha_sim=df_common$alpha_sim_fixed,
						   scale_sim=df_common$scale_sim,
						   nerrs=sim_options$nerrs,
						   cond_error=sim_options$cond_error,
						   draw_mlhs=sim_options$draw_mlhs,
						   algo_gen=sim_options$algo_gen,
						   model_num=sim_options$model_num,
						   tol = sim_options$tol,
						   max_loop = sim_options$max_loop)
	} else if (sim_options$price_change_only == TRUE){

		if(!is.null(df_common$gamma_sim_fixed) &
		   !is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p=df_common$price_p_list,
								  gamma_sim=df_common$gamma_sim_fixed,
								  alpha_sim=df_common$alpha_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop)

		} else if(!is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p=df_common$price_p_list,
								  gamma_sim=df_common$gamma_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  !is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p=df_common$price_p_list,
								  alpha_sim=df_common$alpha_sim_fixed,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop)

		} else if(is.null(df_common$gamma_sim_fixed) &
				  is.null(df_common$alpha_sim_fixed)){

			demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
								  price_p=df_common$price_p_list,
								  scale_sim=df_common$scale_sim,
								  nerrs=sim_options$nerrs,
								  cond_error=sim_options$cond_error,
								  draw_mlhs=sim_options$draw_mlhs,
								  algo_gen=sim_options$algo_gen,
								  model_num=sim_options$model_num,
								  tol = sim_options$tol,
								  max_loop = sim_options$max_loop)
		}
	}
	return(demand)
}
