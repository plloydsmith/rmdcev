//Code for MDCEV Simulation Functions <- "
functions {
vector[] DrawError2_rng(real quant_num, vector quant_j, vector price_j,
					vector psi_j, vector gamma_j, vector alpha, real scale,
					int ngoods, int nerrs, int cond_error){

	vector[ngoods+1] error[nerrs];

	if (cond_error == 0){	// Unconditional
		for(err in 1:nerrs){
			error[err] = rep_vector(0, ngoods+1);
			for(j in 2:ngoods+1) // keep outside good error to zero
				error[err, j] = -log(-log(uniform_rng(0, 1))) * scale; //uniform(0,1) draws
		}
	} else if (cond_error == 1){	// Conditional
		// Cacluate the demand function, g
		vector[ngoods + 1] cond_demand = append_row(quant_num, quant_j);
		// compute vk and v1
		real v1 = (alpha[1] - 1) * log(quant_num);
		vector[ngoods] vk = psi_j + (alpha[2:ngoods+1] - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
		// ek = v1 - vk and assume error term is zero for outside good
		vector[ngoods + 1] ek = append_row(0, (v1 - vk) / scale);

		// Calculate errors
		// For unvisited alternatives, draw from truncated multivariate logistic distribution
		for(err in 1:nerrs)
			for(j in 1:ngoods+1)
				error[err, j] = cond_demand[j] > 0 ? ek[j] * scale :
							-log(-log(uniform_rng(0, 1) * exp(-exp(-ek[j])))) * scale;
		}
return(error);
}

row_vector Shuffle_rng(row_vector inv, int nerrs){

	row_vector[nerrs] out;
	row_vector[nerrs] temp1 = rep_row_vector(0, nerrs) ;
	row_vector[nerrs] temp2 = to_row_vector(uniform_rng(temp1, 1));
	out = inv[sort_indices_asc(temp2)];

return(out);
}

/**Generate random draws using Modified Latin Hypercube Sampling algorithm or uniform
#' Algorithm described in
#' Hess, S., Train, K., and Polak, J. (2006) Transportation Research 40B, 147 - 163.
**/
row_vector DrawMlhs_rng(int nerrs, int draw_mlhs){

	row_vector[nerrs] error;
	row_vector[nerrs] temp0 = rep_row_vector(0, nerrs);

	if(draw_mlhs == 0){
		error = to_row_vector(uniform_rng(temp0, 1));
	} else if(draw_mlhs == 1){
		int temp1[nerrs];
		row_vector[nerrs] temp;

		for (err in 1:nerrs)
			temp1[err] = err - 1;

		temp = to_row_vector(temp1) / nerrs;
		error = Shuffle_rng(temp + to_row_vector(uniform_rng(temp0, 1))/ nerrs, nerrs);
	}
return(error);
}

vector[] DrawError_rng(real quant_num, vector quant_j, vector price_j,
					vector psi_j, vector gamma_j, vector alpha, real scale,
					int ngoods, int nerrs, int cond_error, int draw_mlhs){

	vector[ngoods+1] out[nerrs];
	matrix[ngoods+1, nerrs] error;
	matrix[nerrs, ngoods+1] error_t;
	vector[nerrs] temp0 = rep_vector(0, nerrs);

	if (cond_error == 0){	// Unconditional
		error[1] = rep_row_vector(0, nerrs);
		for(j in 2:ngoods+1)
			error[j] = -log(-log(DrawMlhs_rng(nerrs, draw_mlhs))) * scale;

	} else if (cond_error == 1){	// Conditional
		// Cacluate the demand function, g
		vector[ngoods + 1] cond_demand = append_row(quant_num, quant_j);
		// compute vk and v1
		real v1 = (alpha[1] - 1) * log(quant_num);
		vector[ngoods] vk = psi_j + (alpha[2:ngoods+1] - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
		// ek = v1 - vk and assume error term is zero for outside good
		vector[ngoods + 1] ek = append_row(0, (v1 - vk) / scale);

		// Calculate errors
		// For unvisited alternatives, draw from truncated multivariate logistic distribution
		for(j in 1:ngoods+1)
			error[j] = cond_demand[j] > 0 ? rep_row_vector(ek[j] * scale, nerrs) :
					-log(-log(DrawMlhs_rng(nerrs, draw_mlhs) * exp(-exp(-ek[j])))) * scale;
		}
	error_t = error';
	for(err in 1:nerrs)
		out[err] = error_t[err]';
return(out);
}

/**
 * Returns an ranking vector of non-numeraire goods by MUzero
 * @return integer vector of ordered goods by MUzero
 */
int[] CalcGoodOrder(vector MUzero, int ngoods) {

	int order_x[ngoods+1];
	vector[ngoods] ord_goods;
	int order_MU[ngoods] = sort_indices_desc(MUzero[2:ngoods+1]);

	for (j in 1:ngoods)
		ord_goods[j] = j;

	order_x = sort_indices_asc(append_row(1.0, to_vector(ord_goods[order_MU] + 1)));
return(order_x);
}
/**
 * Sort matrix of parameters/prices
 * @return Matrix of ngoods x 5
 */
matrix SortParmMatrix(vector MUzero, vector price, vector gamma, vector alpha, int ngoods) {

	matrix[ngoods+1, 4] parm_matrix;
	vector[ngoods] MU_j = MUzero[2:ngoods+1];
	vector[ngoods] price_j = price[2:ngoods+1];
	vector[ngoods] gamma_j = gamma[2:ngoods+1];
	vector[ngoods] alpha_j = alpha[2:ngoods+1];
	int order_MU[ngoods] = sort_indices_desc(MU_j);	// find ranking of non-numeraire goods by MUzero

	parm_matrix = append_col(append_row(MUzero[1], MU_j[order_MU]),
	append_col(append_row(price[1], price_j[order_MU]),
	append_col(append_row(gamma[1], gamma_j[order_MU]), append_row(alpha[1], alpha_j[order_MU]))));
return(parm_matrix);
}

real ComputeE(int M, real lambda, vector g_price, vector b, vector c, vector d){
	real output;
	vector[M] temp;
	for (m in 1:M)
		temp[m] = g_price[m] * (lambda ^ b[m] / c[m] -d[m]);
	output =  sum(temp);
return(output);
}
/**
 * Calculate MarshDemands using general or hybrid approach
 * @return vector of ngood demands
 */
vector MarshallianDemand(real inc, vector price, vector MUzero, vector gamma, vector alpha,
							int ngoods, int algo_gen, real tol_e, int max_loop) {

	vector[ngoods+1] mdemand;
	real lambda;
	int M = 1; // Indicator of which ordered alteratives (<=M) are being considered
	int exit = 0;
	real E;
	int order_x[ngoods+1] = CalcGoodOrder(MUzero, ngoods);
	vector[ngoods+1] X = rep_vector(0, ngoods+1);
	vector[ngoods+1] d = append_row(0, rep_vector(1, ngoods));
	matrix[ngoods+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, ngoods);
	vector[ngoods+1] mu = col(parm_matrix, 1); // obtain mu
	vector[ngoods+1] g = col(parm_matrix, 3); // obtain gamma
	vector[ngoods+1] g_price = g .* col(parm_matrix, 2);

	if (algo_gen == 0) { //Hybrid
		real lambda_num;
		real lambda_den;
		real alpha_1 = alpha[1];
		vector[ngoods+1] b;
		vector[ngoods+1] c;

		for (j in 1:ngoods+1)
			b[j] = pow(mu[j], inv(1 - alpha_1));

		c = g_price .* b;

		while (exit == 0){
			// Calculate lambda equal to MUzero(M+1)
			lambda_num = inc + sum(g_price[1:M]) - 1; // minus one for numeraire
			lambda_den = sum(c[1:M]);
			lambda = pow(lambda_num / lambda_den, alpha_1 - 1);

			//Compare lambda to baseline utility of the next lowest alternative
			//(M+1). If lambda exceeds this value then all lower-valued
			//alternatives have zero demand.

			if (lambda > mu[min(M + 1, ngoods + 1)] || M == ngoods+1){
				// Compute demands (using eq. 12 in Pinjari and Bhat)
				for (m in 1:M)
					X[m] = (pow(lambda / mu[m], 1 / (alpha_1 - 1)) - d[m]) * g[m];
				exit = 1;

			} else if ( M < ngoods + 1)
				M += 1; // adds one to M
			}

	} else if (algo_gen == 1) {//	General
		real lambda_l;
		real lambda_u;
		vector[ngoods+1] c;
		vector[ngoods+1] b_temp = col(parm_matrix, 4);
		vector[ngoods+1] b = inv(append_row(b_temp[1], b_temp[2:ngoods+1]) - 1);

		for (j in 1:ngoods+1)
			c[j] = mu[j] ^ b[j];

		while (exit == 0){
			lambda = mu[M + 1];	// Calculate lambda equal to MUzero(M+1)

			E = ComputeE(M, lambda, g_price, b, c, d);			// Calculate E

			if (E >= inc || M+1 == ngoods+1){
				if(E < inc)
					M += 1;
//				M = E < inc ? M + 1 : M;
				lambda_l = E < inc ? 0 : lambda;
				lambda_u = mu[M];
				lambda = (lambda_l + lambda_u) / 2;

				for (n in 1:max_loop){
					real lambda_mid = (lambda_l + lambda_u) / 2;

					E = ComputeE(M, lambda, g_price, b, c, d);

					// Update lambdas's
					if (E < inc)
						lambda_u = lambda_mid;
					else if (E > inc)
						lambda_l = lambda_mid;
//					lambda_u = E > inc ? lambda_u : (lambda_l + lambda_u) / 2;
//					lambda_l = E < inc ? lambda_l : (lambda_l + lambda_u) / 2;
					lambda = (lambda_l + lambda_u) / 2;

					if (fabs((E - inc) / (E + inc) * 0.5) < tol_e) break;
				}
				// Compute demands (using eq. 12 in Pinjari and Bhat)
				for (m in 1:M)
					X[m] = ((lambda / mu[m]) ^ b[m] -d[m]) * g[m];
				exit = 1;

			} else if (E < inc && M + 1 < ngoods + 1)
				M += 1; // adds one to M
		}
	}
	// This code puts the choices back in their original order and exports demands
	mdemand = X[order_x];
return(mdemand);
}
/**
 * Calculate utility for all J goods
 * @return utility value
 */
real ComputeUtilJ(real inc, vector quant_j, vector price_j,
					vector psi_j, vector gamma_j, vector alpha,
					int ngoods, int model_num) {

	real output;
	real util_num; // numeraire
	vector[ngoods] util_j;

	util_num = 1 / alpha[1] * pow(inc -  price_j' * quant_j, alpha[1]);

	if (model_num == 1){
		util_j = psi_j .* gamma_j .* log(quant_j ./ gamma_j + 1);
	} else if (model_num != 1){
		for (n in 1:ngoods)
			util_j[n] = (psi_j[n] * gamma_j[n]) / alpha[n+1] *
			(pow(quant_j[n] / gamma_j[n] + 1, alpha[n+1]) - 1);
	}

	output = util_num + sum(util_j);
return(output);
}
/**
 * Calculate utility for all M consumed goods
 * @return utility value
 */
real ComputeUtilM(int M, real lambda1, vector g_psi_a, vector a_a_1, vector mu_a_a_1,
					vector psi, vector g, vector price, vector d, int model_num){
	real output;
	vector[M] temp;
	temp[1] = g_psi_a[1] * (pow(lambda1, a_a_1[1]) * mu_a_a_1[1] - d[1]); // compute numeraire util
	if (M > 1){
		for (m in 2:M){
			if (model_num == 1)
				temp[m] = psi[m] * g[m] * log(psi[m] / (lambda1 * price[m]));
			 else if (model_num != 1)
			 	temp[m] = g_psi_a[m] * (pow(lambda1, a_a_1[m]) * mu_a_a_1[m] - d[m]);
		}
	}
	output =  sum(temp);
return(output);
}
/**
 * Calculate HicksDemands using general or hybrid approach
 * @return vector of ngood demands
 */
vector HicksianDemand(real util, vector price,
				vector MUzero, vector gamma, vector alpha,
				int ngoods, int algo_gen, int model_num, real tol_l, int max_loop) {

	vector[ngoods+1] hdemand;
	int M = 1; // Indicator of which ordered alteratives (<=M) are being considered
	int exit = 0;
	real lambda1;
	real util_new;
	int order_x[ngoods+1] = CalcGoodOrder(MUzero, ngoods);
	vector[ngoods+1] X = rep_vector(0, ngoods+1); // vector to hold zero demands
	vector[ngoods+1] d = append_row(0, rep_vector(1, ngoods));
	matrix[ngoods+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, ngoods);
	vector[ngoods+1] mu = col(parm_matrix, 1); // obtain mu
	vector[ngoods+1] g = col(parm_matrix, 3); // obtain gamma

	if (algo_gen == 0) { //Hybrid approach to demand simulation (constant alpha's)
		real lambda_num;
		real lambda_den;
		real lambda;
		real alpha_1 = alpha[1];
		vector[ngoods+1] g_psi = g .* mu .* col(parm_matrix, 2); // obtain gamma_psi
		vector[ngoods+1] b;
		vector[ngoods+1] c;

		for (j in 1:ngoods+1)
			b[j] = pow(mu[j], -alpha_1 / (alpha_1 - 1)); // want price/psi so take negative of exponent

		c = g_psi .* b;

		while (exit == 0){
			// Calculate 1/lambda for a given M
			lambda_num = alpha_1 * util + sum(g_psi[1:M]) - g_psi[1] / g[1]; // utility not expenditure and subtract numeriare psi
			lambda_den = sum(c[1:M]);
			lambda = pow(lambda_num / lambda_den, -(alpha_1 - 1) / alpha_1); // new exponent term
			lambda1 = 1 / lambda;   // Compare 1/lambda to MU

			// Compare 1/lambda to baseline utility of the next lowest alternative
			// (M+1). If lambda exceeds this value then all lower-valued
			// alternatives have zero demand.
			if (lambda1 > mu[min(M + 1, ngoods + 1)] || M == ngoods+1){

				// Compute demands (using eq. 12 in Pinjari and Bhat)
				for (m in 1:M)
					X[m] = (pow(inv(lambda * mu[m]), inv(alpha_1 - 1)) - d[m]) * g[m];
				exit = 1;

			} else if (M < ngoods + 1)
				M += 1; // adds one to M
		}
	} else if (algo_gen == 1) {//General approach to demand simulation (het. alpha's)
		real lambda_l;
		real lambda_u;
		vector[ngoods+1] price_ord = col(parm_matrix, 2); // price
		vector[ngoods+1] a = col(parm_matrix, 4);//	alpha
		vector[ngoods+1] psi = mu .* price_ord;
		vector[ngoods+1] g_psi_a = g .* psi ./ a; // gamma*psi/alpha
		vector[ngoods+1] a_a_1 = a ./ (a - 1); // (alpha/(alpha-1))
		vector[ngoods+1] mu_a_a_1; // (1/MUzero)^(alpha/(alpha-1))

		for (j in 1:ngoods+1)
			mu_a_a_1[j] = pow(inv(mu[j]), a_a_1[j]);

		while (exit == 0){
			lambda1 = mu[M + 1];// Calculate lambda1 equal to MUzero(M+1)

			// Calculate new utility
			util_new = ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi, g, price_ord, d, model_num);

			if (util_new >= util || M+1 == ngoods+1){
				if(util_new < util)
					M += 1;
				lambda_l = util_new < util ? 0 : lambda1;
				lambda_u = mu[M];
				lambda1 = (lambda_l + lambda_u) / 2;

				for (n in 1:max_loop){
					real lambda_mid = (lambda_l + lambda_u) / 2;

					util_new = ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi, g, price_ord, d, model_num);

					// Update lambdas's
					if (util_new < util)
						lambda_u = lambda_mid;
					else if (util_new > util)
						lambda_l = lambda_mid;

					lambda1 = (lambda_l + lambda_u) / 2;

					if (fabs((lambda_l - lambda_u) / (lambda_l + lambda_u) * 0.5) < tol_l) break;
				}

			// Compute demands (using eq. 12 in Pinjari and Bhat)
			for (m in 1:M)
				X[m] = ((lambda1 / mu[m]) ^ (1 / (a[m] - 1))  -d[m]) * g[m];

			exit = 1;

			} else if (util_new < util && M+1 < ngoods+1)
				M += 1; // adds one to M
		}
	}
	// This code puts the choices back in their original order and exports demands
	hdemand = X[order_x];
return(hdemand);
}

// Overall code
matrix CalcmdemandOne_rng(real inc, vector price,
						vector psi_sims, vector gamma_sims, vector alpha_sims, real scale_sims,
						int nerrs, int algo_gen, real tol, int max_loop){

	int ngoods = num_elements(gamma_sims);
	int nsims = 1;
	matrix[nsims, ngoods+1] mdemand_out;

	for (sim in 1:nsims){
		vector[ngoods] psi_j = psi_sims;
		vector[ngoods + 1] gamma = append_row(1, gamma_sims);
		vector[ngoods + 1] alpha = alpha_sims;
		real scale = scale_sims;
		vector[ngoods + 1] error[nerrs];
		matrix[nerrs, ngoods + 1] mdemand;
		row_vector[ngoods + 1] mdemand_sims;
		matrix[ngoods + 1, nerrs] mdemand_trans;

		for(err in 1:nerrs)
			for(j in 1:ngoods+1)
				error[err, j] = -log(-log(uniform_rng(0, 1))) * scale; //uniform(0,1) draws

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[ngoods + 1] MUzero_b;
			vector[ngoods + 1] psi_b_err;
			psi_b_err = exp(append_row(0, psi_j) + error[err]);
			MUzero_b = psi_b_err ./ price;

			mdemand[err] = MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
			ngoods, algo_gen, tol, max_loop)';
		}
		mdemand_trans = mdemand';

		for(g in 1:ngoods+1)
  			mdemand_sims[g] = mean(mdemand_trans[g]);

		mdemand_out[sim] = mdemand_sims;
	}

return(mdemand_out);
}

/**
 * Calculate WTP for each individual, simulation, and policy
 * @param inc income
 * @param quant_j non-numeraire consumption amounts
 * @param price full price vector
 * @param price_p_policy list of length npols containing vectors of length ngood
 * @param psi_p_sims list of length sims constaining npols X ngoods matrices of policy psi
 * @param psi_sims matrix of length nsims x ngoods containing baseline psi
 * @param gamma_sims list of length nsims containing vectors of length ngoods
 * @param alpha_sims list of length nsims containing vectors of length ngoods
 * @param scale_sims vector of length nsims
 * @param rest of model parameter options
 * @return Matrix of nsims x npols wtp
 */
matrix CalcWTP_rng(real inc, vector quant_j, vector price,
						vector[] price_p_policy, matrix[] psi_p_sims,
						matrix psi_sims, vector[] gamma_sims, vector[] alpha_sims, vector scale_sims,
						int nerrs, int cond_error, int draw_mlhs,
						int algo_gen, int model_num, real tol, int max_loop){

	int ngoods = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[nsims, npols] wtp;
	real quant_num = inc - quant_j' * price[2:ngoods+1];

	for (sim in 1:nsims){
		vector[ngoods] psi_j = psi_sims[sim]';
		matrix[npols, ngoods] psi_p_policy = psi_p_sims[sim]; //change for no psi_p
//		vector[ngoods + 1] psi_b_err[nerrs]; // keep psi for use in policies
		vector[ngoods + 1] gamma = append_row(1, gamma_sims[sim]);
		vector[ngoods + 1] alpha = alpha_sims[sim];
		real scale = scale_sims[sim];
		vector[ngoods + 1] error[nerrs];
		vector[npols] wtp_policy;
		vector[nerrs] util;

		error = DrawError_rng(quant_num, quant_j, price[2:ngoods+1],
							psi_j, gamma[2:ngoods+1], alpha, scale,
							ngoods, nerrs, cond_error, draw_mlhs);

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[ngoods + 1] mdemand;
			vector[ngoods + 1] MUzero_b;
//			psi_b_err[err] = exp(append_row(0, psi_j) + error[err]);
//			MUzero_b = psi_b_err[err] ./ price;
			vector[ngoods + 1] psi_b_err;
			psi_b_err = exp(append_row(0, psi_j) + error[err]); //change for no psi_p
			MUzero_b = psi_b_err ./ price; //change for no psi_p

			if (cond_error == 1){
				mdemand = append_row(quant_num, quant_j);
			} else if(cond_error == 0)
				mdemand = MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
				ngoods, algo_gen, tol, max_loop);

			util[err] = ComputeUtilJ(inc, mdemand[2:ngoods+1], price[2:ngoods+1],
									psi_b_err[2:ngoods+1], gamma[2:ngoods+1], alpha,
									ngoods, model_num); // add err to psi_b if no psi_p
		}

		for (policy in 1:npols){
			vector[ngoods + 1] price_p = price + price_p_policy[policy]; // add price increase
			vector[ngoods] psi_p = psi_p_policy[policy]'; //change for no psi_p
			vector[nerrs] wtp_err;

			for (err in 1:nerrs){
//				vector[ngoods + 1] MUzero_p = psi_b_err[err] ./ price_p;	// change to psi_p_err for policy
				vector[ngoods + 1] MUzero_p = exp(append_row(0, psi_p) + error[err]) ./ price_p; //change for no psi_p
				vector[ngoods + 1] hdemand;

				hdemand = HicksianDemand(util[err], price_p, MUzero_p, gamma, alpha,
										ngoods, algo_gen, model_num, tol, max_loop);

				wtp_err[err] = inc - price_p' * hdemand;
			}
			wtp_policy[policy] = mean(wtp_err);
		}
	wtp[sim] = wtp_policy';
	}
return(wtp);
}

/**
 * Calculate WTP for each individual, simulation, and policy using only price changes
 * @param inc income
 * @param quant_j non-numeraire consumption amounts
 * @param price full price vector
 * @param psi_p_sims list of length sims constaining npols X ngoods matrices of policy psi
 * @param psi_sims matrix of length nsims x ngoods containing baseline psi
 * @param gamma_sims list of length nsims containing vectors of length ngoods
 * @param alpha_sims list of length nsims containing vectors of length ngoods
 * @param scale_sims vector of length nsims
 * @param rest of model parameter options
 * @return Matrix of nsims x npols wtp
 */
matrix CalcWTPPriceOnly_rng(real inc, vector quant_j, vector price,
						vector[] price_p_policy,
						matrix psi_sims, vector[] gamma_sims, vector[] alpha_sims, vector scale_sims,
						int nerrs, int cond_error, int draw_mlhs, int algo_gen, int model_num, real tol, int max_loop){

	int ngoods = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[nsims, npols] wtp;
	real quant_num = inc - quant_j' * price[2:ngoods+1];

	for (sim in 1:nsims){
		vector[ngoods] psi_j = psi_sims[sim]';
		vector[ngoods + 1] psi_b_err[nerrs]; // keep psi for use in policies
		vector[ngoods + 1] gamma = append_row(1, gamma_sims[sim]);
		vector[ngoods + 1] alpha = alpha_sims[sim];
		real scale = scale_sims[sim];
		vector[ngoods + 1] error[nerrs];
		vector[npols] wtp_policy;
		vector[nerrs] util;

		error = DrawError_rng(quant_num, quant_j, price[2:ngoods+1],
							psi_j, gamma[2:ngoods+1], alpha, scale,
							ngoods, nerrs, cond_error, draw_mlhs);

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[ngoods + 1] mdemand;
			vector[ngoods + 1] MUzero_b;
			psi_b_err[err] = exp(append_row(0, psi_j) + error[err]); //change for no psi_p
			MUzero_b = psi_b_err[err] ./ price; //change for no psi_p

			if (cond_error == 1){
				mdemand = append_row(quant_num, quant_j);
			} else if(cond_error == 0)
				mdemand = MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
				ngoods, algo_gen, tol, max_loop);

			util[err] = ComputeUtilJ(inc, mdemand[2:ngoods+1], price[2:ngoods+1],
									psi_b_err[err, 2:ngoods+1], gamma[2:ngoods+1], alpha,
									ngoods, model_num);
		}

		for (policy in 1:npols){
			vector[ngoods + 1] price_p = price + price_p_policy[policy]; // add price increase
			vector[nerrs] wtp_err;

			for (err in 1:nerrs){
				vector[ngoods + 1] MUzero_p = psi_b_err[err] ./ price_p;
				vector[ngoods + 1] hdemand;

				hdemand = HicksianDemand(util[err], price_p, MUzero_p, gamma, alpha,
										ngoods, algo_gen, model_num, tol, max_loop);

				wtp_err[err] = inc - price_p' * hdemand;
			}
			wtp_policy[policy] = mean(wtp_err);
		}
	wtp[sim] = wtp_policy';
	}
return(wtp);
}

matrix[] CalcMarshallianDemand_rng(real inc, vector quant_j, vector price,
						vector[] price_p_policy, matrix[] psi_p_sims,
						matrix psi_sims, vector[] gamma_sims, vector[] alpha_sims, vector scale_sims,
						int nerrs, int cond_error, int draw_mlhs,
						int algo_gen, int model_num, real tol, int max_loop){

	int ngoods = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[npols, ngoods+1] mdemand_out[nsims];
	real quant_num = inc - quant_j' * price[2:ngoods+1];

	for (sim in 1:nsims){
		vector[ngoods] psi_j = psi_sims[sim]';
		matrix[npols, ngoods] psi_p_policy = psi_p_sims[sim]; //change for no psi_p
//		vector[ngoods + 1] psi_b_err[nerrs]; // keep psi for use in policies
		vector[ngoods + 1] gamma = append_row(1, gamma_sims[sim]);
		vector[ngoods + 1] alpha = alpha_sims[sim];
		real scale = scale_sims[sim];
		vector[ngoods + 1] error[nerrs];
		matrix[npols, ngoods + 1] mdemand_pols;
		vector[nerrs] util;

		error = DrawError_rng(quant_num, quant_j, price[2:ngoods+1],
							psi_j, gamma[2:ngoods+1], alpha, scale,
							ngoods, nerrs, cond_error, draw_mlhs);

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[ngoods + 1] mdemand_util;
			vector[ngoods + 1] MUzero_b;
//			psi_b_err[err] = exp(append_row(0, psi_j) + error[err]);
//			MUzero_b = psi_b_err[err] ./ price;
			vector[ngoods + 1] psi_b_err;
			psi_b_err = exp(append_row(0, psi_j) + error[err]); //change for no psi_p
			MUzero_b = psi_b_err ./ price; //change for no psi_p

			if (cond_error == 1){
				mdemand_util= append_row(quant_num, quant_j);
			} else if(cond_error == 0)
				mdemand_util = MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
				ngoods, algo_gen, tol, max_loop);

			util[err] = ComputeUtilJ(inc, mdemand_util[2:ngoods+1], price[2:ngoods+1],
									psi_b_err[2:ngoods+1], gamma[2:ngoods+1], alpha,
									ngoods, model_num); // add err to psi_b if no psi_p
		}

		for (policy in 1:npols){
			vector[ngoods + 1] price_p = price + price_p_policy[policy]; // add price increase
			vector[ngoods] psi_p = psi_p_policy[policy]'; //change for no psi_p
			row_vector[ngoods + 1] mdemand_g;
			matrix[nerrs, ngoods + 1] mdemand_p;
			matrix[ngoods + 1, nerrs] mdemand_trans;

			for (err in 1:nerrs){
				vector[ngoods + 1] MUzero_p = exp(append_row(0, psi_p) + error[err]) ./ price_p;

				mdemand_p[err] = MarshallianDemand(inc, price, MUzero_p, gamma, alpha,
					ngoods, algo_gen, tol, max_loop)';
			}
		mdemand_trans = mdemand_p';

		for(g in 1:ngoods+1)
  			mdemand_g[g] = mean(mdemand_trans[g]);

		mdemand_pols[npols] = mdemand_g;
		}
	mdemand_out[sim] = mdemand_pols;
	}
return(mdemand_out);
}
}

