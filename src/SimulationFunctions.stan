//Code for MDCEV Simulation Functions <- "
functions {

vector Shuffle_rng(vector inv_temp, int nerrs){

	vector[nerrs] out;
	vector[nerrs] temp1 = to_vector(uniform_rng(rep_vector(0, nerrs), 1));
	out = inv_temp[sort_indices_asc(temp1)];

return(out);
}

/**Generate random draws using Modified Latin Hypercube Sampling algorithm or uniform
#' Algorithm described in
#' Hess, S., Train, K., and Polak, J. (2006) Transportation Research 40B, 147 - 163.
**/
vector DrawMlhs_rng(int nerrs, int draw_mlhs){

	vector[nerrs] error;
	vector[nerrs] temp0 = rep_vector(0, nerrs);

	if(draw_mlhs == 0){
		error = to_vector(uniform_rng(temp0, 1));
	} else if(draw_mlhs == 1){
		int temp1[nerrs];
		vector[nerrs] temp;

		for (err in 1:nerrs)
			temp1[err] = err - 1;

		temp = to_vector(temp1) / nerrs;
		error = Shuffle_rng(temp + to_vector(uniform_rng(temp0, 1))/ nerrs, nerrs);
	}
return(error);
}

vector[] DrawError_rng(real quant_num, vector quant_j, vector price_j,
					vector psi_j, vector phi_j, vector gamma_j, vector alpha, real scale,
					int model_num, int nalts, int nerrs, int cond_error, int draw_mlhs){

	vector[nalts+1] out[nerrs];
	matrix[nerrs, nalts+1] error;
	vector[nerrs] temp0 = rep_vector(0, nerrs);

	if (cond_error == 0){	// Unconditional
	//	error[1] = rep_row_vector(0, nerrs);
		for(j in 1:nalts+1)
			error[ , j] = -log(-log(DrawMlhs_rng(nerrs, draw_mlhs))) * scale;

	} else if (cond_error == 1){	// Conditional
		// Cacluate the demand function, g
		vector[nalts + 1] cond_demand = append_row(quant_num, quant_j);
		vector[nalts + 1] ek;
		vector[nalts] vk;
		real v1 = (alpha[1] - 1) * log(quant_num);
		// compute vk and v1
			if (model_num != 5)
				vk = psi_j + (alpha[2:nalts+1] - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
			else if (model_num == 5)
				vk = psi_j - log(phi_j .* quant_j + gamma_j) + log(phi_j) - log(price_j);
			// ek = v1 - vk and assume error term is zero for numeraire
			ek = append_row(0, (v1 - vk) / scale);
		// Calculate errors
		// For unvisited alternatives, draw from truncated multivariate logistic distribution
		for(j in 1:nalts+1)
			error[ , j] = cond_demand[j] > 0 ? rep_vector(ek[j] * scale, nerrs) :
					-log(-log(DrawMlhs_rng(nerrs, draw_mlhs) * exp(-exp(-ek[j])))) * scale;
	}

	for(err in 1:nerrs)
		out[err] = error[err]';
return(out);
}

/**
 * Returns an ranking vector of non-numeraire alts by MUzero
 * @return integer vector of ordered alts by MUzero
 */
int[] CalcAltOrder(vector MUzero, int nalts) {

	int order_x[nalts+1];
	vector[nalts] ord_alts;
	int order_MU[nalts] = sort_indices_desc(MUzero[2:nalts+1]);

	for (j in 1:nalts)
		ord_alts[j] = j;

	order_x = sort_indices_asc(append_row(1.0, to_vector(ord_alts[order_MU] + 1)));
return(order_x);
}
/**
 * Sort matrix of parameters/prices
 * @return Matrix of nalts x 4
 */
matrix SortParmMatrix(vector MUzero, vector price, vector gamma, vector alpha_phi, int nalts) {

	matrix[nalts+1, 4] parm_matrix;
	vector[nalts] MU_j = MUzero[2:nalts+1];
	vector[nalts] price_j = price[2:nalts+1];
	vector[nalts] gamma_j = gamma[2:nalts+1];
	vector[nalts] alpha_phi_j = alpha_phi[2:nalts+1];
	int order_MU[nalts] = sort_indices_desc(MU_j);	// find ranking of non-numeraire alts by MUzero

	parm_matrix = append_col(append_row(MUzero[1], MU_j[order_MU]),
	append_col(append_row(price[1], price_j[order_MU]),
	append_col(append_row(gamma[1], gamma_j[order_MU]), append_row(alpha_phi[1], alpha_phi_j[order_MU]))));
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

real ComputeKtE(int M, real lambda, vector mu, vector g_price__phi, real alpha_1){
	real output;
	vector[M] temp;
	temp[1] = pow(lambda, inv(alpha_1 - 1));
	if(M > 1){
		for (m in 2:M)
			temp[m] = (mu[m] / lambda - 1) * g_price__phi[m];
	}
	output =  sum(temp);
return(output);
}

/**
 * Calculate MarshDemands using general or hybrid approach
 * @return vector of nalt demands
 */
vector MarshallianDemand(real income, vector price, vector MUzero, vector phi, vector gamma, vector alpha,
							int nalts, int algo_gen, int model_num, real tol_e, int max_loop) {

	vector[nalts+1] mdemand;
	real lambda;
	int M = 1; // Indicator of which ordered alternatives (<=M) are being considered
	int exit = 0;
	real E;
	int order_x[nalts+1] = CalcAltOrder(MUzero, nalts);
	vector[nalts+1] X = rep_vector(0, nalts+1);
	vector[nalts+1] d = append_row(0, rep_vector(1, nalts));

	if (algo_gen == 0) { //Hybrid
		matrix[nalts+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
		vector[nalts+1] mu = col(parm_matrix, 1); // obtain mu
		vector[nalts+1] g = col(parm_matrix, 3); // obtain gamma
		vector[nalts+1] g_price = g .* col(parm_matrix, 2);
		real alpha_1 = alpha[1];
		vector[nalts+1] b;

		for (j in 1:nalts+1)
			b[j] = g_price[j] * pow(mu[j], inv(1 - alpha_1));

		while (exit == 0){
			// Calculate lambda equal to MUzero(M+1)
			real lambda_num = income + sum(g_price[1:M]) - 1; // minus one for numeraire
			real lambda_den = sum(b[1:M]);
			lambda = pow(lambda_num / lambda_den, alpha_1 - 1);

			//Compare lambda to baseline utility of the next lowest alternative
			//(M+1). If lambda exceeds this value then all lower-valued
			//alternatives have zero demand.
			if (lambda > mu[min(M + 1, nalts + 1)] || M == nalts+1){
				// Compute demands (using eq. 12 in Pinjari and Bhat)
				for (m in 1:M)
					X[m] = (pow(lambda / mu[m], inv(alpha_1 - 1)) - d[m]) * g[m];
				exit = 1;

			} else if ( M < nalts + 1)
				M += 1; // adds one to M
			}

	} else if (algo_gen == 1) {//	General
		real lambda_l;
		real lambda_u;
		matrix[nalts+1, 4] parm_matrix;
		vector[nalts+1] mu; // obtain mu
		vector[nalts+1] c;
		vector[nalts+1] b;
		real alpha_1;
		vector[nalts+1] g__phi; // obtain gamma/phi
		vector[nalts+1] g_price__phi; // obtain gamma*price/phi

		if (model_num < 5){
			parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
			mu = col(parm_matrix, 1); // obtain mu = psi/price
			g__phi = col(parm_matrix, 3); // obtain gamma
			b = inv(col(parm_matrix, 4) - 1);
			for (j in 1:nalts+1)
				c[j] = mu[j] ^ b[j];
		} else if (model_num == 5){
			parm_matrix = SortParmMatrix(MUzero, price, gamma, phi, nalts);
			alpha_1 = alpha[1];
			mu = col(parm_matrix, 1); // obtain mu = psi/price  or= psi*phi/(gamma*price) for kt_ee
			g__phi = col(parm_matrix, 3) ./ col(parm_matrix, 4); // obtain gamma/phi
		}
		g_price__phi = g__phi .* col(parm_matrix, 2); // obtain gamma*price/phi

		while (exit == 0){
			lambda = mu[M + 1];	// Calculate lambda equal to MUzero(M+1)

			if (model_num < 5)
				E = ComputeE(M, lambda, g_price__phi, b, c, d);

			else if (model_num == 5)
				E = ComputeKtE(M, lambda, mu, g_price__phi, alpha_1);

			if (E >= income || M+1 == nalts+1){
				if (E < income)
					M += 1;

				lambda_l = E < income ? 0 : lambda;
				lambda_u = mu[M];
				lambda = (lambda_l + lambda_u) / 2;

				for (n in 1:max_loop){
					real lambda_mid = (lambda_l + lambda_u) / 2;

					if (model_num < 5)
						E = ComputeE(M, lambda, g_price__phi, b, c, d);
					else if (model_num == 5)
						E = ComputeKtE(M, lambda, mu, g_price__phi, alpha_1);

					// Update lambdas's
					if (E < income)
						lambda_u = lambda_mid;
					else if (E > income)
						lambda_l = lambda_mid;

					lambda = (lambda_l + lambda_u) / 2;

					if (fabs((E - income) / (E + income) * 0.5) < tol_e) break;
				}
				// Compute demands (using eq. 12 in Pinjari and Bhat)
				if (model_num < 5){
					for (m in 1:M)
						X[m] = ((lambda / mu[m]) ^ b[m] -d[m]) * g__phi[m];
				}
				else if (model_num == 5){
					X[1] = pow(lambda, inv(alpha_1 - 1));

					if(M > 1){
						for (m in 2:M)
							X[m] = (mu[m] / lambda - 1) * g__phi[m];
					}
				}

				exit = 1;

			} else if (E < income && M + 1 < nalts + 1)
				M += 1; // adds one to M
		}
	}
	// This code puts the choices back in their original order and exports demands
	mdemand = X[order_x];
return(mdemand);
}

/**
 * Calculate utility for all J alts
 * @return utility value
 */
real ComputeUtilJ(real income, vector quant_j, vector price_j,
					vector psi_j, vector phi_j, vector gamma_j, vector alpha,
					int nalts, int model_num) {

	real output;
	real util_num; // numeraire
	vector[nalts] util_j;

	if (model_num == 4)
		util_num = log(income -  price_j' * quant_j);
	else
		util_num = pow(income -  price_j' * quant_j, alpha[1]) / alpha[1];

	if (model_num == 1 || model_num == 4){
		util_j = psi_j .* gamma_j .* log(quant_j ./ gamma_j + 1);
	} else if (model_num == 5){
		util_j = psi_j .* log((phi_j .* quant_j ./ gamma_j) +1);
	} else {
		for (n in 1:nalts)
			util_j[n] = (psi_j[n] * gamma_j[n]) / alpha[n+1] *
			(pow(quant_j[n] / gamma_j[n] + 1, alpha[n+1]) - 1);
	}

	output = util_num + sum(util_j);
return(output);
}

/**
 * Calculate utility for all M consumed alts
 * @return utility value
 */
real ComputeUtilM(int M, real lambda1, vector g_psi_a, vector a_a_1, vector mu_a_a_1,
					vector psi, vector g, vector price, vector d, int model_num){
	real output;
	vector[M] temp = rep_vector(0, M);
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

real ComputeKtUtilM(int M, real lambda1, vector psi, vector mu, real alpha_1){
	// use mu in logs
	real output;
	vector[M] temp = rep_vector(0, M);
	temp[1] = pow(lambda1, alpha_1 / (alpha_1 - 1)) / alpha_1; // compute numeraire util
	if (M > 1){
		for (m in 2:M)
			temp[m] = psi[m] * log(mu[m] / lambda1);
	}
	output =  sum(temp);
return(output);
}

/**
 * Calculate HicksDemands using general or hybrid approach
 * @return vector of nalt demands
 */
vector HicksianDemand(real util, vector price,
				vector MUzero, vector phi, vector gamma, vector alpha,
				int nalts, int algo_gen, int model_num, real tol_l, int max_loop) {

	vector[nalts+1] hdemand;
	int M = 1; // Indicator of which ordered alteratives (<=M) are being considered
	int exit = 0;
	real lambda1;
	real util_new;
	int order_x[nalts+1] = CalcAltOrder(MUzero, nalts);
	vector[nalts+1] X = rep_vector(0, nalts+1); // vector to hold zero demands
	vector[nalts+1] d = append_row(0, rep_vector(1, nalts));

	if (algo_gen == 0) { //Hybrid approach to demand simulation (constant alpha's)
		matrix[nalts+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
		real alpha_1 = alpha[1]; // all alpha's are equal
		vector[nalts+1] mu = col(parm_matrix, 1); // obtain mu
		vector[nalts+1] g = col(parm_matrix, 3); // obtain gamma
		vector[nalts+1] g_psi = g .* mu .* col(parm_matrix, 2); // obtain gamma_psi
		vector[nalts+1] c;

		if (model_num == 3){
			vector[nalts+1] b;
			for (j in 1:nalts+1)
				b[j] = pow(mu[j], -alpha_1 / (alpha_1 - 1)); // want price/psi so take negative of exponent

			c = g_psi .* b;

		} if (model_num == 4){
			c = g_psi .* log(mu);
		}

		while (exit == 0){
			// Calculate 1/lambda for a given M
			if (model_num == 3){
				real lambda_num = alpha_1 * util + sum(g_psi[1:M]) - g_psi[1]; // utility not expenditure and subtract numeriare psi
				real lambda_den = sum(c[1:M]);
				lambda1 = pow(lambda_num / lambda_den, (alpha_1 - 1) / alpha_1); // create 1/lambda term
			} else if (model_num == 4){
				real lambda_num = util - sum(c[1:M]);
				real lambda_den = sum(g_psi[1:M]);
				lambda1 = inv(exp(lambda_num / lambda_den)); // create 1/lambda term = 1/exp(expression)
			}
			// Compare 1/lambda to baseline utility of the next lowest alternative
			// (M+1). If lambda exceeds this value then all lower-valued
			// alternatives have zero demand.
			if (lambda1 > mu[min(M + 1, nalts + 1)] || M == nalts+1){

				// Compute demands (using eq. 12 in Pinjari and Bhat)
				for (m in 1:M)
					X[m] = (pow(lambda1 / mu[m], inv(alpha_1 - 1)) - d[m]) * g[m];
				exit = 1;

			} else if (M < nalts + 1)
				M += 1; // adds one to M
		}
	} else if (algo_gen == 1) {//General approach to demand simulation (het. alpha's)
		real lambda_l;
		real lambda_u;
		if (model_num < 5){
			matrix[nalts+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
			vector[nalts+1] mu = col(parm_matrix, 1); // obtain mu
			vector[nalts+1] g = col(parm_matrix, 3); // obtain gamma
			vector[nalts+1] price_ord = col(parm_matrix, 2); // price
			vector[nalts+1] a = col(parm_matrix, 4);//	alpha
			vector[nalts+1] psi = mu .* price_ord;
			vector[nalts+1] g_psi_a = g .* psi ./ a; // gamma*psi/alpha
			vector[nalts+1] a_a_1 = a ./ (a - 1); // (alpha/(alpha-1))
			vector[nalts+1] mu_a_a_1; // (1/MUzero)^(alpha/(alpha-1))

			for (j in 1:nalts+1)
				mu_a_a_1[j] = pow(inv(mu[j]), a_a_1[j]);

			while (exit == 0){
				lambda1 = mu[M + 1];// Calculate lambda1 equal to MUzero(M+1)

				// Calculate new utility
				util_new = ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi, g, price_ord, d, model_num);

				if (util_new >= util || M+1 == nalts+1){
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
					X[m] = (pow(lambda1 / mu[m], inv(a[m] - 1)) - d[m]) * g[m];
				exit = 1;

				} else if (util_new < util && M+1 < nalts+1)
					M += 1; // adds one to M
			}
		} else if (model_num == 5){
			matrix[nalts+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, phi, nalts);
			real alpha_1 = alpha[1];
			vector[nalts+1] mu = col(parm_matrix, 1); // obtain mu = psi*phi/(gamma*price)
			vector[nalts+1] g__phi = col(parm_matrix, 3) ./ col(parm_matrix, 4); // gamma/phi
			vector[nalts+1] psi_ord = mu .* col(parm_matrix, 2) .* g__phi;

			while (exit == 0){
				lambda1 = mu[M + 1];// Calculate lambda1 equal to MUzero(M+1)

				// Calculate new utility
				util_new = ComputeKtUtilM(M, lambda1, psi_ord, mu, alpha_1);

				if (util_new >= util || M+1 == nalts+1){
					if(util_new < util)
						M += 1;
					lambda_l = util_new < util ? 0 : lambda1;
					lambda_u = mu[M];
					lambda1 = (lambda_l + lambda_u) / 2;

					for (n in 1:max_loop){
						real lambda_mid = (lambda_l + lambda_u) / 2;

						util_new = ComputeKtUtilM(M, lambda1, psi_ord, mu, alpha_1);

						// Update lambdas's
						if (util_new < util)
							lambda_u = lambda_mid;
						else if (util_new > util)
							lambda_l = lambda_mid;

						lambda1 = (lambda_l + lambda_u) / 2;

						if (fabs((lambda_l - lambda_u) / (lambda_l + lambda_u) * 0.5) < tol_l) break;
					}
				// Compute demands (using modified version of eq. 12 in Pinjari and Bhat)
					X[1] = pow(lambda1, inv(alpha_1 - 1));
					if(M > 1){
						for (m in 2:M)
							X[m] = (mu[m] / lambda1 - 1) * g__phi[m];
					}
				exit = 1;

				} else if (util_new < util && M+1 < nalts+1)
					M += 1; // adds one to M
			}
		}
	}
	// This code puts the choices back in their original order and exports demands
	hdemand = X[order_x];

return(hdemand);
}

/**
 * Calculate WTP for each individual, simulation, and policy
 * @param income income
 * @param quant_j non-numeraire consumption amounts
 * @param price full price vector
 * @param price_p_policy list of length npols containing vectors of length nalt
 * @param psi_p_sims list of length sims constaining npols X nalts matrices of policy psi
 * @param psi_sims matrix of length nsims x nalts containing baseline psi
 * @param gamma_sims list of length nsims containing vectors of length nalts
 * @param alpha_sims list of length nsims containing vectors of length nalts
 * @param scale_sims vector of length nsims
 * @param rest of model parameter options
 * @return Matrix of nsims x npols wtp
 */
matrix CalcWTP_rng(real income, vector quant_j, vector price,
						vector[] price_p_policy, matrix[] psi_p_sims, matrix[] phi_p_sims,
						matrix psi_sims, matrix phi_sims, matrix gamma_sims, matrix alpha_sims, vector scale_sims,
						int nerrs, int cond_error, int draw_mlhs,
						int algo_gen, int model_num, int price_change_only, real tol, int max_loop){

	int nalts = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[nsims, npols] wtp;
	real quant_num = income - quant_j' * price[2:nalts+1];

	for (sim in 1:nsims){
		vector[nalts] psi_j = psi_sims[sim]';
		vector[nalts + 1] phi;
		matrix[npols, nalts] psi_p_policy;
		vector[nalts + 1] gamma = append_row(1, gamma_sims[sim]');
		vector[nalts + 1] alpha = alpha_sims[sim]';
		real scale = scale_sims[sim];
		vector[nalts + 1] error[nerrs];
		vector[npols] wtp_policy;
		vector[nerrs] util;

		if (price_change_only == 0)
			psi_p_policy = psi_p_sims[sim];

		if (model_num < 5)
			phi = rep_vector(1, nalts + 1);
		else if (model_num == 5)
			phi = append_row(1,phi_sims[sim]');

		error = DrawError_rng(quant_num, quant_j, price[2:nalts+1],
							psi_j, phi[2:nalts+1], gamma[2:nalts+1], alpha, scale,
							model_num, nalts, nerrs, cond_error, draw_mlhs);

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[nalts + 1] mdemand;
			vector[nalts + 1] MUzero_b;
			vector[nalts + 1] psi_b_err = exp(append_row(0, psi_j) + error[err]);
			MUzero_b = psi_b_err ./ price;
			if (model_num == 5)
				MUzero_b = MUzero_b .* phi ./ gamma; // specific form for kt_ee model

			if (cond_error == 1){
				mdemand = append_row(quant_num, quant_j);
			} else if(cond_error == 0)
				mdemand = MarshallianDemand(income, price, MUzero_b, phi, gamma, alpha,
											nalts, algo_gen, model_num, tol, max_loop);

			util[err] = ComputeUtilJ(income, mdemand[2:nalts+1], price[2:nalts+1],
									psi_b_err[2:nalts+1], phi[2:nalts+1], gamma[2:nalts+1], alpha,
									nalts, model_num);
		}

		for (policy in 1:npols){
			vector[nalts + 1] price_p = price + price_p_policy[policy]; // add price increase
			vector[nalts] psi_p;
			vector[nerrs] wtp_err;

			if (price_change_only == 0)
				psi_p = psi_p_policy[policy]';
			else
				psi_p = psi_j;

			for (err in 1:nerrs){
				vector[nalts + 1] hdemand;
				vector[nalts + 1] MUzero_p = exp(append_row(0, psi_p) + error[err]) ./ price_p; //change for no psi_p;
				if (model_num == 5)
					MUzero_p = MUzero_p .* phi ./ gamma; // specific form for kt_ee model

				hdemand = HicksianDemand(util[err], price_p, MUzero_p, phi, gamma, alpha,
											nalts, algo_gen, model_num, tol, max_loop);

				wtp_err[err] = income - price_p' * hdemand;
			}
			wtp_policy[policy] = mean(wtp_err);
		}
	wtp[sim] = wtp_policy';
	}
return(wtp);
}

matrix[] CalcMarshallianDemand_rng(real income, vector quant_j, vector price,
						vector[] price_p_policy, matrix[] psi_p_sims, matrix[] phi_p_sims,
						matrix psi_sims, matrix phi_sims, matrix gamma_sims, matrix alpha_sims, vector scale_sims,
						int nerrs, int cond_error, int draw_mlhs,
						int algo_gen, int model_num, int price_change_only, real tol, int max_loop){

	int nalts = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[npols, nalts+1] mdemand_out[nsims];
	real quant_num = income - quant_j' * price[2:nalts+1];

	for (sim in 1:nsims){
		vector[nalts] psi_j = psi_sims[sim]';
		vector[nalts + 1] phi;
		matrix[npols, nalts] psi_p_policy;
		vector[nalts + 1] gamma = append_row(1, gamma_sims[sim]');
		vector[nalts + 1] alpha = alpha_sims[sim]';
		real scale = scale_sims[sim];
		vector[nalts + 1] error[nerrs];
		matrix[npols, nalts + 1] mdemand_pols;
		vector[nerrs] util;

		if (price_change_only == 0)
			psi_p_policy = psi_p_sims[sim];

		if (model_num < 5)
			phi = rep_vector(1, nalts + 1);
		else if (model_num == 5)
			phi = append_row(1,phi_sims[sim]');

		error = DrawError_rng(quant_num, quant_j, price[2:nalts+1],
							psi_j, phi[2:nalts+1], gamma[2:nalts+1], alpha, scale,
							model_num, nalts, nerrs, cond_error, draw_mlhs);

		// Compute Marshallian demands and baseline utility
		for (err in 1:nerrs){
			vector[nalts + 1] mdemand_util;
			vector[nalts + 1] MUzero_b;
			vector[nalts + 1] psi_b_err = exp(append_row(0, psi_j) + error[err]);
			MUzero_b = psi_b_err ./ price;
			if (model_num == 5)
				MUzero_b = MUzero_b .* phi ./ gamma; // specific form for kt_ee model

			if (cond_error == 1){
				mdemand_util= append_row(quant_num, quant_j);
			} else if(cond_error == 0)
				mdemand_util = MarshallianDemand(income, price, MUzero_b, phi, gamma, alpha,
											nalts, algo_gen, model_num, tol, max_loop);

			util[err] = ComputeUtilJ(income, mdemand_util[2:nalts+1], price[2:nalts+1],
									psi_b_err[2:nalts+1], phi[2:nalts+1], gamma[2:nalts+1], alpha,
									nalts, model_num); // add err to psi_b if no psi_p
		}

		for (policy in 1:npols){
			vector[nalts + 1] price_p = price + price_p_policy[policy]; // add price increase
			vector[nalts] psi_p; //change for no psi_p
			vector[nalts + 1] mdemand_p = rep_vector(0, nalts + 1);

			if (price_change_only == 0)
				psi_p = psi_p_policy[policy]';
			else
				psi_p = psi_j;

			for (err in 1:nerrs){
				vector[nalts + 1] MUzero_p = exp(append_row(0, psi_p) + error[err]) ./ price_p;
				if (model_num == 5)
					MUzero_p = MUzero_p .* phi ./ gamma; //change for no psi_p

				mdemand_p = mdemand_p + MarshallianDemand(income, price_p, MUzero_p, phi, gamma, alpha,
								nalts, algo_gen, model_num, tol, max_loop) / nerrs; // take average
			}
		mdemand_pols[policy] = mdemand_p';
		}
	mdemand_out[sim] = mdemand_pols;
	}
return(mdemand_out);
}

// CalcmdemandOne_rng
// Calculates baseline Marshallian demand for with only one simulations and no policies
// Used to simulate data
vector CalcmdemandOne_rng(real income, vector price,
						vector psi_j, vector phi_j, vector gamma_j, vector alpha, real scale,
						int nerrs, int model_num, int algo_gen, real tol, int max_loop){

	int nalts = num_elements(price) - 1; // subtract numeraire
	vector[nalts + 1] mdemand = rep_vector(0, nalts + 1);
	vector[nalts + 1] gamma = append_row(1, gamma_j);
	vector[nalts + 1] phi;
	vector[nalts + 1] error[nerrs];

	if (model_num < 5)
		phi = rep_vector(1, nalts + 1);
	else if (model_num == 5)
		phi = append_row(1, phi_j);

	for(err in 1:nerrs)
		for(g in 1:nalts+1)
			error[err, g] = -log(-log(uniform_rng(0, 1))) * scale; //uniform(0,1) draws

	// Compute Marshallian demands and baseline utility
	for (err in 1:nerrs){
		vector[nalts + 1] MUzero_b = exp(append_row(0, psi_j) + error[err]) ./ price;
		if (model_num == 5)
			MUzero_b = MUzero_b .* phi ./ gamma; //change for no psi_p

		mdemand = mdemand + MarshallianDemand(income, price, MUzero_b, phi, gamma, alpha,
						nalts, algo_gen, model_num, tol, max_loop) / nerrs; // take average
	}

return(mdemand);
}

// end of functions
}
