// MDCEV Simulation Functions
//
// Functions (internal helpers):
//   Shuffle_rng          - shuffle a vector using a random permutation
//   DrawMlhs_rng         - draw errors via MLHS or iid uniform sampling
//   DrawError_rng        - draw unconditional or conditional Gumbel errors
//   CalcAltOrder         - rank non-numeraire alts by MUzero
//   SortParmMatrix       - sort parameter matrix by descending MUzero
//   ComputeE             - compute expenditure sum for bisection
//   ComputeUtilM         - compute utility over M consumed alternatives
//   ComputeKtUtilM       - compute KT utility over M consumed alternatives
//   MarshallianDemand    - compute Marshallian demands (hybrid or general)
//   ComputeUtilJ         - compute utility over all J alternatives
//   HicksianDemand       - compute Hicksian demands (hybrid or general)
//
// Functions exported to R via Rcpp:
//   CalcWTP_rng              - WTP matrix (nsims x npols)
//   CalcMarshallianDemand_rng - demand array (nsims x npols x nalts+1)
//   CalcmdemandOne_rng       - baseline Marshallian demand for data simulation
//
// model_num encoding:
//   1 = gamma   (alpha = 0 for non-numeraire)
//   2 = alpha   (gamma = 1)
//   3 = hybrid  (single alpha, alternative-specific gammas)
//   4 = hybrid0 (alpha = 0; most common MDCEV)
//   5 = kt_ee   (environmental economics KT specification; adds phi parameters)
//
// algo_gen encoding:
//   0 = hybrid closed-form (constant alpha; Pinjari-Bhat algorithm)
//   1 = general bisection  (heterogeneous alpha or kt_ee)
functions {

/**
 * Shuffle a vector using a random permutation
 * @param inv_temp vector to shuffle
 * @param nerrs    length of vector
 * @return permuted vector
 */
vector Shuffle_rng(vector inv_temp, int nerrs){

	vector[nerrs] out;
	vector[nerrs] temp1 = to_vector(uniform_rng(rep_vector(0, nerrs), 1));
	out = inv_temp[sort_indices_asc(temp1)];

return out;
}

/**
 * Generate random draws using Modified Latin Hypercube Sampling or iid uniform.
 * Algorithm described in Hess, S., Train, K., and Polak, J. (2006)
 * Transportation Research 40B, 147-163.
 * @param nerrs     number of draws
 * @param draw_mlhs 0 = iid uniform, 1 = MLHS
 * @return vector of nerrs draws in [0,1)
 */
vector DrawMlhs_rng(int nerrs, int draw_mlhs){

	vector[nerrs] error;
	vector[nerrs] temp0 = rep_vector(0, nerrs);

	if(draw_mlhs == 0){
		error = to_vector(uniform_rng(temp0, 1));
	} else {
		vector[nerrs] temp = to_vector(linspaced_int_array(nerrs, 0, nerrs-1)) / nerrs;
		error = Shuffle_rng(temp + to_vector(uniform_rng(temp0, 1))/ nerrs, nerrs);
	}
return error;
}

/**
 * Draw unconditional or conditional Gumbel errors for all alternatives.
 * @param quant_num  observed numeraire quantity
 * @param quant_j    observed non-numeraire quantities (length nalts)
 * @param price_j    non-numeraire prices (length nalts)
 * @param psi_j      non-numeraire psi parameters (length nalts)
 * @param phi_j      non-numeraire phi parameters (kt_ee only; length nalts)
 * @param gamma_j    non-numeraire gamma parameters (length nalts)
 * @param alpha      full alpha vector (length nalts+1)
 * @param scale      scale parameter
 * @param model_num  model type (1-5)
 * @param nalts      number of non-numeraire alternatives
 * @param nerrs      number of error draws
 * @param cond_error 0 = unconditional, 1 = conditional on observed quantities
 * @param draw_mlhs  0 = iid uniform, 1 = MLHS
 * @return array of nerrs vectors of length nalts+1
 */
array[] vector DrawError_rng(real quant_num, vector quant_j, vector price_j,
				vector psi_j, vector phi_j, vector gamma_j, vector alpha, real scale,
				int model_num, int nalts, int nerrs, int cond_error, int draw_mlhs){

	array[nerrs] vector[nalts+1] out;
	matrix[nerrs, nalts+1] error;

	if (cond_error == 0){	// Unconditional
		for(j in 1:(nalts+1))
			error[ , j] = -log(-log(DrawMlhs_rng(nerrs, draw_mlhs))) * scale;

	} else {	// Conditional
		// Calculate the demand function, g
		vector[nalts + 1] cond_demand = append_row(quant_num, quant_j);
		vector[nalts + 1] ek;
		vector[nalts] vk;
		real v1 = (alpha[1] - 1) * log(quant_num);
		// compute vk and v1
			if (model_num != 5)
				vk = psi_j + (alpha[2:(nalts+1)] - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
			else
				vk = psi_j - log(phi_j .* quant_j + gamma_j) + log(phi_j) - log(price_j);
			// ek = v1 - vk and assume error term is zero for numeraire
			ek = append_row(0, (v1 - vk) / scale);
		// Calculate errors
		// For unvisited alternatives, draw from truncated multivariate logistic distribution
		for(j in 1:(nalts+1))
			error[ , j] = cond_demand[j] > 0 ? rep_vector(ek[j] * scale, nerrs) :
					-log(-log(DrawMlhs_rng(nerrs, draw_mlhs) * exp(-exp(-ek[j])))) * scale;
	}

	for(err in 1:nerrs)
		out[err] = error[err]';
return out;
}

/**
 * Returns a ranking vector of non-numeraire alts by MUzero
 * @param MUzero vector of length nalts+1
 * @param nalts  number of non-numeraire alternatives
 * @return integer vector of ordered alts by MUzero (length nalts+1)
 */
array[] int CalcAltOrder(vector MUzero, int nalts) {

	array[nalts+1] int order_x;
	vector[nalts] ord_alts;
	array[nalts] int order_MU = sort_indices_desc(MUzero[2:(nalts+1)]);

	for (j in 1:nalts)
		ord_alts[j] = j;

	order_x = sort_indices_asc(append_row(1.0, to_vector(ord_alts[order_MU] + 1)));
return order_x;
}

/**
 * Sort matrix of parameters/prices by descending MUzero
 * @param MUzero    full MUzero vector (length nalts+1)
 * @param price     full price vector (length nalts+1)
 * @param gamma     full gamma vector (length nalts+1)
 * @param alpha_phi full alpha (or phi for kt_ee) vector (length nalts+1)
 * @param nalts     number of non-numeraire alternatives
 * @return matrix (nalts+1) x 4: columns are mu, price, gamma, alpha_phi sorted by mu desc
 */
matrix SortParmMatrix(vector MUzero, vector price, vector gamma, vector alpha_phi, int nalts) {

	matrix[nalts+1, 4] parm_matrix;
	vector[nalts] MU_j = MUzero[2:(nalts+1)];
	vector[nalts] price_j = price[2:(nalts+1)];
	vector[nalts] gamma_j = gamma[2:(nalts+1)];
	vector[nalts] alpha_phi_j = alpha_phi[2:(nalts+1)];
	array[nalts] int order_MU = sort_indices_desc(MU_j);	// find ranking of non-numeraire alts by MUzero

	parm_matrix = append_col(append_row(MUzero[1], MU_j[order_MU]),
	append_col(append_row(price[1], price_j[order_MU]),
	append_col(append_row(gamma[1], gamma_j[order_MU]), append_row(alpha_phi[1], alpha_phi_j[order_MU]))));
return parm_matrix;
}

/**
 * Compute expenditure sum for bisection: sum_m g_price__phi[m] * (lambda^b[m] / c[m] - d[m])
 * @param M            number of consumed alternatives
 * @param lambda       current lambda value
 * @param g_price__phi gamma*price/phi vector (length >= M)
 * @param b            exponent vector (length >= M)
 * @param c            mu^b vector (length >= M)
 * @param d            indicator vector: 0 for numeraire, 1 for non-numeraire (length >= M)
 * @return scalar expenditure
 */
real ComputeE(int M, real lambda, vector g_price__phi, vector b, vector c, vector d){
	return sum(g_price__phi[1:M] .* (pow(lambda, b[1:M]) ./ c[1:M] - d[1:M]));
}

/**
 * Calculate Marshallian demands using the hybrid (constant-alpha) or general (bisection) approach
 * @param income    individual income
 * @param price     full price vector (length nalts+1)
 * @param MUzero    marginal utility of zero vector (length nalts+1)
 * @param phi       phi vector (kt_ee) or ones (length nalts+1)
 * @param gamma     full gamma vector (length nalts+1)
 * @param alpha     full alpha vector (length nalts+1)
 * @param nalts     number of non-numeraire alternatives
 * @param algo_gen  0 = hybrid closed-form, 1 = general bisection
 * @param model_num model type (1-5)
 * @param tol_e     convergence tolerance for bisection on expenditure
 * @param max_loop  maximum bisection iterations
 * @return vector of nalts+1 demands
 */
vector MarshallianDemand(real income, vector price, vector MUzero, vector phi, vector gamma, vector alpha,
						int nalts, int algo_gen, int model_num, real tol_e, int max_loop) {

	vector[nalts+1] mdemand;
	real lambda;
	int M = 1; // Indicator of which ordered alternatives (<=M) are being considered
	int exit = 0;
	real E;
	array[nalts+1] int order_x = CalcAltOrder(MUzero, nalts);
	vector[nalts+1] X = rep_vector(0, nalts+1);
	vector[nalts+1] d = append_row(0, rep_vector(1, nalts));

	if (algo_gen == 0) { //Hybrid
		matrix[nalts+1, 4] parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
		vector[nalts+1] mu = col(parm_matrix, 1); // obtain mu
		vector[nalts+1] g = col(parm_matrix, 3); // obtain gamma
		vector[nalts+1] g_price = g .* col(parm_matrix, 2);
		real alpha_1 = alpha[1];
		vector[nalts+1] b = g_price .* pow(mu, inv(1 - alpha_1));

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
				X[1:M] = (pow(lambda ./ mu[1:M], inv(alpha_1 - 1)) - d[1:M]) .* g[1:M];
				exit = 1;

			} else if ( M < nalts + 1)
				M += 1; // adds one to M
			}

	} else { //	General
		real lambda_l;
		real lambda_u;
		matrix[nalts+1, 4] parm_matrix;
		vector[nalts+1] mu; // obtain mu
		vector[nalts+1] c;
		vector[nalts+1] b;
		vector[nalts+1] g__phi; // obtain gamma/phi
		vector[nalts+1] g_price__phi; // obtain gamma*price/phi

		if (model_num < 5){
			parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
			b = inv(col(parm_matrix, 4) - 1);
			g__phi = col(parm_matrix, 3); // obtain gamma

		} else {
			parm_matrix = SortParmMatrix(MUzero, price, gamma, phi, nalts);
			b = inv(alpha - 1);
			g__phi = col(parm_matrix, 3) ./ col(parm_matrix, 4); // obtain gamma/phi
		}

		mu = col(parm_matrix, 1); // obtain mu = psi/price
		g_price__phi = g__phi .* col(parm_matrix, 2); // obtain gamma*price/phi
		c = pow(mu, b);

		while (exit == 0){
			lambda = mu[M + 1];	// Calculate lambda equal to MUzero(M+1)

			E = ComputeE(M, lambda, g_price__phi, b, c, d);

			if (E >= income || M+1 == nalts+1){
				if (E < income)
					M += 1;

				lambda_l = E < income ? 0 : lambda;
				lambda_u = mu[M];
				lambda = (lambda_l + lambda_u) / 2;

				for (n in 1:max_loop){
					E = ComputeE(M, lambda, g_price__phi, b, c, d);

					// Update lambda bounds
					if (E < income)
						lambda_u = lambda;
					else if (E > income)
						lambda_l = lambda;

					lambda = (lambda_l + lambda_u) / 2;

					if (abs((E - income) / (E + income) * 0.5) < tol_e) break;
				}
				// Compute demands (using eq. 12 in Pinjari and Bhat)
				X[1:M] = (pow(lambda ./ mu[1:M], b[1:M]) - d[1:M]) .* g__phi[1:M];

				exit = 1;

			} else if (E < income && M + 1 < nalts + 1)
				M += 1; // adds one to M
		}
	}
	// This code puts the choices back in their original order and exports demands
	mdemand = X[order_x];
return mdemand;
}

/**
 * Calculate utility for all J alternatives
 * @param income    individual income
 * @param quant_j   non-numeraire quantities (length nalts)
 * @param price_j   non-numeraire prices (length nalts)
 * @param psi       full psi vector (length nalts+1)
 * @param phi_j     non-numeraire phi parameters (kt_ee only; length nalts)
 * @param gamma_j   non-numeraire gamma parameters (length nalts)
 * @param alpha     full alpha vector (length nalts+1)
 * @param nalts     number of non-numeraire alternatives
 * @param model_num model type (1-5)
 * @return scalar utility
 */
real ComputeUtilJ(real income, vector quant_j, vector price_j,
				vector psi, vector phi_j, vector gamma_j, vector alpha,
				int nalts, int model_num) {
	real output;
	real util_num; // numeraire
	vector[nalts] util_j;

	if (model_num == 4)
		util_num = psi[1] * log(income -  price_j' * quant_j);
	else
		util_num = psi[1] * pow(income -  price_j' * quant_j, alpha[1]) / alpha[1];

	if (model_num == 1 || model_num == 4){
		util_j = psi[2:(nalts+1)] .* gamma_j .* log(quant_j ./ gamma_j + 1);
	} else if (model_num == 5){
		util_j = psi[2:(nalts+1)] .* log((phi_j .* quant_j ./ gamma_j) +1);
	} else {
		for (n in 1:nalts)
			util_j[n] = (psi[n+1] * gamma_j[n]) / alpha[n+1] *
			(pow(quant_j[n] / gamma_j[n] + 1, alpha[n+1]) - 1);
	}

	output = util_num + sum(util_j);
return output;
}

/**
 * Calculate utility for M consumed alternatives (sorted order; used in Hicksian bisection)
 * @param M          number of consumed alternatives
 * @param lambda1    current 1/lambda value
 * @param g_psi_a    gamma*psi/alpha vector (length >= M)
 * @param a_a_1      alpha/(alpha-1) vector (length >= M)
 * @param mu_a_a_1   (1/mu)^(alpha/(alpha-1)) vector (length >= M)
 * @param psi        psi vector (sorted; length >= M)
 * @param g          gamma vector (sorted; length >= M)
 * @param price      price vector (sorted; length >= M)
 * @param d          indicator vector: 0 numeraire, 1 non-numeraire (length >= M)
 * @param model_num  model type (1-5)
 * @return scalar utility
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
			else
			 	temp[m] = g_psi_a[m] * (pow(lambda1, a_a_1[m]) * mu_a_a_1[m] - d[m]);
		}
	}
	output = sum(temp);
return output;
}

/**
 * Calculate KT utility for M consumed alternatives (kt_ee model; used in Hicksian bisection)
 * @param M       number of consumed alternatives
 * @param lambda1 current 1/lambda value
 * @param psi     psi vector (sorted; psi = mu * price * gamma / phi; length >= M)
 * @param mu      MUzero vector (sorted; length >= M)
 * @param alpha_1 common alpha parameter
 * @return scalar utility
 */
real ComputeKtUtilM(int M, real lambda1, vector psi, vector mu, real alpha_1){
	real util_num = psi[1] * pow(lambda1 / mu[1], alpha_1 / (alpha_1 - 1)) / alpha_1;
	return M == 1 ? util_num : util_num + sum(psi[2:M] .* log(mu[2:M] / lambda1));
}

/**
 * Calculate HicksianDemands using general or hybrid approach
 * @param util      target utility level
 * @param price     full price vector (length nalts+1)
 * @param MUzero    marginal utility of zero vector (length nalts+1)
 * @param phi       phi vector (kt_ee) or ones (length nalts+1)
 * @param gamma     full gamma vector (length nalts+1)
 * @param alpha     full alpha vector (length nalts+1)
 * @param nalts     number of non-numeraire alternatives
 * @param algo_gen  0 = hybrid closed-form, 1 = general bisection
 * @param model_num model type (1-5)
 * @param tol_l     convergence tolerance for bisection on lambda
 * @param max_loop  maximum bisection iterations
 * @return vector of nalts+1 demands
 */
vector HicksianDemand(real util, vector price,
			vector MUzero, vector phi, vector gamma, vector alpha,
			int nalts, int algo_gen, int model_num, real tol_l, int max_loop) {

	vector[nalts+1] hdemand;
	int M = 1; // Indicator of which ordered alternatives (<=M) are being considered
	int exit = 0;
	real lambda1;
	real util_new;
	array[nalts+1] int order_x = CalcAltOrder(MUzero, nalts);
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
			c = g_psi .* pow(mu, -alpha_1 / (alpha_1 - 1)); // want price/psi so take negative of exponent

		} else if (model_num == 4){
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
				X[1:M] = (pow(lambda1 ./ mu[1:M], inv(alpha_1 - 1)) - d[1:M]) .* g[1:M];
				exit = 1;

			} else if (M < nalts + 1)
				M += 1; // adds one to M
		}
	} else { //General approach to demand simulation (het. alpha's)
		real lambda_l;
		real lambda_u;
		matrix[nalts+1, 4] parm_matrix;
		vector[nalts+1] psi_ord;
		vector[nalts+1] mu;
		vector[nalts+1] price_ord; // price
		vector[nalts+1] a;//	alpha
		vector[nalts+1] g_psi_a; // gamma*psi/alpha
		vector[nalts+1] a_a_1; // (alpha/(alpha-1))
		vector[nalts+1] mu_a_a_1; // (1/MUzero)^(alpha/(alpha-1))
		vector[nalts+1] g__phi; // obtain gamma
		real alpha_1;

		if (model_num < 5){
			parm_matrix = SortParmMatrix(MUzero, price, gamma, alpha, nalts);
			mu = col(parm_matrix, 1); // obtain mu
			g__phi = col(parm_matrix, 3); // obtain gamma
			price_ord = col(parm_matrix, 2); // price
			a = col(parm_matrix, 4);//	alpha
			a_a_1 = a ./ (a - 1); // (alpha/(alpha-1))
			psi_ord = mu .* price_ord;
			g_psi_a = g__phi .* psi_ord ./ a; // gamma*psi/alpha
			mu_a_a_1 = pow(inv(mu), a_a_1);

		} else {
			alpha_1 = alpha[1];
			a = alpha;
			parm_matrix = SortParmMatrix(MUzero, price, gamma, phi, nalts);
			mu = col(parm_matrix, 1); // obtain mu = psi*phi/(gamma*price)
			g__phi = col(parm_matrix, 3) ./ col(parm_matrix, 4); // gamma/phi
			psi_ord = mu .* col(parm_matrix, 2) .* g__phi;
		}

		while (exit == 0){
			lambda1 = mu[M + 1];// Calculate lambda1 equal to MUzero(M+1)

			// Calculate new utility
			if (model_num < 5)
				util_new = ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi_ord, g__phi, price_ord, d, model_num);
			else
				util_new = ComputeKtUtilM(M, lambda1, psi_ord, mu, alpha_1);

			if (util_new >= util || M+1 == nalts+1){
				if(util_new < util)
					M += 1;
				lambda_l = util_new < util ? 0 : lambda1;
				lambda_u = mu[M];
				lambda1 = (lambda_l + lambda_u) / 2;

				for (n in 1:max_loop){
					if (model_num < 5)
						util_new = ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi_ord, g__phi, price_ord, d, model_num);
					else
						util_new = ComputeKtUtilM(M, lambda1, psi_ord, mu, alpha_1);

					// Update lambda bounds
					if (util_new < util)
						lambda_u = lambda1;
					else if (util_new > util)
						lambda_l = lambda1;

					lambda1 = (lambda_l + lambda_u) / 2;

					if (abs((lambda_l - lambda_u) / (lambda_l + lambda_u) * 0.5) < tol_l) break;
				}

			// Compute demands (using eq. 12 in Pinjari and Bhat)
			X[1:M] = (pow(lambda1 ./ mu[1:M], 1 ./ (a[1:M]-1)) - d[1:M]) .* g__phi[1:M];

			exit = 1;

			} else if (util_new < util && M+1 < nalts+1)
				M += 1; // adds one to M
		}
	}
	// This code puts the choices back in their original order and exports demands
	hdemand = X[order_x];

return hdemand;
}

/**
 * Calculate WTP for each individual, simulation, and policy
 * @param income            individual income
 * @param quant_j           non-numeraire consumption amounts (length nalts)
 * @param price             full baseline price vector (length nalts+1)
 * @param price_p_policy    array[npols] of price-change vectors (length nalts+1)
 * @param psi_p_sims        array[nsims] of npols x nalts policy psi matrices
 * @param phi_p_sims        array[nsims] of npols x nalts policy phi matrices
 * @param psi_sims          matrix nsims x nalts of baseline psi
 * @param phi_sims          matrix nsims x nalts of baseline phi
 * @param gamma_sims        matrix nsims x nalts of baseline gamma
 * @param alpha_sims        matrix nsims x (nalts+1) of baseline alpha
 * @param scale_sims        vector nsims of scale parameters
 * @param nerrs             number of error draws per simulation
 * @param cond_error        0 = unconditional, 1 = conditional on observed quantities
 * @param draw_mlhs         0 = iid uniform, 1 = MLHS
 * @param algo_gen          0 = hybrid, 1 = general bisection
 * @param model_num         model type (1-5)
 * @param price_change_only 1 = price change only (psi unchanged across policies)
 * @param tol               convergence tolerance
 * @param max_loop          maximum bisection iterations
 * @return matrix nsims x npols of WTP values
 */
matrix CalcWTP_rng(real income, vector quant_j, vector price,
					array[] vector price_p_policy, array[] matrix psi_p_sims, array[] matrix phi_p_sims,
					matrix psi_sims, matrix phi_sims, matrix gamma_sims, matrix alpha_sims, vector scale_sims,
					int nerrs, int cond_error, int draw_mlhs,
					int algo_gen, int model_num, int price_change_only, real tol, int max_loop){

	int nalts = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	matrix[nsims, npols] wtp;
	real quant_num = income - quant_j' * price[2:(nalts+1)];

	for (sim in 1:nsims){
		vector[nalts] psi_j = psi_sims[sim]';
		vector[nalts + 1] phi;
		matrix[npols, nalts] psi_p_policy;
		vector[nalts + 1] gamma = append_row(1, gamma_sims[sim]');
		vector[nalts + 1] alpha = alpha_sims[sim]';
		real scale = scale_sims[sim];
		array[nerrs] vector[nalts + 1] error;
		vector[npols] wtp_policy;
		vector[nerrs] util;

		if (price_change_only == 0)
			psi_p_policy = psi_p_sims[sim];

		if (model_num < 5)
			phi = rep_vector(1, nalts + 1);
		else
			phi = append_row(1,phi_sims[sim]');

		error = DrawError_rng(quant_num, quant_j, price[2:(nalts+1)],
						psi_j, phi[2:(nalts+1)], gamma[2:(nalts+1)], alpha, scale,
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

			util[err] = ComputeUtilJ(income, mdemand[2:(nalts+1)], price[2:(nalts+1)],
								psi_b_err, phi[2:(nalts+1)], gamma[2:(nalts+1)], alpha,
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
return wtp;
}

/**
 * Calculate Marshallian demand for each individual, simulation, and policy
 * @param income            individual income
 * @param quant_j           non-numeraire consumption amounts (length nalts)
 * @param price             full baseline price vector (length nalts+1)
 * @param price_p_policy    array[npols] of price-change vectors (length nalts+1)
 * @param psi_p_sims        array[nsims] of npols x nalts policy psi matrices
 * @param phi_p_sims        array[nsims] of npols x nalts policy phi matrices
 * @param psi_sims          matrix nsims x nalts of baseline psi
 * @param phi_sims          matrix nsims x nalts of baseline phi
 * @param gamma_sims        matrix nsims x nalts of baseline gamma
 * @param alpha_sims        matrix nsims x (nalts+1) of baseline alpha
 * @param scale_sims        vector nsims of scale parameters
 * @param nerrs             number of error draws per simulation
 * @param cond_error        0 = unconditional, 1 = conditional on observed quantities
 * @param draw_mlhs         0 = iid uniform, 1 = MLHS
 * @param algo_gen          0 = hybrid, 1 = general bisection
 * @param model_num         model type (1-5)
 * @param price_change_only 1 = price change only (psi unchanged across policies)
 * @param tol               convergence tolerance
 * @param max_loop          maximum bisection iterations
 * @return array[nsims] of matrix[npols, nalts+1] demand values
 */
array[] matrix CalcMarshallianDemand_rng(real income, vector quant_j, vector price,
					array[] vector price_p_policy, array[] matrix psi_p_sims, array[] matrix phi_p_sims,
					matrix psi_sims, matrix phi_sims, matrix gamma_sims, matrix alpha_sims, vector scale_sims,
					int nerrs, int cond_error, int draw_mlhs,
					int algo_gen, int model_num, int price_change_only, real tol, int max_loop){

	int nalts = num_elements(quant_j);
	int nsims = num_elements(scale_sims);
	int npols = size(price_p_policy);
	array[nsims] matrix[npols, nalts+1] mdemand_out;
	real quant_num = income - quant_j' * price[2:(nalts+1)];

	for (sim in 1:nsims){
		vector[nalts] psi_j = psi_sims[sim]';
		vector[nalts + 1] phi;
		matrix[npols, nalts] psi_p_policy;
		vector[nalts + 1] gamma = append_row(1, gamma_sims[sim]');
		vector[nalts + 1] alpha = alpha_sims[sim]';
		real scale = scale_sims[sim];
		array[nerrs] vector[nalts + 1] error;
		matrix[npols, nalts + 1] mdemand_pols;
		vector[nerrs] util;

		if (price_change_only == 0)
			psi_p_policy = psi_p_sims[sim];

		if (model_num < 5)
			phi = rep_vector(1, nalts + 1);
		else
			phi = append_row(1,phi_sims[sim]');

		error = DrawError_rng(quant_num, quant_j, price[2:(nalts+1)],
						psi_j, phi[2:(nalts+1)], gamma[2:(nalts+1)], alpha, scale,
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

			util[err] = ComputeUtilJ(income, mdemand_util[2:(nalts+1)], price[2:(nalts+1)],
								psi_b_err, phi[2:(nalts+1)], gamma[2:(nalts+1)], alpha,
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
return mdemand_out;
}

/**
 * Calculate baseline Marshallian demand with a single simulation and no policies.
 * Used to simulate data.
 * Note: error draws here are iid uniform (not MLHS), intentionally.
 * @param income    individual income
 * @param price     full price vector (length nalts+1)
 * @param psi_j     non-numeraire psi parameters (length nalts)
 * @param phi_j     non-numeraire phi parameters (kt_ee only; length nalts)
 * @param gamma_j   non-numeraire gamma parameters (length nalts)
 * @param alpha     full alpha vector (length nalts+1)
 * @param scale     scale parameter
 * @param nerrs     number of error draws
 * @param model_num model type (1-5)
 * @param algo_gen  0 = hybrid, 1 = general bisection
 * @param tol       convergence tolerance
 * @param max_loop  maximum bisection iterations
 * @return vector of nalts+1 average demands
 */
vector CalcmdemandOne_rng(real income, vector price,
					vector psi_j, vector phi_j, vector gamma_j, vector alpha, real scale,
					int nerrs, int model_num, int algo_gen, real tol, int max_loop){

	int nalts = num_elements(price) - 1; // subtract numeraire
	vector[nalts + 1] mdemand = rep_vector(0, nalts + 1);
	vector[nalts + 1] gamma = append_row(1, gamma_j);
	vector[nalts + 1] phi;
	array[nerrs] vector[nalts + 1] error;

	if (model_num < 5)
		phi = rep_vector(1, nalts + 1);
	else
		phi = append_row(1, phi_j);

	for(err in 1:nerrs)
		for(g in 1:(nalts+1))
			error[err, g] = -log(-log(uniform_rng(0, 1))) * scale; //uniform(0,1) draws

	// Compute Marshallian demands and baseline utility
	for (err in 1:nerrs){
		vector[nalts + 1] MUzero_b = exp(append_row(0, psi_j) + error[err]) ./ price;
		if (model_num == 5)
			MUzero_b = MUzero_b .* phi ./ gamma; //change for no psi_p

		mdemand = mdemand + MarshallianDemand(income, price, MUzero_b, phi, gamma, alpha,
					nalts, algo_gen, model_num, tol, max_loop) / nerrs; // take average
	}

return mdemand;
}

// end of functions
}
