
matrix gamma_ll(vector gamma, int I, int J, int Gamma) {

	matrix[I, J] gamma_j;

	if (Gamma == 0)
	  gamma_j = rep_matrix(1, I, J);
	else if (Gamma == J)
	  gamma_j = rep_matrix(gamma', I);
	else if (Gamma == 1)
	  gamma_j = rep_matrix(gamma[1], I, J);

return(gamma_j);
}

matrix alpha_j_ll(vector alpha, int I, int J, int model_num) {

	matrix[I, J] alpha_j;

	if (model_num == 1 || model_num == 5)
	  alpha_j = rep_matrix(0, I, J);
	else if (model_num == 2)
	  alpha_j = rep_matrix(alpha[2:(J+1)]', I);
	else if (model_num == 3)
	  alpha_j = rep_matrix(alpha[1], I, J);
	else if (model_num == 4)
	  alpha_j = rep_matrix(1e-03, I, J);

return(alpha_j);
}

vector alpha_1_ll(vector alpha, int I, int model_num) {

	vector[I] alpha_1;

	if (model_num == 4)
	  alpha_1 = rep_vector(1e-03, I);
	else
	  alpha_1 = rep_vector(alpha[1], I);

return(alpha_1);
}

vector mdcev_ll(matrix quant_j, matrix price_j, vector log_num, vector income,
				vector M, vector log_M_fact, // data
				matrix lpsi, matrix gamma_j, vector alpha1, matrix alpha_j, real scale_full,// parameters
				int I, int J, matrix nonzero, int trunc_data)  {//options

	vector[I] log_like;
	vector[J] ones_j = rep_vector(1, J);
	matrix[I, J] v_j= (lpsi + (alpha_j - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j)) / scale_full;
	vector[I] v1 = (alpha1 - 1) .* log_num / scale_full;
	matrix[I, J] logf = log1m(alpha_j) - log(quant_j + gamma_j);
	vector[I] logf1 = log1m(alpha1) - log_num;

	log_like = (1 - M) * log(scale_full) + logf1 + v1 + rows_dot_product(nonzero, logf + v_j) +
		log(inv(exp(logf1)) + rows_dot_product(nonzero, price_j ./ exp(logf))) -
		M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact;

	if (trunc_data == 1){
		matrix[I, J+1] v_1 = exp(append_col((alpha1 - 1) .* log(income), lpsi - log(price_j)) / scale_full);
		vector[I] like_trunc = v_1[ , 1] ./ (v_1 * rep_vector(1, J+1));

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log_like - log1m(like_trunc);
	}

return(log_like);
}

real DeterminJacob(vector phi_quant_term, real alpha, vector price_j_num,
                  vector nonzero, int J){
// older version, use Bhat (2008) form below
  real log_j_det;
  matrix[J, J] jacobian = add_diag(rep_matrix((1 - alpha) * price_j_num, J), inv(phi_quant_term));
  jacobian = add_diag(diag_post_multiply(jacobian, nonzero), 1 - nonzero);
  log_j_det = log_determinant(jacobian);

return(log_j_det);
}

vector kt_ll(matrix quant_j, matrix price_j, vector log_num, vector income,
          matrix lpsi, matrix phi_ij, matrix gamma, vector alpha, real scale_full,
          int I, int J, matrix nonzero, int trunc_data, int jacobian_analytical_grad){

    vector[I] log_like;
    vector[J] ones_j = rep_vector(1, J);
    matrix[I, J] g; // demand function
    matrix[I, J] phi_quant_term =  (phi_ij .* quant_j + gamma) ./ phi_ij;
    vector[I] log_j_det;

	if (jacobian_analytical_grad == 0){
	  for(i in 1:I){
		vector[J] price_j_num = price_j[i]' ./ exp(log_num[i]);
         log_j_det[i] = DeterminJacob(phi_quant_term[i]', alpha[i],
                    price_j_num, nonzero[i]', J);
    	}
	} else if (jacobian_analytical_grad == 1){
		log_j_det = log1m(alpha) - log_num + rows_dot_product(nonzero, log(inv(phi_quant_term))) +
         		log(exp(log_num) ./ (1 - alpha) + rows_dot_product(nonzero, phi_quant_term .* price_j));
	}
	 // Calculate the demand function, g
  	g =  (-lpsi + log(phi_quant_term .* price_j) - rep_matrix((1 - alpha) .* log_num, J)) / scale_full;

  	// Calculate the liklihood log(j_det*exp(x))=log(j_det)+x
  	log_like = log_j_det + (nonzero .*(-g - log(scale_full)) + (-exp(-g))) * ones_j ;

	// adjust for truncation
	if(trunc_data == 1){
	  vector[I] like_trunc = exp(-exp(-(-lpsi  + log(price_j) - log(phi_ij) + log(gamma) -
	                      rep_matrix((1 - alpha) .* log(income), J)) / scale_full) * ones_j);
      for(i in 1:I)
		  like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

   	      log_like = log_like - log1m(like_trunc);
	}

return(log_like);
}
