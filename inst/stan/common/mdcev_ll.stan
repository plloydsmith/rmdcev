
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

matrix alpha_ll(vector alpha, int I, int J, int model_num) {

	matrix[I, J+1] alpha_full;

	if (model_num == 1)
	  alpha_full = append_col(rep_vector(alpha[1], I), rep_matrix(0, I, J));
	else if (model_num == 2)
	  alpha_full = rep_matrix(alpha', I);
	else if (model_num == 3)
	  alpha_full = rep_matrix(alpha[1], I, J+1);
	else if (model_num == 4)
	  alpha_full = rep_matrix(1e-03, I, J+1);

return(alpha_full);
}

vector mdcev_ll(matrix quant_j, matrix price_j, vector log_num, vector income,
				vector M, vector log_M_fact, // data
				matrix lpsi, matrix gamma_j, vector alpha1, matrix alpha_j, real scale_full,// parameters
				int I, int J, matrix nonzero, int trunc_data)  {//options

	vector[I] log_like;
	vector[J] ones_j = rep_vector(1, J);
	matrix[I, J] v_j= lpsi + (alpha_j - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
	vector[I] v1 = (alpha1 - 1) .* log_num / scale_full;
	matrix[I, J] logf = log(1 - alpha_j) - log(quant_j + gamma_j);
	vector[I] logf1 = log(1 - alpha1) - log_num;
	v_j = v_j / scale_full;

	if (trunc_data == 0){
		log_like = (1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* price_j ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact;

	} else if (trunc_data == 1){
		matrix[I, J+1] v_1;
		vector[I] like_cond;
		vector[I] like_trunc;
		vector[I] sumv;

		like_cond = exp((1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* price_j ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact);

		v_1 = append_col((alpha1 - 1) .* log(income), lpsi - log(price_j));
		v_1 = exp(v_1 / scale_full);
		sumv = v_1 * rep_vector(1, J+1);

		like_trunc = col(v_1, 1) ./ sumv;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log(like_cond ./ (1 - like_trunc));
	}

return(log_like);
}

real DeterminJacob(vector phi_quant_gamma, real alpha, vector phi_j, vector price_j_num,
                  vector nonzero, int J){

  matrix[J, J] jacobian;
  real j_det;

  jacobian = rep_matrix((1 - alpha) * price_j_num, J);
  jacobian = jacobian + diag_matrix(phi_j  ./ phi_quant_gamma);
//    jacobian = add_diag(jacobian, phi_j  ./ phi_quant_gamma); // Add row_vector d to the diagonal of matrix m.
  jacobian = diag_post_multiply(jacobian, nonzero) + diag_matrix(1 - nonzero);
  j_det = fabs(determinant(jacobian));

return(j_det);

}

vector kt_ll(vector income, vector log_num, matrix quant_j, matrix price_j,
          matrix psi_i, matrix phi_ij, matrix gamma, vector alpha, real scale_full,
          int I, int J, matrix nonzero, int trunc_data){

    vector[I] log_like;
    vector[J] ones_j = rep_vector(1, J);
    matrix[I, J] g; // demand function
    vector[I] like;
    matrix[I, J] phi_quant_gamma =  phi_ij .* quant_j + gamma;
    vector[I] log_j_det;

//	  for(i in 1:I){
//		vector[J] price_j_num = price_j[i]' ./ exp(log_num[i]);
//         j_det[i] = DeterminJacob(phi_quant_gamma[i]', alpha[i], phi_ij[i]',
//                    price_j_num, nonzero[i]', J);
//    }
	log_j_det = log(1 - alpha) - log_num + (nonzero .* (log(phi_ij) - log(phi_quant_gamma))) * ones_j +
         		log(exp(log_num) ./ (1 - alpha) + (nonzero .* phi_quant_gamma .* price_j ./ phi_ij) * ones_j);

	 // Calculate the demand function, g
  	g =  (-psi_i  + log(price_j) - log(phi_ij) +
  		log(phi_quant_gamma) - rep_matrix((1 - alpha) .* log_num, J)) ./ scale_full;

  	// Calculate the liklihood
  	like = (nonzero .*(-g - log(scale_full)) + (-exp(-g))) * ones_j;

	// adjust for truncation
	if(trunc_data == 1){
	  matrix[I, J] g_t = (-psi_i  + log(price_j) - log(phi_ij) + log(gamma) -
	                      rep_matrix((1 - alpha) .* log(income), J)) ./ scale_full;
  	  vector[I] like_trunc = exp(-exp(-g_t) * ones_j);

      for(i in 1:I)
		  like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

   	      log_like = log_j_det + like - log(1 - like_trunc);
	}
	else
  		log_like = log_j_det + like ; // Calculate the liklihood log(j_det*exp(x))=log(j_det)+x

return(log_like);
}
