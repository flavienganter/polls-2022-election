data {
  
  // Poll data
  int N;
  int tot_eff[N];
  int vote_eff[N];
  int id_cand;
  int id_poll[N];
  int P;
  int id_date[N];
  int id_month[N];
  int M;
  int id_house[N];
  int F;
  matrix[N,3] X;
  
  // EZ and CT adjustement
  int isn_z[N];
  int isn_t[N];
  
  // Splines
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_basis;
  int D;
  matrix[num_basis,D] B;
  
}
parameters {
  
  // Splines
  real alpha0;
  real alpha_raw[num_basis];
  real<lower=0> tau_alpha;
  
  // Covariates
  real mu[P];
  real lambda[F];
  real tau_lambda_tilde;
  real beta[3];
  real nu[2];
  real<lower=0> tau_mu;
  real<lower=0> tau_lambda;
  
  // EZ and CT adjustment
  real gamma_z_tilde;
  real gamma_t_tilde;
  
}
transformed parameters {
  
  real alpha[num_basis];
  real gamma_z;
  real gamma_t;
  
  // Spline coefficients, specified as a random walk
  // to avoid overfit
  alpha[1] = alpha_raw[1];
  for (i in 2:num_basis)
    alpha[i] = alpha[i-1] + alpha_raw[i] * tau_alpha;
  
  // EZ and CT adjustment
  if (id_cand < 11) {
    gamma_z = gamma_z_tilde;
    gamma_t = gamma_t_tilde;
  } else if (id_cand == 11) {
    gamma_z = 0;
    gamma_t = gamma_t_tilde;
  } else {
    gamma_z = 0;
    gamma_t = 0;
  }
  
}
model {
  
  
  // Priors
  
    // Splines
    alpha0 ~ normal(0, 1);
    alpha_raw ~ normal(0, 2);
    tau_alpha ~ student_t(3, 0, 1);
    
    // Poll effect
    mu ~ normal(0, 1);
    tau_mu ~ student_t(3, 0, 1);
    
    // House effect
    lambda ~ normal(0, 1);
    tau_lambda ~ student_t(3, 0, 1);
    
    // Population definition and poll type effects
    beta ~ normal(0, 1);
    nu ~ normal(0, 1);
    
    // EZ and CT adjustment
    gamma_z_tilde ~ normal(0, 1);
    gamma_t_tilde ~ normal(0, 1);
  
  
  // Likelihood
  for (i in 1:N) {
    if (id_month[i] == 1) {
      target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  alpha0 * id_date[i] + to_row_vector(alpha) * B[,id_date[i]] + // Spline
                                  tau_mu * mu[id_poll[i]] + // Poll effect
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3] +
                                  gamma_z * isn_z[i]); // EZ adjustment
    } else if (id_month[i] < 4 || id_month[i] > 5) {
      target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  alpha0 * id_date[i] + to_row_vector(alpha) * B[,id_date[i]] + // Spline
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3]);
    } else {
      target += binomial_logit_lpmf(vote_eff[i] | tot_eff[i],
                                  alpha0 * id_date[i] + to_row_vector(alpha) * B[,id_date[i]] + // Spline
                                  tau_mu * mu[id_poll[i]] + // Poll effect
                                  tau_lambda * lambda[id_house[i]] + // House effect
                                  X[i,1] * (beta[1] + nu[1] * (id_date[i] - 1)) + // Population definition and poll type effects
                                  X[i,2] * (beta[2] + nu[2] * (id_date[i] - 1)) +
                                  X[i,3] * beta[3] +
                                  gamma_t * isn_t[i]); // CT adjustment
    }
  }
                                       
}
generated quantities {
  
  // Estimate each candidate's score for every day between Sep 1
  // and the date of the most recent poll
  real<lower=0,upper=1> prob[D];
  for (d in 1:D)
      prob[d] = inv_logit(alpha0 * d + to_row_vector(alpha) * B[,d]);
      
}
