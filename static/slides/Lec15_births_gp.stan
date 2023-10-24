// Based on gpbf2.stan from https://avehtari.github.io/casestudies/Birthdays/birthdays.html

functions {
  vector diagSPD_EQ(real sigma, real ls, real L, int M) {
    return sigma * sqrt(sqrt(2*pi()) * ls) * exp(-0.25*(ls*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }
  vector diagSPD_periodic(real sigma, real ls, int M) {
    real a = 1/ls^2;
    vector[M] q = exp(log(sigma) + 0.5 * (log(2) - a + to_vector(log_modified_bessel_first_kind(linspaced_int_array(M, 1, M), a))));
    return append_row(q,q);
  }
  matrix PHI(int N, int M, real L, vector x) {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }
  matrix PHI_periodic(int N, int M, real w0, vector x) {
    matrix[N,M] mw0x = diag_post_multiply(rep_matrix(w0*x, M), linspaced_vector(M, 1, M));
    return append_col(cos(mw0x), sin(mw0x));
  }
}

data {
  int<lower=1> N;      // number of observations
  vector[N] x;         // univariate covariate
  vector[N] y;         // target variable
  int<lower=1> M_f1;   // number of basis functions for smooth function
  int<lower=1> M_f2;   // number of cos and sin functions for periodic
  real period;
}

transformed data {
  // Basis functions for f1
  real L_f1 = 1.5 * max(x);
  matrix[N,M_f1] PHI_f1 = PHI(N, M_f1, L_f1, x);
  
  // Basis functions for f2
  //real period = 7;
  matrix[N,2*M_f2] PHI_f2 = PHI_periodic(N, M_f2, 2*pi()/period, x);
  
  // Concatenated basis functions
  matrix[N,M_f1+2*M_f2] PHI_f = append_col(PHI_f1, PHI_f2);
}

parameters {
  vector[M_f1] beta_f1;         // the basis functions coefficients for f1
  vector[2*M_f2] beta_f2;       // the basis functions coefficients for f2
  real<lower=0> lengthscale_f1;
  real<lower=0> lengthscale_f2;
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma_f2;       // scale of f2
  real<lower=0> sigma;          // residual scale
}

model {
  // spectral densities for f1 and f2
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  vector[2*M_f2] diagSPD_f2 = diagSPD_periodic(sigma_f2, lengthscale_f2, M_f2);
  
  // priors
  beta_f1 ~ normal(0, 1);
  beta_f2 ~ normal(0, 1);
  lengthscale_f1 ~ normal(0, 1); //lognormal(log(700/xsd), 1);
  lengthscale_f2 ~ normal(0, 1);
  sigma_f1 ~ normal(0, 1);
  sigma_f2 ~ normal(0, 1);
  sigma ~ normal(0, .5);
  
  // model
  y ~ normal_id_glm(
    PHI_f,
    0.0,
    append_row(diagSPD_f1 .* beta_f1, diagSPD_f2 .* beta_f2),
    sigma
  ); 
}

generated quantities {
  vector[N] f1;
  vector[N] f2;
  {
    // spectral densities for f1
    vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
    vector[2*M_f2] diagSPD_f2 = diagSPD_periodic(sigma_f2, lengthscale_f2, M_f2);
    // functions scaled back to original scale
    f1 = PHI_f1 * (diagSPD_f1 .* beta_f1);
    f2 = PHI_f2 * (diagSPD_f2 .* beta_f2);
  }
}
