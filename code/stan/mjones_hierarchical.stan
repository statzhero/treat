//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data{
  int<lower=1> N;                  // num obs
  int<lower=1> J;                  // num industry-years
  int<lower=1> K;                  // num coefficients
  int<lower=1, upper=J> IndYearID[N];  // Industry-year for obs
  vector[N] TA;                    // Outcome total accruals
  matrix[N, K] x;                  // Predictors InAt, ChRev, PPE
}
parameters{
  matrix[K, J] z;                  // standard normal sampler
  cholesky_factor_corr[K] L_Omega; // hypprior coefficient correlation
  vector<lower=0>[K] tau;          // hypprior coefficient scales
  vector[K] mu_b;                  // hypprior mean coefficients
  real<lower=0> sigma;             // error-term scale
}
transformed parameters{
  matrix[J, K] b;                  // coefficient vector
  // The multivariate non-centered version:
  b = (rep_matrix(mu_b, J) + diag_pre_multiply(tau,L_Omega) * z)';
}
model{
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  mu_b  ~ normal(0, 2.5);
  sigma ~ exponential(1);
  tau ~ exponential(1);
  TA ~ normal(rows_dot_product(b[IndYearID] , x), sigma);
}
generated quantities {
  vector[N] y_fit_wo;
  vector[N] y_fit;
  vector[N] log_lik;
  for ( i in 1:N ) {
    y_fit[i] = normal_rng(b[IndYearID[i], 1] * x[i,1] +
                          b[IndYearID[i], 2] * x[i,2] +
                          b[IndYearID[i], 3] * x[i,3], sigma);
    y_fit_wo[i] = b[IndYearID[i], 1] * x[i,1] +
                  b[IndYearID[i], 2] * x[i,2] +
                  b[IndYearID[i], 3] * x[i,3];
    log_lik[i] = normal_lpdf(TA[i] | b[IndYearID[i], 1] * x[i,1] +
                                     b[IndYearID[i], 2] * x[i,2] +
                                     b[IndYearID[i], 3] * x[i,3], sigma);
  }
}