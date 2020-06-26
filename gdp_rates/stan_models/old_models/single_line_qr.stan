data {
  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> N;
  matrix[T, N] X;
  vector[T] y;
}

transformed data {
  int T1;
  matrix[T, N] X_centered;
  vector[T] y_centered;
  matrix[T, N] Q;
  matrix[N, N] R;
  matrix[N, N] R_inv;
  

  T1 = T0 + 1;

  for (n in 1:N) {
    X_centered[:, n] = X[:, n] - mean(X[:, n]);
  }
  y_centered = y - mean(y);

  Q = qr_Q(X_centered)[:, 1:N] * T;
  R = qr_R(X_centered)[1:N, :] / T;
  R_inv = inverse(R);
}

parameters {
  vector[T] u_err; //Slope innovation
  vector[T] v_err; //Level innovation
  vector[N] beta_tilde;
  real<lower=0> s_obs;
  real<lower=0> s_slope;
  real<lower=0> s_level;

  real effect;
}

transformed parameters {
  vector[T] u; //Level
  vector[T] v; //Slope
  vector[N] beta;

  
  // state space
  u[1] = u_err[1];
  v[1] = v_err[1];
  for (t in 2:T) {
    u[t] = u[t-1] + v[t-1] + s_level * u_err[t];
    v[t] = v[t-1] + s_slope * v_err[t];
  }

  beta = R_inv * beta_tilde;

}

model {
  // state space priors
  u_err ~ normal(0, 1);
  v_err ~ normal(0, 1);
  beta ~ normal(0, 1);
  s_obs ~ normal(0, 5);
  s_slope ~ normal(0, 0.1);
  s_level ~ normal(0, 0.1);

  effect ~ normal(0, 1);

  // fitting state-space model to pre-intervention data
  y[1:T0] ~ normal(u[1:T0] + Q[1:T0, :]*beta_tilde, s_obs);

  // fitting effect to make state-space prediction match observed data
  y[T0+1:T] ~ normal(u[T0+1:T] + Q[T0+1:T, :]*beta_tilde + effect, s_obs);

}


