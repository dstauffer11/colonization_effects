data {
  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> M;
  int<lower=1> N;
  int<lower=1> J;
  matrix[T, M*J] X;
  vector[T] y[J];
  int<lower=1, upper=M*J> indices[2*J];
  int<lower=1, upper=N> indices_short[2*J];


//  int<lower=1> C;
//  int<lower=0, upper=1> colonizer[J, T];
}

transformed data {
  int T1;
  matrix[T, M*J] X_centered;
  vector[T] y_centered[J];
  int n;
  matrix[T, M*J] Q = rep_matrix(0, T, M*J);
  matrix[M, M*J] R = rep_matrix(0, M, M*J);
  matrix[M, M*J] R_inv = rep_matrix(0, M, M*J);
  

  T1 = T0 + 1;

  for (mj in 1:M*J) {
    X_centered[:, mj] = X[:, mj] - mean(X[:, mj]);
  }



  for (j in 1:J) {
    y_centered[j] = y[j] - mean(y[j]);

    n = indices[2*j] - indices[2*j-1] + 1;
    // for each case, do QR decomposition and store in wide Q, R, R_inv matrices
    Q[:, indices[2*j-1]:indices[2*j]] = qr_Q(X_centered[:, indices[2*j-1]:indices[2*j]])[:, 1:n] * T;
    R[1:n, indices[2*j-1]:indices[2*j]] = qr_R(X_centered[:, indices[2*j-1]:indices[2*j]])[1:n, :] / T;
    R_inv[1:n, indices[2*j-1]:indices[2*j]] = inverse(R[1:n, indices[2*j-1]:indices[2*j]]);
  }

}


parameters {
  vector[T] u_err[J]; //Slope innovation
  vector[T] v_err[J]; //Level innovation
  vector[N] beta_tilde;
  real<lower=0> s_obs[J];
  real<lower=0> s_slope[J];
  real<lower=0> s_level[J];

  real m_effect;
  real effect[J];

  real<lower=0> s_n_beta;
  real<lower=0> s_u_err;
  real<lower=0> s_v_err;
  real<lower=0> s_effect;
  real<lower=0> s_s_obs;
//  real<lower=0> s_s_slope;
//  real<lower=0> s_s_level;
}

transformed parameters {
  vector[T] u[J]; //Level
  vector[T] v[J]; //Slope
  vector[N] neighbor_beta; // country regressio coefficients

  
  // state space
  for (j in 1:J) {
    u[j][1] = u_err[j][1];
    v[j][1] = v_err[j][1];
    for (t in 2:T) {
      u[j][t] = u[j][t-1] + v[j][t-1] + s_level[j] * u_err[j][t];
      v[j][t] = v[j][t-1] + s_slope[j] * v_err[j][t];
    }
  }

  for (j in 1:J) {
    // for each case, neighbor_beta = R_inv*beta_tilde
    neighbor_beta[indices_short[2*j-1]:indices_short[2*j]] = R_inv[1:indices[2*j] - indices[2*j-1] + 1, indices[2*j-1]:indices[2*j]] * beta_tilde[indices_short[2*j-1]:indices_short[2*j]];
  }
  
}

model {
  m_effect ~ std_normal();
  s_n_beta ~ std_normal();
  s_u_err ~ std_normal();
  s_v_err ~ std_normal();
  s_effect ~ cauchy(0, 5);
  s_s_obs ~ cauchy(0, 5);
//  s_s_slope ~ normal(0, 0.03);
//  s_s_level ~ normal(0, 0.03);


  effect ~ normal(m_effect, s_effect);
  neighbor_beta ~ normal(0, s_n_beta);

  // state space priors
  for (j in 1:J) {
    u_err[j] ~ normal(0, s_u_err);
    v_err[j] ~ normal(0, s_v_err);
    s_obs[j] ~ cauchy(0, s_s_obs);
    s_slope[j] ~ normal(0, 0.03); // s_s_slope
    s_level[j] ~ normal(0, 0.03); // s_s_level

    // fitting state-space model to pre-intervention data
    y_centered[j][1:T0] ~ normal(u[j][1:T0] + Q[1:T0, indices[2*j-1]:indices[2*j]]*beta_tilde[indices_short[2*j-1]:indices_short[2*j]], s_obs[j]);

    // fitting effect to make state-space prediction match observed data
    y_centered[j][T1:T] ~ normal(u[j][T1:T] + Q[T1:T, indices[2*j-1]:indices[2*j]]*beta_tilde[indices_short[2*j-1]:indices_short[2*j]] + effect[j], s_obs[j]);
  }
}


