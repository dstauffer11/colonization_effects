data {

  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> N;
  int<lower=1> J;
  matrix[T, N] X;
  vector[T] y[J];
  int indices[2*J];

}

transformed data {

  int T1;
  matrix[T, N] X_norm;
  vector[J] norms;
  vector[T] y_norm[J];
  vector[J] post_growth;

  T1 = T0 + 1;

  X_norm = 10*X;

  for (j in 1:J) {
    y_norm[j] = 10*y[j];
  }
  
  for (j in 1:J) {
    post_growth[j] = sum(y[j][T1:T]);
  }
}


parameters {
  vector[T] u_err[J]; //Slope innovation
  vector[T] v_err[J]; //Level innovation
  vector[N] neighbor_beta;
  real<lower=0> s_obs[J];
  real<lower=0> s_slope[J];
  real<lower=0> s_level[J];
  real<lower=0> s_s_obs;

  real mean_effect;
  // real<lower=0, upper=pi()/2> s_effect_unif;
  real<lower=0> s_effect;
  vector[J] raw_effects;

}

transformed parameters {
  vector[T] u[J]; //Level
  vector[T] v[J]; //Slope
  // real<lower=0> s_effect;
  vector[J] real_effects;
  
  
  // state space
  for (j in 1:J) {
    u[j][1] = u_err[j][1];
    v[j][1] = v_err[j][1];
    for (t in 2:T) {
      u[j][t] = u[j][t-1] + v[j][t-1] + s_level[j] * u_err[j][t];
      v[j][t] = v[j][t-1] + s_slope[j] * v_err[j][t];
    }
  }

  real_effects = (mean_effect + s_effect * raw_effects);

}

model {
  mean_effect ~ normal(0, 0.2);
  s_effect ~ normal(0, 0.2);
  raw_effects ~ std_normal();

  s_s_obs ~ cauchy(0, 1);
  s_obs ~ normal(0, s_s_obs);

  neighbor_beta ~ std_normal();

  for (j in 1:J) {
    // state space priors
    u_err[j] ~ std_normal();
    v_err[j] ~ std_normal();

    s_slope[j] ~ normal(0, 0.03);
    s_level[j] ~ normal(0, 0.03);

    // fitting state-space model to pre-intervention data
    y_norm[j][1:T0] ~ normal(u[j][1:T0] + X_norm[1:T0, indices[2*j-1]:indices[2*j]]*neighbor_beta[indices[2*j-1]:indices[2*j]], s_obs[j]);


    post_growth[j] ~ normal(sum(u[j][T1:T] + X_norm[T1:T, indices[2*j-1]:indices[2*j]]*neighbor_beta[indices[2*j-1]:indices[2*j]]) / 10 + real_effects[j], sqrt(T-T0)*s_obs[j]);
  }

}


