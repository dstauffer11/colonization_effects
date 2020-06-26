data {
  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> N;
  int<lower=1> J;
  matrix[T, N] X;
  vector[T] y[J];
  int indices[2*J];


  int<lower=1> C;
  int<lower=1, upper=C> colonizer[J];
}

transformed data {
  int T1;
  matrix[T, N] X_norm;
  vector[J] norms;
  vector[T] y_norm[J];

  T1 = T0 + 1;

  for (i in 1:N) {
    X_norm[:, i] = X[:, i] / sd(X[:, i]);
  }

  for (j in 1:J) {
    norms[j] = sd(y[j]);
    y_norm[j] = y[j] / norms[j];
  }
}


parameters {
  vector[T] u_err[J]; //Slope innovation
  vector[T] v_err[J]; //Level innovation
  vector[N] neighbor_beta;
  real<lower=0> s_obs[J];
  real<lower=0> s_slope[J];
  real<lower=0> s_level[J];

  real mean_effect;
  real<lower=0> s_effect;
  vector[J] raw_effects;

//  real mean_colonizer;
  real<lower=0> s_colonizer;
  vector[C] raw_colonizer_effects;
}

transformed parameters {
  vector[T] u[J]; //Level
  vector[T] v[J]; //Slope
  vector[J] real_effects;
  vector[C] colonizer_effects;

  
  
  // state space
  for (j in 1:J) {
    u[j][1] = u_err[j][1];
    v[j][1] = v_err[j][1];
    for (t in 2:T) {
      u[j][t] = u[j][t-1] + v[j][t-1] + s_level[j] * u_err[j][t];
      v[j][t] = v[j][t-1] + s_slope[j] * v_err[j][t];
    }
  }

//  colonizer_effects = mean_colonizer + s_colonizer*raw_colonizer_effects;
//  real_effects = mean_effect + s_effect*raw_effects;
  colonizer_effects = mean_effect + s_effect*raw_colonizer_effects;
  for (j in 1:J) {
    real_effects[j] = colonizer_effects[colonizer[j]] + s_colonizer*raw_effects[j];
  }

}

model {
  mean_effect ~ normal(0, 0.1);
  s_effect ~ normal(0, 0.1);
  raw_effects ~ std_normal();

//  mean_colonizer ~ normal(0, 0.1);
  s_colonizer ~ normal(0, 0.1);
  raw_colonizer_effects ~ std_normal();
  

  neighbor_beta ~ std_normal();


  for (j in 1:J) {
    // state space priors
    u_err[j] ~ std_normal();
    v_err[j] ~ std_normal();
    s_obs[j] ~ std_normal();
    s_slope[j] ~ normal(0, 0.03);
    s_level[j] ~ normal(0, 0.03);


    // fitting state-space model to pre-intervention data
    y_norm[j][1:T0] ~ normal(u[j][1:T0] + X_norm[1:T0, indices[2*j-1]:indices[2*j]]*neighbor_beta[indices[2*j-1]:indices[2*j]], s_obs[j]);

    // fitting effect to make state-space prediction match observed data
    y_norm[j][T1:T] ~ normal(u[j][T1:T] + X_norm[T1:T, indices[2*j-1]:indices[2*j]]*neighbor_beta[indices[2*j-1]:indices[2*j]] + real_effects[j] / norms[j], s_obs[j]);
  }
}


