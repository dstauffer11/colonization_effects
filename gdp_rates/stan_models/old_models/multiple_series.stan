data {
  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> N;
  matrix[T, N] X;
  vector[T] y;
}

parameters {
  vector[T] u_err; //Slope innovation
  vector[T] v_err; //Level innovation
  vector[N] beta;
  real<lower=0> s_obs;
  real<lower=0> s_slope;
  real<lower=0> s_level;
}

transformed parameters {
  vector[T0] u; //Level
  vector[T0] v; //Slope
  real[T-T0] benefit;
  u[1] = u_err[1];
  v[1] = v_err[1];
  for (t in 2:T0) {
    u[t] = u[t-1] + v[t-1] + s_level * u_err[t];
    v[t] = v[t-1] + s_slope * v_err[t];
  }


}

model {
  u_err ~ normal(0,1);
  v_err ~ normal(0,1);
  beta ~ normal(0, 1);
  s_obs ~ normal(0, 5);
  s_slope ~ normal(0, 5);
  s_level ~ normal(0, 5);

  y[1:T0] ~ normal(u + X[1:T0, :]*beta, s_obs);
}

generated quantities {
  vector[T] y_rep;
  vector[T-T0] u_rep;
  vector[T-T0] v_rep;

  for (t in 1:T0) {
    y_rep[t] = normal_rng(u[t] + X[t, :]*beta, s_obs);
  }
  

  
  u_rep[1] = u[T0] + v[T0] + s_level*u_err[T0];
  v_rep[1] = v[T0] + s_slope * v_err[T0];
  y_rep[T0+1] = normal_rng(u_rep[1] + X[T0+1, :]*beta, s_obs);
  for (t in 2:T-T0) {
    u_rep[t] = u_rep[t-1] + v_rep[t-1] + s_level * u_err[t+T0];
    v_rep[t] = v_rep[t-1] + s_slope * v_err[t+T0];
    y_rep[t+T0] = normal_rng(u_rep[t] + X[t+T0, :]*beta, s_obs);
  }
}
