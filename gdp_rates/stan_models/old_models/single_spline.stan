functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}



data {
  int<lower=2> T;
  int<lower=1> T0;
  int<lower=1> N;
  matrix[T, N] X;
  vector[T] y;

  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  real spline_points[T];
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, T] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(spline_points, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, T] = 1;
}


parameters {
  vector[T] u_err; //Slope innovation
  vector[T] v_err; //Level innovation
  vector[N] beta;
  real<lower=0> s_obs;
  real<lower=0> s_slope;
  real<lower=0> s_level;

  row_vector[num_basis] a_raw;
  real a0;  // intercept
  real<lower=0> tau;
  real<lower=0> s_spline;
}

transformed parameters {
  vector[T] u; //Level
  vector[T] v; //Slope

  row_vector[num_basis] a; // spline coefficients
  vector[T] spline_hat;

  
  // state space
  u[1] = u_err[1];
  v[1] = v_err[1];
  for (t in 2:T) {
    u[t] = u[t-1] + v[t-1] + s_level * u_err[t];
    v[t] = v[t-1] + s_slope * v_err[t];
  }

  // spline values
  a[1] = a_raw[1];
  for (i in 2:num_basis)
    a[i] = a[i-1] + a_raw[i]*tau;
  spline_hat = a0*to_vector(spline_points) + to_vector(a*B);
}

model {
  // state space priors
  u_err ~ normal(0, 1);
  v_err ~ normal(0, 1);
  beta ~ normal(0, 1);
  s_obs ~ normal(0, 5);
  s_slope ~ normal(0, 0.1);
  s_level ~ normal(0, 0.1);

  // spline priors
  a_raw ~ normal(0, 1);
  a0 ~ normal(0, 1);
  tau ~ normal(0, 1);
  s_spline ~ normal(0, 5);

  // fitting state-space model to pre-intervention data
  y[1:T0] ~ normal(u[1:T0] + X[1:T0, :]*beta, s_obs);

  // fitting spline as effect to make state-space prediction match observed data
  y ~ normal(u + X*beta - spline_hat, s_spline);


}


