library(rstan)
library(bayesplot)
library(ggplot2)
library(CausalImpact)
library(gridExtra)
library(splines)
options(mc.cores = parallel::detectCores())



J <- 40
mean_effect = 0.04
effect_sizes <- rnorm(J, mean_effect, 0.02)
deviations <- runif(J, 0.01, 0.05)
T <- 20
T0 <- 10
FS = 3
N <- FS*J

X = matrix(0, nrow = T, ncol = N)
Y <- c()

indices_long <- c()

for (j in 1:J) {
  rn <- rnorm(T, 0, deviations[j])
  y <- cumsum(rn)
  x1 <- y + rnorm(T, 0, 0.5*deviations[j]) 
  x2 <- y + rnorm(T, 0, 0.3*deviations[j]) 
  x3 <- y + rnorm(T, 0, deviations[j]) 
  offset <- c(rep(0, T0), rep(effect_sizes[j], T - T0))
  x <- cbind(x1, x2, x3)
  y <- y + offset
  X[, (FS*(j-1)+1):(FS*j)] = x
  Y[[j]] <- y

  indices_long = c(indices_long, c(FS*(j-1)+1,FS*j))
}

X = apply(X,2 , norm <- function(x){return (x - mean(x))})
for (i in 1:J) {
  Y[[i]] = Y[[i]] - mean(Y[[i]][1:T0])
}


data <- list(
	T=T,
	T0=T0,
	N=N,
  J=J,
	X=X,
	y=Y,
  indices=indices_long
)

data2 <- list(
  T=T,
  T0=T0,
  N=N*J,
  M=N,
  J=J,
  X=X,
  y=Y,
  indices=indices_long,
  indices_short=indices_long
)


sm <- stan_model('stan_models/multiple_lines.stan')
fit <- sampling(sm, data=data, iter=2000, control=list(adapt_delta=0.9))

mean_effect_draws = apply(as.array(fit, pars = c('real_effects')), 3, mean)
effect_df = data.frame(cbind(mean_effect_draws, effect_sizes))
ggplot(effect_df, aes(x=effect_sizes, y=mean_effect_draws)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  xlab('Effect') + 
  ylab('Modeled Effect') + 
  ggtitle('Validating Model Effect Inference')

individual_effect_draws = as.matrix(fit, pars = c('real_effects'))
apply(sweep(individual_effect_draws, 2, effect_sizes), 2, mean) - (apply(individual_effect_draws, 2, mean) - effect_sizes)
mcmc_intervals(sweep(individual_effect_draws, 2, effect_sizes))


traceplot(fit, pars=c('mean_effect', 'real_effects[1]', 's_obs[1]', 's_effect', 's_level[1]'))

pairs(fit, pars=c('mean_effect', 'effects[1]', 's_obs[1]', 's_effect', 's_level[1]', 'neighbor_beta[1]'))
pairs(fit, pars=c('mean_effect', 's_effect', 'effects[1]', 'effects[2]', 'effects[10]', 'effects[20]'))
pairs(fit, pars = c('neighbor_beta[1]', 'neighbor_beta[2]', 'neighbor_beta[4]', 'neighbor_beta[3]', 'neighbor_beta[10]', 'neighbor_beta[30]'))
pairs(fit, pars = c('s_obs[1]', 's_obs[2]', 's_obs[4]', 's_obs[3]', 's_obs[10]', 's_obs[20]'))
pairs(fit, pars=c('s_slope[1]', 's_level[1]', 's_slope[2]', 's_level[2]', 's_slope[20]', 's_level[20]'))
pairs(fit, pars=c('u_err[1,1]', 'v_err[1,1]', 'u_err[1,2]', 'v_err[1,2]', 'u_err[1,20]', 'v_err[1,20]'))
pairs(fit, pars=c('u_err[2,1]', 'v_err[2,1]', 'u_err[2,2]', 'v_err[2,2]', 'u_err[2,3]', 'v_err[2,3]', 'u_err[2,4]', 'v_err[2,4]', 'u_err[2,20]', 'v_err[2,20]'))
pairs(fit, pars=c('u_err[10,1]', 'v_err[10,1]', 'u_err[10,2]', 'v_err[10,2]', 'u_err[10,3]', 'v_err[10,3]', 'u_err[10,4]', 'v_err[10,4]', 'u_err[10,20]', 'v_err[10,20]'))


draws <- as.array(fit, pars = c('mean_effect', 's_level[1]', 's_level[2]', 's_effect'))

draws <- as.array(fit, pars = c('mean_effect', 's_effect', 'real_effects[2]', 'real_effects[10]', 'real_effects[1]'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)


mean_effect_draws = as.array(fit, pars = c('mean_effect'))
mcmc_dens(mean_effect_draws, pars=c('mean_effect')) + 
  geom_vline(xintercept=mean_effect, color='black') + 
  xlim(mean_effect - 0.1, mean_effect + 0.1) +
  annotate("text", x=mean_effect-0.005, y=30, label='True Value', colour='black', angle=90, size=5) +
  xlab('Effect Size') + 
  ylab('Density') + 
  ggtitle('Validating Mean Effect Inference')



dev.new()

sm2 <- stan_model('stan_models/multiple_lines_qr.stan')
fit2 <- sampling(sm2, data=data2, iter=2000, control=list(adapt_delta=0.9))
mean_effect_draws = as.array(fit2, pars = c('m_effect'))
mcmc_dens(mean_effect_draws, pars=c('m_effect'))


individual_effect_draws = as.matrix(fit, pars = c('effect', 'm_effect'))
mcmc_intervals(individual_effect_draws)

dev.new()

individual_effect_draws = as.matrix(fit2, pars = c('effect', 'm_effect'))
mcmc_intervals(individual_effect_draws)


mcmc_scatter(
  as.matrix(fit),
  pars = c('neighbor_beta[1]', 'neighbor_beta[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

dev.new()

mcmc_scatter(
  as.matrix(fit2),
  pars = c('neighbor_beta[1]', 'neighbor_beta[2]'),
  np = nuts_params(fit2), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

dev.new()

mcmc_scatter(
  as.matrix(fit2),
  pars = c('beta_tilde[1]', 'beta_tilde[2]'),
  np = nuts_params(fit2), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)




mcmc_scatter(
  as.matrix(fit),
  pars = c('neighbor_beta[3]', 'neighbor_beta[4]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

dev.new()

mcmc_scatter(
  as.matrix(fit2),
  pars = c('neighbor_beta[3]', 'neighbor_beta[4]'),
  np = nuts_params(fit2), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

dev.new()

mcmc_scatter(
  as.matrix(fit2),
  pars = c('beta_tilde[3]', 'beta_tilde[4]'),
  np = nuts_params(fit2), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)



plot(fit, pars=c('beta', 's_obs', 's_slope', 's_level', 'effect'))
plot(fit, pars=c('effect', 'm_effect'))
traceplot(fit, pars=c('beta', 's_obs', 's_slope', 's_level', 'effect'))
pairs(fit, pars=c('beta', 's_obs', 's_slope', 's_level', 'effect'))
stan_diag(fit, info = 'sample')


y_rep <- extract(fit)$y_rep
y_rep_m <- apply(y_rep, 2, mean)
y_quantiles <- t(apply( y_rep, 2 , quantile, probs = c(0.05, 0.95)))
effect <- y - y_rep_m
results <- as.data.frame(cbind(y, effect, y_rep_m, y_quantiles))
results$idu <- as.numeric(row.names(results))
colnames(results) <- c('y', 'effect', 'y_rep', 'y_min', 'y_max', 'idu')
results$effect_min <- results$y - results$y_min
results$effect_max <- results$y - results$y_max
results$cum_effect <- cumsum(results$effect)
results$cum_effect_min <- c(rep(0, T0), cumsum(results[T0+1:T, 'effect_min'])[1:T0])
results$cum_effect_max <- c(rep(0, T0), cumsum(results[T0+1:T, 'effect_max'])[1:T0])


p_model <- ggplot(results, aes(idu, y_rep)) +
    geom_line(linetype='dashed', color='blue') +
    geom_line(aes(idu, y), color='black') +
    geom_ribbon(aes(x=idu, ymin=y_min,ymax=y_max),alpha=0.3, color='blue')

p_effect <- ggplot(results, aes(idu, effect)) +
    geom_line(linetype='dashed', color='blue') +
    geom_ribbon(aes(x=idu, ymin=effect_min,ymax=effect_max),alpha=0.3, color='blue')

p_cum_effect <- ggplot(results, aes(idu, cum_effect)) +
    geom_line(linetype='dashed', color='blue') +
    geom_ribbon(aes(x=idu, ymin=cum_effect_min,ymax=cum_effect_max),alpha=0.3, color='blue')

grid.arrange(p_model, p_effect, p_cum_effect, nrow=3)


color_scheme_set('darkgray')
f <- extract(fit)
ppc_dens_overlay(y = y, yrep = f$effect[1:T0, ])

draws <- as.array(fit, pars = c('beta', 's_obs', 's_slope', 's_level', 'effect'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('u'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('effect'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('spline_hat'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

mcmc_scatter(
  as.matrix(fit),
  pars = c('s_obs', 's_slope'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)




