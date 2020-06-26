library(rstan)
library(bayesplot)
library(ggplot2)
library(CausalImpact)
library(gridExtra)
library(splines)
options(mc.cores = parallel::detectCores())


effect_size <- 10

T <- 100
T0 <- 70
N <- 1

rn <- rnorm(T, 0, 1)
y <- cumsum(rn)
x1 <- y + rnorm(T, 0, 0.1) 
# x2 <- y + rnorm(T, 0, 0.6) 
offset <- c(rep(0, T0), rep(effect_size, T - T0))
# x <- cbind(x1, x2)
x <- matrix(x1)
y <- y + offset

num_knots <- 10
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
spline_points <- seq(from=1, to=T, by=1)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))


data <- list(
	T=T,
	T0=T0,
	N=N,
	X=x,
	y=y,
  num_knots=num_knots,
  knots=knots,
  spline_degree=spline_degree,
  spline_points=spline_points
)


sm <- stan_model('stan_models/single_spline.stan')

fit <- sampling(sm, data=data, iter=2000)
plot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
traceplot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
pairs(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
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

draws <- as.array(fit, pars = c('beta', 's_obs', 's_slope', 's_level'))
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




