library(rstan)
library(bayesplot)
library(ggplot2)
library(CausalImpact)
library(gridExtra)
options(mc.cores = parallel::detectCores())


effect_size <- 2

T <- 100
T0 <- 50
N <- 2

rn <- rnorm(100, 0, 1)
y <- cumsum(rn)
x1 <- y + rnorm(100, 0, 0.1) 
x2 <- y + rnorm(100, 0, 0.1) 
offset <- c(rep(0, 50), rep(effect_size, 50))
x <- cbind(x1, x2)
y <- y + offset


data <- list(
	T=T,
	T0=T0,
	N=N,
	X=x,
	y=y
)

sm <- stan_model('stan_models/single_line.stan')
fit <- sampling(sm, data=data, iter=2000, control=list(adapt_delta=0.9))
mean_effect_draws = as.array(fit, pars = c('effect'))
mcmc_dens(mean_effect_draws, pars=c('effect'))

dev.new()

sm2 <- stan_model('stan_models/single_line_qr.stan')
fit2 <- sampling(sm2, data=data, iter=2000, control=list(adapt_delta=0.9))
mean_effect_draws = as.array(fit2, pars = c('effect'))
mcmc_dens(mean_effect_draws, pars=c('effect'))



plot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
dev.new()
plot(fit2, pars=c('beta', 's_obs', 's_slope', 's_level'))

traceplot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
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
ppc_dens_overlay(y = y, yrep = f$y_rep[1:T0, ])

draws <- as.array(fit, pars = c('beta', 's_obs', 's_slope', 's_level'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

mcmc_scatter(
  as.matrix(fit),
  pars = c('beta[1]', 'beta[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)
dev.new()
mcmc_scatter(
  as.matrix(fit2),
  pars = c('beta[1]', 'beta[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)
dev.new()
mcmc_scatter(
  as.matrix(fit2),
  pars = c('beta_tilde[1]', 'beta_tilde[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)











set.seed(1)
x1 <- 100 + arima.sim(model = list(ar = 0.999), n = 100)
y <- 1.2 * x1 + rnorm(100)
y[71:100] <- y[71:100] + 10
data_impact <- cbind(y, x1)
pre.period <- c(1, 70)
post.period <- c(71, 100)
impact <- CausalImpact(data_impact, pre.period, post.period, model.args=list(prior.level.sd=0.1))
plot(impact)

dev.new()


T <- 100
T0 <- 70
data <- list(
	T=100,
	T0=70,
	N=1,
	X=matrix(x1),
	y=y
)

sm <- stan_model('single_series.stan')

fit <- sampling(sm, data=data, iter=2000)
plot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
traceplot(fit, pars=c('beta', 's_obs', 's_slope', 's_level'))
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
results$cum_effect_min <- c(rep(0, T0), cumsum(results[T0+1:T, 'effect_min'])[1:(T-T0)])
results$cum_effect_max <- c(rep(0, T0), cumsum(results[T0+1:T, 'effect_max'])[1:(T-T0)])


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