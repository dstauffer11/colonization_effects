library(CausalImpact)
library(ggplot2)
library(rjson)
library(rstan)
library(bayesplot)
library(ggplot2)
library(CausalImpact)
library(gridExtra)
library(splines)
options(mc.cores = parallel::detectCores())



country.year.neighbors.list <- fromJSON(file='data/country_year_neighbors_map.json')
gdp.data <- read.csv('data/annual_growth.csv')

T <- 20
T0 <- 10
J <- length(country.year.neighbors.list)

Y <- c()
indices_long <- c()
indices_short = c()
codes <- c()
M = 19;


X = matrix(0, nrow = 20, ncol = M*J)

index = 0
tl = 0
for (country_code in names(country.year.neighbors.list)) {

	index = index + 1

	yin <- country.year.neighbors.list[[country_code]]
	codes <- c(codes, country_code)
	

	first.year <- strtoi(yin[[1]])
	independence.year <- strtoi(yin[[2]])
	last.year <- strtoi(yin[[3]])

	if ((first.year > independence.year - 9) | (last.year < independence.year + 10)) next

	first.year <- independence.year - 9
	last.year <- independence.year + 10

	print(country_code)

	countries <- c(country_code, yin[[4]])
	local.gdp.data <- gdp.data[gdp.data$country_code %in% countries, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year >= first.year, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year <= last.year, ]
	
	y <- local.gdp.data[local.gdp.data$country_code == country_code, 'annual_growth']
	Y <- append(Y, list(y))

	x.long <- local.gdp.data[local.gdp.data$country_code %in% yin[[4]], ]
	x.long <- subset(x.long, select=c('year', 'country_code', 'annual_growth'))
	x.wide <- reshape(x.long, direction = "wide", idvar = 'year', timevar = 'country_code')

	x <- data.matrix(x.wide[ , !(names(x.wide) %in% c('year'))])

	n = min(length(yin[[4]]), M)
	si = M*(index - 1) + 1
	ei = M*(index - 1) + n
	indices_long = c(indices_long, c(si, ei))
	indices_short = c(indices_short, c(tl+1, tl+n))
	tl = tl+n
	X[1:T, si:ei] <- x[, 1:n]

}

X = apply(X,2 , norm <- function(x){return ((x - mean(x)) / sd(x))})
X[is.na(X)] = 0
for (i in 1:J) {
	Y[[i]] = (Y[[i]] - mean(Y[[i]][1:T0])) / sd(Y[[i]][1:T0])
}


for (i in 1:(M*J-1)) {
	print(i)
	if (sum(X[, i]) > 0 & sum(X[, i+1]) > 0) {
		print(cor(X[, i], X[, i+1]))
	}
}



data <- list(
	T=T,
	T0=T0,
	J=J,
	X=X,
	y=Y,
	N=tl,
	M=M,
	indices=indices_long,
	indices_short=indices_short
)

print(codes)


sm <- stan_model('stan_models/multiple_lines.stan')

fit <- sampling(sm, data=data, iter=4000, control=list(adapt_delta=0.9))

# check tree depths
check_treedepth(fit)
# check energy
check_energy(fit)


# Histogram of tree depths
breaks = 0:10
sampler_params = get_sampler_params(fit, inc_warmup=FALSE)
treedepths = do.call(rbind, sampler_params)[, 'treedepth__']
treedepths_hist = hist(treedepths, breaks=breaks, plot=FALSE)
par(mar=c(4, 4, 0.5, 0.5))
plot(treedepths_hist, main='', xlab='theta.1', yaxt='n', ann=FALSE)

# Estimate distribution of mean effect from posterior samples
mean_effect_draws = as.array(fit, pars = c('m_effect'))
mcmc_dens(mean_effect_draws, pars=c('m_effect'))

# Estimate distribution of each indiviadual country's effects
individual_effect_draws = as.matrix(fit, pars = c('effect', 'm_effect'))
mcmc_intervals(individual_effect_draws) + scale_y_discrete(labels=c(codes, 'Overall'))

# check step sizes of sample
sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
stepsizes <- sapply(sampler_params, function(x) x[1,'stepsize__'])
names(stepsizes) <- list("Chain 1", "Chain 2", "Chain 3" ,"Chain 4")
stepsizes
mean(stepsizes)
# 2000 iterations, centered: 0.009013763
# 2000 iteratoins, centered, qr: 0.01270642

# check gradient evaluations
sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients
sum(n_gradients)
# 2000 ierations, centered: 2318880
# 2000 iteratons, centered, qr: 2244971

traceplot(fit, pars=c('beta_tilde[4]', 's_obs[30]', 's_slope[2]', 's_level[3]', 'effect[20]'))

pairs(fit, pars=c('s_obs[1]', 's_obs[2]', 's_s_obs', 's_n_beta', 's_u_err', 's_v_err', 's_effect', 's_s_slope', 's_s_level'))


pairs(fit, pars=c('beta_tilde[4]', 's_s_obs', 's_s_slope', 's_s_level', 'effect[35]'))

stan_diag(fit, info = 'sample')


color_scheme_set('darkgray')
f <- extract(fit)
ppc_dens_overlay(y = y, yrep = f$effect[1:T0, ])

draws <- as.array(fit, pars = c('s_obs'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('s_slope'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('s_level'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('u_err'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('v_err'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)


draws <- as.array(fit, pars = c('u'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('effect'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('beta_tilde'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('neighbor_beta'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)


mcmc_scatter(
  as.matrix(fit),
  pars = c('s_obs', 's_slope'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)


mcmc_scatter(
  as.matrix(fit),
  pars = c('neighbor_beta[1]', 'neighbor_beta[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

dev.new()

mcmc_scatter(
  as.matrix(fit),
  pars = c('beta_tilde[1]', 'beta_tilde[2]'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)


mean_effect_draws = as.array(fit, pars = c('s_effect'))
mcmc_dens(mean_effect_draws, pars=c('s_effect'))



