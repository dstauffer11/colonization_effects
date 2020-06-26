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

N <- 0
X <- matrix(, nrow = 20, ncol = 0)
Y <- c()
indices <- c()
codes <- c()

neighbor_map = list()
colonizers = c()

for (country_code in names(country.year.neighbors.list)) {

	yin = country.year.neighbors.list[[country_code]]

	codes = c(codes, country_code)
	neighbor_map[[country_code]] = yin[[4]]
	colonizers = c(colonizers, yin[[5]])


	first.year <- strtoi(yin[[1]])
	independence.year <- strtoi(yin[[2]])
	last.year <- strtoi(yin[[3]])

	if ((first.year > independence.year - 9) | (last.year < independence.year + 10)) next

	first.year <- independence.year - 9
	last.year <- independence.year + 10

	countries <- c(country_code, yin[[4]])
	local.gdp.data <- gdp.data[gdp.data$country_code %in% countries, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year >= first.year, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year <= last.year, ]
	
	y <- local.gdp.data[local.gdp.data$country_code == country_code, 'annual_growth']
	Y <- append(Y, list(y))

	indices <- c(indices, c(N+1, N+length(yin[[4]])))
	N <- N + length(yin[[4]])

	x.long <- local.gdp.data[local.gdp.data$country_code %in% yin[[4]], ]
	x.long <- subset(x.long, select=c('year', 'country_code', 'annual_growth'))
	x.wide <- reshape(x.long, direction = "wide", idvar = 'year', timevar = 'country_code')

	x <- x.wide[ , !(names(x.wide) %in% c('year'))]

	X <- cbind(X, data.matrix(x))

}

X = apply(X, 2 , norm <- function(x){return (x - mean(x))})
for (i in 1:J) {
	Y[[i]] = Y[[i]] - mean(Y[[i]][1:T0])
}

colonizers_numeric = as.numeric(factor(colonizers))
C = max(colonizers_numeric)

data = list(
	T=T,
	T0=T0,
	N=N,
	J=J,
	X=X,
	y=Y,
	indices=indices
)

sm <- stan_model('stan_models/multiple_lines.stan')
fit <- sampling(sm, data=data, iter=20000, control=list(adapt_delta=0.995, max_treedepth=13), seed=1, chains=1, warmup=5000)

sm <- stan_model('stan_models/multiple_lines_simp.stan')
fit <- sampling(sm, data=data, iter=2000, control=list(adapt_delta=0.9, max_treedepth=10), seed=2)

print(fit, pars=c('mean_effect', 's_effect', 'raw_effects', 'real_effects', 's_s_obs', 's_obs'))

data_colonizers = list(
	T=T,
	T0=T0,
	N=N,
	J=J,
	X=X,
	y=Y,
	indices=indices,
	C=C,
	colonizer=colonizers_numeric
)



smc <- stan_model('stan_models/multiple_lines_colonizers2.stan')
fitc <- sampling(smc, data=data_colonizers, iter=2000, control=list(adapt_delta=0.9, max_treedepth=10), seed=5)

# check tree depths
check_treedepth(fit)
# check energy
check_energy(fit)

# check chain mixing
traceplot(fit, pars=c('mean_effect', 's_effect', 'real_effects[1]', 's_obs[1]', 's_effect', 's_level[1]', 'neighbor_beta[10]'))
traceplot(fit, pars=c('neighbor_beta[1]', 'neighbor_beta[10]', 'neighbor_beta[50]', 'neighbor_beta[100]', 'neighbor_beta[200]', 'neighbor_beta[300]'))
traceplot(fit, pars=c('u[1,1]', 'u[1,15]', 'u[10,3]', 'u[10,5]', 'u[30,5]', 'u[40,20]'))
traceplot(fit, pars=c('mean_effect', 's_colonizer_effect_unif', 's_colonizer_effect', 'raw_colonizer_effects[1]', 
	'raw_colonizer_effects[2]', 'raw_colonizer_effects[3]', 'raw_colonizer_effects[4]', 'raw_colonizer_effects[5]'))

pairs(fit, pars=c('mean_effect', 's_effect', 'real_effects[1]', 'real_effects[10]', 'real_effects[30]'))
pairs(fit, pars=c('mean_effect', 's_effect', 's_obs[1]', 's_obs[10]', 's_obs[30]', 's_s_obs'))

pairs(fitc, pars=c('mean_effect', 's_effect', 's_s_obs', 's_colonizer_effect', 'raw_colonizer_effects'))


# Histogram of tree depths
breaks = 0:13
sampler_params = get_sampler_params(fit, inc_warmup=FALSE)
treedepths = do.call(rbind, sampler_params)[, 'treedepth__']
treedepths_hist = hist(treedepths, breaks=breaks, plot=FALSE)
par(mar=c(4, 4, 0.5, 0.5))
plot(treedepths_hist, main='', xlab='theta.1', yaxt='n', ann=FALSE)

# Estimate distribution of mean effect from posterior samples
mean_effect_draws = as.array(fit, pars = c('mean_effect'))
mcmc_dens(mean_effect_draws, pars=c('mean_effect')) +
	ggtitle('Estimated Effect of Indepence on Economic Growth') +
	xlab('Change in Growth Rate of GDP per Capita (2011 US$)') + 
	ylab('Density')

# Estimate effect of each colonizer
colonizers_reducted = c('BEL', 'DEU', 'FRA', 'GBR', 'PRT', 'RUS')
mean_effect_draws = as.array(fitc, pars = c('real_colonizer_effects'))
mcmc_areas(-mean_effect_draws, regex_pars = "real_colonizer_effects\\[[1-6]\\]", prob=0.8) +
	ggtitle('Estimated Effect of Colonizer its Colonies') +
	xlab('Depression of Growth Rate') + 
	ylab('Density') + 
	xlim(-0.05, 0.05) +
	scale_y_discrete(labels=colonizers_reducted)

# Estimate distribution of each indiviadual country's effects
individual_effect_draws = as.matrix(fit, pars = c('real_effects', 'mean_effect'))
mcmc_intervals(individual_effect_draws) + 
	scale_y_discrete(labels=c(codes, 'Overall')) +
	ggtitle('Estimated Independence Effect by Country') +
	xlab('Growth Rate Change') + 
	theme(
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                colour = "gray"),
  panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                colour = "gray")
  )


# check step sizes of sample
sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
stepsizes <- sapply(sampler_params, function(x) x[1,'stepsize__'])
names(stepsizes) <- list("Chain 1", "Chain 2", "Chain 3" ,"Chain 4")
stepsizes
mean(stepsizes)
# 2000 iterations, centered: 0.009013763

# check gradient evaluations
n_gradients <- sapply(sampler_params, function(x) sum(x[,'n_leapfrog__']))
n_gradients
sum(n_gradients)
# 2000 ierations, centered: 2318880



stan_diag(fit, info = 'sample')


color_scheme_set('darkgray')
f <- extract(fit)
ppc_dens_overlay(y = y, yrep = f$effect[1:T0, ])

draws <- as.array(fit, pars = c('neighbor_beta', 's_obs', 's_slope', 's_level', 'mean_effect'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('u'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('mean_effect', 's_effect', 'raw_effects'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)

draws <- as.array(fit, pars = c('mean_effect', 's_effect', 's_obs'))
np <- nuts_params(fit)
mcmc_parcoord(draws, np=np)



mcmc_scatter(
  as.matrix(fit),
  pars = c('s_obs', 's_slope'),
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)











measure.impact <- function(country_code, yin) {
	first.year <- strtoi(yin[[1]])
	independence.year <- strtoi(yin[[2]])
	last.year <- strtoi(yin[[3]]) - 1
	countries <- c(country_code, yin[[4]])
	local.gdp.data <- gdp.data[gdp.data$country_code %in% countries, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year >= first.year, ]
	local.gdp.data <- local.gdp.data[local.gdp.data$year <= last.year, ]
	t <- as.Date(ISOdate(local.gdp.data[local.gdp.data$country_code == country_code, 'year'], 1, 1))
	
	y <- local.gdp.data[local.gdp.data$country_code == country_code, 'annual_growth']

	x.long <- local.gdp.data[local.gdp.data$country_code %in% yin[[4]], ]
	x.long <- subset(x.long, select=c('year', 'country_code', 'annual_growth'))
	x.wide <- reshape(x.long, direction = "wide", idvar = 'year', timevar = 'country_code')

	x <- x.wide[ , !(names(x.wide) %in% c('year'))]

	# print(x)
	
	simple.data <- zoo(cbind(y, x), t)


	pre.period <- as.Date(ISOdate(c(first.year, independence.year), 1, 1))
	post.period <- as.Date(ISOdate(c(independence.year+1, last.year), 1, 1))
	impact <- CausalImpact(simple.data, pre.period, post.period)


	plot(impact)
	ggsave(sprintf('results/causal_impact/%s.png', country_code))

	return(impact)
}

# measure.impact('Canada', c(1901, 1931, c('United States', 'Mexico'))) 
for (name in names(country.year.neighbors.list)) {
	print(name)
	measure.impact(name, country.year.neighbors.list[[name]])
}
