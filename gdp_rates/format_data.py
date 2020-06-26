# map each country with neighbors - done
# plot model coefficients - done
# improve model if necessary:
#   include gap over time as benefit? - no
#   model correlation of countries' benefits by comparing independence dates, overlap of shared countries - maybe... not
#   normalize/decorrelate predictors - done
#   QR decomposition of predictors - done
#   aggregate effects across all countries and by former colonizer (heirarchical) - done
#   fix beta vector to not have extra zeros - done
#   sample standard deviation for effects with transformed cauchy

# plot gdp growth time series - done



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import seaborn as sns
sns.color_palette("muted")
from neighbors import Neighbors, CountryCodes
world = Neighbors()
cc = CountryCodes()




independence_df = pd.read_csv('data/independence_data.csv', header=0, delimiter=', ')
independence_df['country_code'] = independence_df.country.map(cc.name_to_code)
independence_df['colonizer_code'] = independence_df.colonizer.map(cc.name_to_code)



gdp_wide = pd.read_excel('data/mpd2018.xlsx', sheet_name='rgdpnapc', header=1)
gdp = gdp_wide.set_index('year').stack().reset_index(name='gdp_per_capita').rename(columns={'level_1':'country_code'})
gdp = gdp.rename(columns={'rgdpnapc': 'year'})


gdp = gdp.sort_values(['country_code', 'year'])
gdp_rate = gdp.groupby(['country_code']) \
	.apply(lambda g: pd.DataFrame(data={
		'year': g.year[1:], 
		'gpc_growth': (np.log(g.gdp_per_capita.values[1:]) - np.log(g.gdp_per_capita.values[:-1])),
		'year_gap': g.year.values[1:] - g.year.values[:-1],
		})) \
	.reset_index()
gdp_rate['annual_growth'] = gdp_rate.gpc_growth / gdp_rate.year_gap
gdp_rate_wide = pd.pivot_table(gdp_rate, index='year', columns='country_code', values='annual_growth').reset_index()


def last_year(df, country_code):
	return int(df[(df.year_gap == 1) & (df.country_code == country_code)].year.max())

def first_year(df, country_code):
	ly = last_year(df, country_code)
	cdf = df[df.country_code == country_code]
	for year in cdf[(cdf.year_gap == 1)].sort_values('year').year:
		if (ly - year) - len(cdf[(cdf.year >= year) & (cdf.year <= ly)]) < 5:
			return year

def interpolate_gdp_rates(df, country_code, fy, ly):
	cdf = df[df.country_code == country_code]
	df = df[df.country_code != country_code]
	cdf_yr = cdf.set_index('year').reindex(list(range(fy, ly+1))).reset_index().sort_values('year')
	cdf_yr['annual_growth'] = cdf_yr.annual_growth.interpolate(method='linear')
	# cdf = pd.concat([cdf[cdf.year < fy], cdf_yr])
	df = pd.concat([df, cdf])
	return df

for country_code in set(gdp_rate.country_code):
	fy = first_year(gdp_rate, country_code)
	ly = last_year(gdp_rate, country_code)
	gdp_rate = interpolate_gdp_rates(gdp_rate, country_code, fy, ly)


country_year_neighbors_map = {}
for _, row in independence_df.iterrows():
	country_code = row.country_code
	year = row.independence_year
	if country_code not in set(gdp_rate.country_code):
		# print('{} not found in GDP data.'.format(row.country))
		continue
	elif first_year(gdp_rate, country_code) > year - 9:
		# print('Not enough GDP data for {}: GDP starts in {}, independence in {}'.format(row.country, first_year(gdp_rate, country_code), year))
		continue
	else:
		neighbors = world.gdp_neighbors(country_code)
		gdp_neighbors = []
		fy = first_year(gdp_rate, country_code)
		ly = last_year(gdp_rate, country_code)

		for neighbor in neighbors:

			if neighbor in independence_df.country_code and abs(independence_df.set_index('country_code').loc[neighbor].independence_year - year) < 5:
				continue
			if neighbor in set(gdp_rate.country_code) and first_year(gdp_rate, neighbor) <= year - 9:
				fy = max(fy, first_year(gdp_rate, neighbor))
				gdp_neighbors.append(neighbor)


		if len(gdp_neighbors) == 0:
			# print('No good neighbors found for {}'.format(row.country))
			continue
		else:
			# print(country_code, fy, year, ly)
			country_year_neighbors_map[country_code] = (fy, year, ly, gdp_neighbors, row.colonizer_code)

assert(gdp_rate[pd.isnull(gdp_rate.annual_growth)].shape[0] == 0)

def plot_country_growth(country_code, neighbors, fy, year, ly, gdp_df, gdp_col, show=False, log=True):
	local_gdp = gdp_df[gdp_df.country_code == country_code]
	scaling = np.log if log else lambda x: x

	fig, ax = plt.subplots(figsize=(9, 5))
	ax.plot(local_gdp.year, local_gdp[gdp_col].map(scaling), 
		c='r', linestyle='--', marker='.', label=cc.code_to_name(country_code), zorder=5)

	for neighbor in neighbors:
		neighbor_gdp = gdp_df[gdp_df.country_code == neighbor]
		ax.plot(neighbor_gdp.year, neighbor_gdp[gdp_col].map(scaling), 
			linestyle='--', marker='.', label=cc.code_to_name(neighbor), alpha=0.4, zorder=0)

	ax.axvline(year, color='k')
	ax.legend()
	ax.grid(zorder=10, alpha=0.5)
	ax.set_xlim([fy, ly])
	ax.set_title('{} with Neighbors {}'.format(country_code, ', '.join(neighbors)))
	ax.set_xlabel('Year')
	if show:
		plt.show()
	return ax


# def pc(country, gdp_rate):
# 	local_gdp_rate = gdp_rate[gdp_rate.country == country]
# 	fig, ax = plt.subplots(figsize=(9, 5))
# 	ax.plot(local_gdp_rate.year, local_gdp_rate.annual_growth, c='r', linestyle='-', label=country)
# 	ax.set_xlim(1800)
# 	plt.show()

# pc('United States', gdp_rate)

# fig, ax = plt.subplots(figsize=(9, 5))
# plt.hist([ year - fy for (fy, year, ly, _, _) in country_year_neighbors_map.values()], bins=np.linspace(0, 50, 51))
# plt.show()
# fig, ax = plt.subplots(figsize=(9, 5))
# plt.hist([ ly - year for (fy, year, ly, _, _) in country_year_neighbors_map.values()], bins=np.linspace(0, 50, 51))
# plt.show()

for country_code, (fy, year, ly, gdp_neighbors, _) in country_year_neighbors_map.items():
	# http://www.futurile.net/2016/02/27/matplotlib-beautiful-plots-with-style/
	ax = plot_country_growth(country_code, gdp_neighbors, fy, year, ly, gdp, 'gdp_per_capita')
	ax.set_ylabel('sLog(GDP per Capita)')
	plt.savefig('results/time_series/gdp/{}.png'.format(country_code))
	plt.close(ax.figure)

	ax = plot_country_growth(country_code, gdp_neighbors, fy, year, ly, gdp_rate, 'annual_growth', log=False)
	ax.set_ylabel('First Difference in Log(GDP per Capita)')
	plt.savefig('results/time_series/gdp_growth/{}.png'.format(country_code))
	plt.close(ax.figure)

	ax = world.map_countries(country_code, gdp_neighbors, gdp_rate_wide, fy, year)
	plt.savefig('results/maps/{}.png'.format(country_code))
	plt.close(ax.figure)

print('Countries to Analyze: {}'.format(', '.join(country_year_neighbors_map.keys())))

with open('data/country_year_neighbors_map.json', 'w') as f:
	json.dump(country_year_neighbors_map, f)

gdp_rate.to_csv('data/annual_growth.csv')
