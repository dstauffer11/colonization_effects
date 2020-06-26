import osgeo.ogr
import shapely.wkt
import pyproj
import country_converter as coco
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
plt.style.use('ggplot')

class CountryCodes:

	def name_to_code(self, name):
		return coco.convert(names=[name], to='ISO3')

	def code_to_name(self, code):
		return coco.convert(names=[code], to='name_short')


class Neighbors:

	countries = {} # Maps country name to Shapely geometry.
	geod = pyproj.Geod(ellps='WGS84')
	cc = CountryCodes()

	def __init__(self):
		shapefile = osgeo.ogr.Open("world_map_data/TM_WORLD_BORDERS-0.3.shp")
		layer = shapefile.GetLayer(0)

		for i in range(layer.GetFeatureCount()):
			feature = layer.GetFeature(i)
			country = self.cc.name_to_code(feature.GetField("NAME"))
			outline = shapely.wkt.loads(feature.GetGeometryRef().ExportToWkt())

			self.countries[country] = outline

	def max_polygon(self, polygons):
		max_polygon = polygons[0]
		# print(max_polygon)
		max_size = max_polygon.area
		for polygon in polygons:
			if polygon.area > max_size:
				max_size = polygon.area
				max_polygon = polygon
		return max_polygon

	def neighbors(self, country_code):
		found_neighbors = []
		outline = self.countries[country_code]

		for other_country in sorted(self.countries.keys()):
			if country_code == other_country: 
				continue

			other_outline = self.countries[other_country]
			if outline.touches(other_outline):
				found_neighbors.append(other_country)

		return found_neighbors

	def double_neighbors(self, country_code):
		found_neighbors = self.neighbors(country_code)
		found_2_neighbors = [country for neighbor in found_neighbors for country in self.neighbors(neighbor)]
		found_neighbors += found_2_neighbors
		found_neighbors = [c for c in found_neighbors if c != country_code]
		return list(set(found_neighbors))

	def distance_neighbors(self, country_code, distance):
		near_countries = []
		if country_code not in self.countries:
			return []
		country_pt = self.countries[country_code].centroid
		for other_country in self.countries:
			if other_country == country_code:
				continue

			near = False
			q_polygons = self.countries[other_country]
			polygons = q_polygons if type(q_polygons) == shapely.geometry.multipolygon.MultiPolygon else [q_polygons]
			max_polygon = self.max_polygon(polygons)

			polygons = [max_polygon] + [polygon for polygon in polygons if (polygon.area > 0.5) & (polygon != max_polygon)]
			distances = [self.geod.inv(country_pt.x, country_pt.y, polygon.centroid.x, polygon.centroid.y)[2] for polygon in polygons]
			min_idx = np.argmin(distances)
			closest_polygon = polygons[min_idx]

			for other_pt in closest_polygon.exterior.coords:
				pt_distance = self.geod.inv(country_pt.x, country_pt.y, other_pt[0], other_pt[1])[2] / 1000 / 1.6093
				if pt_distance < distance:
					near_countries.append(other_country)
					break

		return near_countries



	def gdp_neighbors(self, country_code):
		return set(self.double_neighbors(country_code)).union(set(self.distance_neighbors(country_code, 500)))








	def plot_shapely(self, ax, q_polygons, **kwargs):

		polygons = q_polygons if type(q_polygons) == shapely.geometry.multipolygon.MultiPolygon else [q_polygons]
		max_polygon = self.max_polygon(polygons)
		xs, ys = max_polygon.exterior.xy
		ax.fill(xs, ys, **kwargs)	

		for polygon in polygons:
			if polygon.area < 0.5 or polygon == max_polygon:
				continue
			xs, ys = polygon.exterior.xy
			ax.fill(xs, ys, **kwargs)

		return max_polygon

	def map_countries(self, country_code, gdp_neighbors, gdp_data, fy, independence_year, show=False, **kwargs):
		csr = gdp_data[(gdp_data.year >= fy) & (gdp_data.year <= independence_year)].corr().loc[country_code]
		fig, ax = plt.subplots(figsize=(15, 9))
		text_args = {
			'va': 'center', 
			'ha': 'center',
			'zorder': 20,
			'fontsize': 8,
			'bbox': dict(facecolor='white', edgecolor='white', boxstyle='round,pad=.03', alpha=0.7),
		}

		cgeo = self.countries[country_code]
		max_polygon = self.plot_shapely(ax, cgeo, fc='tomato', ec='k', zorder=5, **kwargs)
		ax.text(max_polygon.centroid.x, max_polygon.centroid.y, self.cc.code_to_name(country_code), **text_args)

		neighbors = self.gdp_neighbors(country_code)
		for neighbor_code in neighbors:
			if neighbor_code not in self.countries:
				continue

			ngeo = self.countries[neighbor_code]
			if neighbor_code in gdp_neighbors:
				max_polygon = self.plot_shapely(ax, ngeo, fc='royalblue', ec='k', zorder=2, alpha=0.5, **kwargs)
				text_str = '{}\nCorr: {:.2f}'.format(self.cc.code_to_name(neighbor_code), csr.loc[neighbor_code])
				ax.text(max_polygon.centroid.x, max_polygon.centroid.y, text_str, **text_args)
			else:
				max_polygon = self.plot_shapely(ax, ngeo, fc='gray', ec='k', zorder=2, alpha=0.5, **kwargs)
				ax.text(max_polygon.centroid.x, max_polygon.centroid.y, self.cc.code_to_name(neighbor_code), **text_args)
			
			
		ax.set_aspect('equal', 'datalim')

		if show:
			plt.show()
		return ax


# if __name__=="__main__": 
# 	n = Neighbors()
# 	n.map_countries('USA', ['CAN', 'MEX'])






