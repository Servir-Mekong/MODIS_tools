'''
This program should be run after the "fetch_historic_8d.py" and "calc_avg.py" programs
The program calculates the NDVI of permanent water bodies using the "permanent_water_mask.tif" file.
Next, 1/2 SD is added to the mean NDVI of permanent water, and pixels below this threshold are classified "surface water"
Surface water products are generated for the day of year (avg_SW) and the year of interest

Developer: Aakash Ahamed (aakash.ahamed@nasa.gov)
NASA Goddard Space Flight Center
Applied Sciences Program
Hydrological Sciences Laboratory 
'''

from osgeo import gdal
import numpy as np
import os
import shutil

# Functions 

def read_as_array(raster):
	ds = gdal.Open(raster)
	array = np.array(ds.GetRasterBand(1).ReadAsArray())
	return array

def chunks(r,n):
	for i in range (0,len(r),n):
		yield r[i:i+n]
		
def array2raster(rasterfn,newRasterfn,array):
	raster = gdal.Open(rasterfn)
	geotransform = raster.GetGeoTransform()
	proj = raster.GetProjection()
	originX = geotransform[0]
	originY = geotransform[3]
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]
	cols = raster.RasterXSize
	rows = raster.RasterYSize

	driver = gdal.GetDriverByName('GTiff')
	outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32) # Change dtype here
	outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
	outband = outRaster.GetRasterBand(1)
	outband.WriteArray(array)
		
	outRaster.SetProjection(proj)
	outband.FlushCache()
	
def main():

	# dirs
	dir = os.getcwd()
	ref_dir = os.path.join(dir,"reference_files")
	historic_dir = os.path.join(dir,"historic")
	results_dir = os.path.join(dir,"results")
	products_dir = os.path.join(dir,"products")

	# create "products" dir
	if os.path.exists(products_dir):
		pass
	else:
		os.mkdir(products_dir)
		
	# Read the water mask and remove nodata values
	water = [x for x in os.listdir(os.getcwd()) if "water" in x]
	water = "".join(water)
	water_array = read_as_array(water)
	water_array[water_array<0] = np.nan
	water_array[water_array>1] = np.nan

	# make a new object for the invert 
	inverted_water_array = [x for x in os.listdir(os.getcwd()) if "water" in x]
	inverted_water_array = "".join(inverted_water_array)
	inverted_water_array = read_as_array(inverted_water_array)
	
	water_array_fin = inverted_water_array
	water_array_fin[water_array_fin < 0] = np.nan
	water_array_fin[water_array_fin > 1] = np.nan
	water_array_fin[water_array_fin == 1] = np.nan
	water_array_fin[water_array_fin == 0] = 1
		
	# Read the images of flood years 
	# Terra
	comparison_img = [x for x in os.listdir(historic_dir) if "2013" in x and "Terra" in x]
	comparison_img = os.path.join(historic_dir,"".join(comparison_img))
	terra_comparison_arr = read_as_array(comparison_img)
	
	year = "".join(comparison_img)
	year = year[-20:-16]
	
	doy = "".join(comparison_img)
	doy = doy[-15:-12]
	
	# Aqua
	comparison_img = [x for x in os.listdir(historic_dir) if "2013" in x and "Aqua" in x]
	comparison_img = os.path.join(historic_dir,"".join(comparison_img))
	aqua_comparison_arr = read_as_array(comparison_img)
	
	year_average = np.nanmean(np.array([aqua_comparison_arr,terra_comparison_arr]), axis = 0)
	
	# Write the average NDVI for that year 
	outfile = "MCD_{}_ndvi_{}.tif".format(year,doy)
	outfile = os.path.join(products_dir,outfile)
	array2raster(comparison_img,outfile, year_average)
	
	# Read the average ndvi files and smash them together
	# Terra
	terra = [x for x in os.listdir(ref_dir) if "terra" in x]
	terra = (os.path.join(ref_dir,"".join(terra)))
	terra_average = read_as_array(terra)
		
	# aqua
	aqua = [x for x in os.listdir(ref_dir) if "aqua" in x]
	aqua = (os.path.join(ref_dir,"".join(aqua)))
	aqua_average = read_as_array(aqua)
		
	average = np.nanmean(np.array([aqua_average,terra_average]), axis = 0)
	
	outfile = '''MCD_avg_ndvi_{}.tif'''.format(doy)
	outfile = os.path.join(products_dir,outfile)
	array2raster(comparison_img,outfile, average)

	# CALCULATE THE NDVI VALUES OF PERMANENT WATER using the INVERTED water mask
	h2o = np.multiply(water_array_fin,average)
	h2o[h2o < -1.0] = np.nan
	h2o[h2o > 1.0] = np.nan
	print(np.nanmean(h2o))
	print("H20 average")
	
	f = h2o[np.logical_not(np.isnan(h2o))]
	
	outfile = 'average_h2o_vals.csv'
	outfile = os.path.join(results_dir,outfile)
	np.savetxt(outfile, f, delimiter = ",", fmt = "%f")
	
	outfile = "h2o.tif"
	outfile = os.path.join(results_dir,outfile)
	array2raster(comparison_img,outfile, h2o)
	
	# CLASSIFY SURFACE WATER BASED ON PERMANENT WATER NDVI AVGs for average image. Write out as water extetn
	half_sd = np.nanstd(h2o)/ 2
	upper_bound = np.nanmean(h2o) + half_sd
	print(upper_bound)
	print("upper bound")
	terrestrial_avg = np.multiply(average,water_array) # Mask ocean water	
	a = np.where(terrestrial_avg > upper_bound, 0, terrestrial_avg)
	a = np.where(water_array_fin == 1 , 1, a) # add the SW back in 
	a = np.ma.array(a, mask=np.isnan(a))
	a = np.where(a != 0, 1, a)

	array2raster(comparison_img, os.path.join(products_dir,'''avg_SW_{}.tif'''.format(doy)), a.data) # write 

	# Now do it for year of interest 
	year_average = np.nanmean(np.array([aqua_comparison_arr,terra_comparison_arr]), axis = 0) # smash the aqua and terra NDVIs together
	year_terrestrial = np.multiply(year_average,water_array) #  mask ocean water 
	x = np.where(year_terrestrial > upper_bound, 0, year_terrestrial)
	x = np.where(water_array_fin == 1 , 1, x) # add the SW back in 
	y = np.ma.array(x, mask=np.isnan(x))
	y = np.where(y != 0, 1, y)
		
	outfile = '''{}_SW_{}.tif'''.format(year,doy)
	outfile = os.path.join(products_dir,outfile)
	array2raster(comparison_img, outfile, y.data)
	
	
if __name__ == '__main__':
	main()
