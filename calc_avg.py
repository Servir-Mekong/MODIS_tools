'''
This program should be run after the "fetch_historic_8d.py" program.
The program loops through a directory of NDVI files ("historic"), calculates the average per-pixel NDVI,
and stores one .tif file each for Aqua and Terra in a directory called "reference_files"

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
	
def calc_average(chunk_list):

	averages = []
	variances = []

	for chunk in chunk_list:
		arrays = []
		temp = chunk
		for i in temp:
			arrays.append(read_as_array(i))
			print("reading...." + i)
		
		for array in arrays: # mask
			array[array == 0] = np.nan
		
		time_array = np.array(arrays)
		average = np.nanmean(time_array,axis = 0)
		averages.append(average)
		variance= np.nanvar(time_array,axis = 0)
		variances.append(variance)

	fin_array = np.array(averages)
	fin_variance = np.array(variances)

	true_average = np.nanmean(fin_array, axis = 0)
	
	return (true_average)
	

def main():
	dir = os.getcwd()
	
	# make dir for averages
	ref_dir = os.path.join(dir,"reference_files")
	
	if os.path.exists(ref_dir):
		pass
	else: 
		os.mkdir(ref_dir)
		
	results_dir = os.path.join(dir,"results")

	if os.path.exists(results_dir):
		pass
	else: 
		os.mkdir(results_dir)
		
	# Set dir with historic files (created automatically from fetch_8d scripts)
	historic_dir = os.path.join(dir,"historic")
	os.chdir(historic_dir)

	# Aqua 
	files = [x for x in os.listdir(os.getcwd()) if x.endswith(".tif") and "Aqua" in x] #and "2011" not in x]

	chunk_list = list(chunks(files,2))

	aqua_average = calc_average(chunk_list)
	
	del chunk_list
	
	match_raster = "".join(files[0])
	doy = match_raster[10:13]
	
	array2raster(match_raster, doy+"_aqua_average.tif", aqua_average)
	
	del aqua_average

	file = [x for x in os.listdir(os.getcwd()) if x.endswith("average.tif")]
	file = "".join(file)

	shutil.copy(file,ref_dir)
	os.remove(file)
	
	# Terra
	files = [x for x in os.listdir(os.getcwd()) if x.endswith(".tif") and "Terra" in x]# and "2011" not in x]

	chunk_list = list(chunks(files,2))
	
	terra_average = calc_average(chunk_list)
	
	match_raster = "".join(files[0])
	doy = match_raster[11:14]
	
	array2raster(match_raster, doy+"_terra_average.tif", terra_average)
	
	file = [x for x in os.listdir(os.getcwd()) if x.endswith("average.tif")]
	file = "".join(file)

	shutil.copy(file,ref_dir)
	os.remove(file)
	
if __name__ == '__main__':
	main()