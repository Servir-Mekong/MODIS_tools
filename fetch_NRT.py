'''
This program does the following:
1) Fetch the most recent MOD09GQ Surface Reflectance files from LANCE
2) Fetch the most recent MOD35_L2 cloud mask files from LANCE
3) Calculates NDVI, masks clouds, resamples, clips to SE Asia region. 
4) Deposits files in the "daily_files_nrt" directory

You must sign up for an account and supply credentials at L 63
If this script ever breaks for inexplicable reasons, try changing to nrt4, nrt2, nrt1 in previous line

Developer: Aakash Ahamed (aakash.ahamed@nasa.gov)
NASA Goddard Space Flight Center
Applied Sciences Program
Hydrological Sciences Laboratory 
'''

#############
# Libraries #
#############

import os
import ftplib
import datetime
from time import gmtime, strftime
from osgeo import gdal
import numpy as np
from numpy import *
import shutil
import argparse

#############
# Functions #
#############

class NDVI_NRT():
	
	def __init__ (self):
		pass
	
	def get_current_directory(self):
		return(os.path.dirname(os.path.realpath(__file__)))	
				
	def get_dates(self):
		year = str(datetime.datetime.now().year)
		month = str(datetime.datetime.now().month)
		day = str(datetime.datetime.now().day)
		day_of_year = datetime.datetime.now().timetuple().tm_yday
		hour = str(datetime.datetime.now().hour)
		minute = str(datetime.datetime.now().minute)
		doy = str(day_of_year - 0) # <-- Change this to go back x days 
		
		date= {
		"year" : year,
		"month": month,
		"day" : day,
		"doy" : doy,
		"hour" : hour,
		"minute" : minute}
			
		return date

	def ftp_protocols(self):
		ftp = ftplib.FTP('nrt3.modaps.eosdis.nasa.gov')
		ftp.login("username","password") # ENTER USERNAME HERE
		ftp.cwd('orders/flood_watch')
		ftp.set_pasv(True)
		return(ftp)
	
	def get_latest_files(self, files_list, product):
	
		latest_list = []
		for file in files_list:
			latest_files = []
			for i in file:
				if i.startswith(product):
					latest_files.append(i[24:36]) # This is the time string in the files 
			latest_list.append(latest_files)
					
		latest_images = []
		for list in latest_list:
			list.sort() 
			latest_images.append(list[-1])

		final_prods = []
		for image in latest_images:
			for list in files_list:
				for item in list:
					if image in item:
						final_prods.append(item)

		return(final_prods)
		
		
	def get_latest_cloud_files(self, files_list, product):
	
		latest_list = []
		for file in files_list:
			latest_files = []
			for i in file:
				if i.startswith(product):
					latest_files.append(i[27:39]) # Same as above function 
			latest_list.append(latest_files)
					
		latest_images = []
		for list in latest_list:
			list.sort() 
			latest_images.append(list[-1])

		final_prods = []
		for image in latest_images:
			for list in files_list:
				for item in list:
					if image in item:
						final_prods.append(item)

		return(final_prods)
		
	def build_vrt_string(self, dir, band, product):
		files = os.listdir(dir)
		list = [x for x in files if x.endswith(band+".tif") and product in x]
		string = " ".join(list)
		return string

	def build_vrt_table(self, file_string):
		outfile = file_string[0:14] + "_" + file_string[-10:-4]
		command = '''gdalbuildvrt {}.vrt {}'''.format(outfile, file_string)
		os.system(command)
		return [x for x in os.listdir(os.getcwd()) if file_string[0:16] in x]
		
	def build_cloud_vrt_table(self, file_string):
		outfile = file_string[0:17] + "_" + file_string[-6:-4]
		command = '''gdalbuildvrt {}.vrt {}'''.format(outfile, file_string)
		os.system(command)
		return [x for x in os.listdir(os.getcwd()) if file_string[0:16] in x]
		
	def mosaic(self, vrt):
		mosaic_cmd = '''gdalwarp {} {}.tif -t_srs EPSG:4326 -dstnodata "-999" -overwrite'''.format(vrt,vrt[0:-4])
		os.system(mosaic_cmd)
	
	def get_raster_extent(self, raster):
		r = gdal.Open(raster)
		geoTransform = r.GetGeoTransform()
		ulx = geoTransform[0]
		uly = geoTransform[3]
		lrx = ulx + geoTransform[1]*r.RasterXSize
		lry =  uly + geoTransform[5]*r.RasterYSize
		pixelX=geoTransform[1]
		pixelY=geoTransform[5]
		extent = [ulx,uly,lrx,lry] 
		del geoTransform

		return(extent)
	
	def clip(self, extent, raster, outfiles):
		str_extent = ' '.join(map(str,extent))
		command = '''gdal_translate -projwin {} {} {}'''.format(str_extent, raster, outfiles)
		os.system(command)
		return(outfiles)
		
	def read_as_array(self, raster):
		ds = gdal.Open(raster)
		flood_array = np.array(ds.GetRasterBand(1).ReadAsArray())
		return flood_array

	def calc_NDVI(self, list, sat): 
		if sat == "MOD": 
			b1 = [x for x in list if x.startswith("MOD09") if "Band_1" in x]
			b1 = "".join(b1)			
			b2 = [x for x in list if x.startswith("MOD09") if "Band_2" in x]
			b2 = "".join(b2)

			# Read as array 
			g = gdal.Open(b1)
			red = g.ReadAsArray()

			g = gdal.Open(b2)
			nir = g.ReadAsArray()

			red = array(red, dtype = float32)
			nir = array(nir, dtype = float32)

			check = np.logical_and ( red > 0, nir > 0 )
			ndvi = np.where ( check,  (nir - red ) / ( nir + red ), -999)
			
			ndvi[ndvi>1] = np.nan #0
			ndvi[ndvi < -1] = np.nan #0
			
			return(ndvi)
			
		if sat == "MYD":
			b1 = [x for x in list if x.startswith("MYD09") if "Band_1" in x]
			b1 = "".join(b1)
			b2 = [x for x in list if x.startswith("MYD09")  if "Band_2" in x]
			b2 = "".join(b2)
			
			# Read as array 
			g = gdal.Open(b1)
			red = g.ReadAsArray()
			
			red.astype(float32)

			g = gdal.Open(b2)
			nir = g.ReadAsArray()
			
			nir.astype(float32)
			
			check = np.logical_and ( red > 0, nir > 0 )
			ndvi = np.where ( check,  (nir - red ) / ( nir + red ), -999)
			
			ndvi[ndvi>1] = np.nan #0
			ndvi[ndvi < -1] = np.nan #0
			
			return(ndvi)
			
	def write_raster(self, array, raster, dates, satellite):
		match_raster = gdal.Open(raster)
		cols = array.shape[0]
		rows = array.shape[1]
		trans=match_raster.GetGeoTransform()
		proj=match_raster.GetProjection()
		# nodatav = array.GetNoDataValue()
		
		outfile =  satellite+dates['year']+dates['doy']+".tif"
		
		# Create the file, using the information from the original file
		outdriver = gdal.GetDriverByName("GTiff")
		outdata = outdriver.Create(str(outfile), rows, cols, 1, gdal.GDT_Float32)

		# Write the array to the file
		outdata.GetRasterBand(1).WriteArray(array)
		
		# Georeference the image
		outdata.SetGeoTransform(trans)

		# Write projection information
		outdata.SetProjection(proj)
		return outfile
		
	def array2raster(self, rasterfn,newRasterfn,array):
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
		
	def build_cloud_tiff(self, list):
		command = '''gdal_translate -ot "Byte" -of "Gtiff" -b "5" {} {}.tif '''.format("".join(list),"".join(list[0:26]) + "_b5")
		print(command)
		os.system(command)
		
	def read_bits (self, sat):
		mask = [x for x in os.listdir(os.getcwd()) if x.startswith(sat+"35") if x.endswith("b5.tif")]
		mask = "".join(mask)
		mask = os.path.join(os.getcwd(),mask)
		mask = gdal.Open(mask)
		mask = mask.ReadAsArray()
		
		# read the bits 
		bitlist = []

		for x in np.nditer(mask):
			bitlist.append(bin(x)[2:].zfill(8)) # change to zpad different lengths 
		
		cloudmask = []

		for item in bitlist:
			if item[5:7] == "00" or item[5:7] == "01" or item[5:7] == "10":
				cloudmask.append(0)
			else:
				cloudmask.append(1)

		cloud_array = np.array(cloudmask).reshape(mask.shape[0],mask.shape[1])
		return(cloud_array)
		
	def get_raster_resolution(self, raster):
		r = gdal.Open(raster)
		geoTransform = r.GetGeoTransform()
		ulx = geoTransform[0]
		uly = geoTransform[3]
		lrx = ulx + geoTransform[1]*r.RasterXSize
		lry =  uly + geoTransform[5]*r.RasterYSize
		pixelX=geoTransform[1]
		pixelY=geoTransform[5]
		extent = [ulx,uly,lrx,lry] 
		del geoTransform
		return(pixelX,pixelY)
		
	def resample(self, resolution, raster, outfiles):
		str_resolution = ' '.join(map(str,resolution))
		command = '''gdalwarp -tr {} {} {} -overwrite'''.format(str_resolution, raster, outfiles)
		print(command)
		os.system(command)
		return(outfiles)

	def cleanup(self, directory):
		files = [x for x in os.listdir(directory) if "water_mask" not in x and x.endswith(".tif") or x.endswith(".hdf")]
		for i in files:
			os.remove(i)
		vrts = [x for x in os.listdir(os.getcwd()) if x.endswith(".vrt")]
		for vrt in vrts:
			os.remove(vrt)
			
########
# MAIN #
########

def main():

	#######################
	# NRT NDVI (M*D09) GQ #
	#######################
	
	# Make class instance
	image = NDVI_NRT()

	# Connect to FTP and cd to relevant dir 
	ftp = image.ftp_protocols()
	dates = image.get_dates()
	starttime = strftime("%Y-%m-%d %H:%M:%S", gmtime()) # print start time for logs 
	print(starttime)

	tile_dict = {
		"1" : "*100E010N*",
		"2" : "*100E020N*",
		"3" : "*100E030N*",
		"4" : "*090E010N*",
		"5" : "*090E020N*",
		"6" : "*090E030N*",
		"7" : "*110E010N*",
		"8" : "*110E020N*",
		"9" : "*110E030N*"
		}

	files_list = []	

	# Make list of lists containing files for each tile
	for k,v in tile_dict.items():
		files_list.append(ftp.nlst("*"+dates['year']+dates['doy']+v))

	# Filter those lists for the latest files
	MOD09_prods = image.get_latest_files(files_list, "MOD09")
	MOD09_prods = [x for x in MOD09_prods if "Band_7" not in x]

	for i in MOD09_prods:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)

	MYD09_prods = image.get_latest_files(files_list, "MYD09")
	MYD09_prods = [x for x in MYD09_prods if "Band_7" not in x]

	for i in MYD09_prods:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)
		
	# Build VRTs for b1 and b2 and Mosaic
	MYD_b1 = image.build_vrt_string(os.getcwd(),"Band_1","MYD")
	MYD_b2 = image.build_vrt_string(os.getcwd(),"Band_2","MYD")

	image.build_vrt_table(MYD_b1)
	image.build_vrt_table(MYD_b2)

	MOD_b1 = image.build_vrt_string(os.getcwd(),"Band_1","MOD")
	MOD_b2 = image.build_vrt_string(os.getcwd(),"Band_2","MOD")

	image.build_vrt_table(MOD_b1)
	image.build_vrt_table(MOD_b2)

	vrt_files = [x for x in os.listdir(image.get_current_directory()) if x.endswith(".vrt")]
	for file in vrt_files:
		image.mosaic(file)
		
	# Clean up
	trash = [x for x in os.listdir(os.getcwd()) if "250m" in x or x.endswith(".vrt")]
	for i in trash:
		os.remove(i)

	# Calculate NDVI 
	MYD_mosaics = [x for x in os.listdir(os.getcwd()) if x.startswith("MYD09")] #x.endswith("clip.tif")]	
	MYD_ndvi = image.calc_NDVI(MYD_mosaics, "MYD")

	MOD_mosaics = [x for x in os.listdir(os.getcwd()) if x.startswith("MOD09")] #x.endswith("clip.tif")]	
	MOD_ndvi = image.calc_NDVI(MOD_mosaics, "MOD")

	###############
	# Cloud Masks #
	###############
	
	# Grab the cloud masks
	MOD35_prods = image.get_latest_cloud_files(files_list, "MOD35_L2")

	for i in MOD35_prods:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)
		
	MYD35_prods = image.get_latest_cloud_files(files_list, "MYD35_L2")

	for i in MYD35_prods:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)
	
	# Convert the HDF files to tiffs
	cloud_files_aqua = [x for x in os.listdir(os.getcwd()) if "MYD35_L2" in x]
	for file in cloud_files_aqua:
		image.build_cloud_tiff(file)

	cloud_files_terra = [x for x in os.listdir(os.getcwd()) if "MOD35_L2" in x]	
	for file in cloud_files_terra:
		image.build_cloud_tiff(file)
		
	# Mosaic
	MYD_35 = image.build_vrt_string(os.getcwd(), "b5", "MYD35_L2")
	MOD_35 = image.build_vrt_string(os.getcwd(), "b5", "MOD35_L2")

	image.build_cloud_vrt_table(MYD_35)
	image.build_cloud_vrt_table(MOD_35)

	cloud_vrt = [x for x in os.listdir(os.getcwd()) if "35" in x if x.endswith(".vrt")]
	for i in cloud_vrt:
		image.mosaic(i)
	
	# Clean up 
	cloud_trash = [x for x in os.listdir(os.getcwd()) if "N" in x if "E" in x]
	for i in cloud_trash:
		os.remove(i)
	
	mod_mask = image.read_bits("MOD")
	myd_mask = image.read_bits("MYD")

	mod_mask = np.repeat(mod_mask, 4, axis = 0)
	mod_mask = np.repeat(mod_mask, 4, axis = 1)

	myd_mask = np.repeat(myd_mask, 4, axis = 0)
	myd_mask = np.repeat(myd_mask, 4, axis = 1)
	
	# Apply cloud masks
	mod_product = np.multiply(MOD_ndvi,mod_mask) 
	mod_product[mod_product == 0] = np.nan

	myd_product = np.multiply(MYD_ndvi,myd_mask)
	myd_product[myd_product == 0] = np.nan

	match_raster = [x for x in os.listdir(os.getcwd()) if "MOD09" in x if "Band_1" in x if "E" not in x if "vrt" not in x]
	match_raster = "".join(match_raster)
	match_raster = os.path.join(os.getcwd(),match_raster) 
	
	mod_tiff = image.write_raster(mod_product,match_raster,dates,"MOD_NDVI_temp")
	myd_tiff = image.write_raster(myd_product,match_raster,dates,"MYD_NDVI_temp")
	
	# Resample to the water mask resolution 
	reference_file = [x for x in os.listdir(os.getcwd()) if "water_mask" in x]
	reference_file = "".join(reference_file)
	reference_file = os.path.join(os.getcwd(),reference_file)
	resolution = image.get_raster_resolution(reference_file)
	
	final_prod_mod = image.resample(resolution, mod_tiff, "MOD_resampled.tif")
	final_prod_myd = image.resample(resolution, myd_tiff, "MYD_resampled.tif")
	
	# Clip to water mask extent
	extent = image.get_raster_extent(reference_file)
	image.clip(extent, os.path.join(os.getcwd(),"MOD_resampled.tif"), "MOD_NDVI_"+dates['year']+dates['doy']+".tif")
	image.clip(extent, os.path.join(os.getcwd(),"MYD_resampled.tif"), "MYD_NDVI_"+dates['year']+dates['doy']+".tif")
	
	# Select the clipped products 
	final_prods = [x for x in os.listdir(os.getcwd()) if dates['year'] in x and "temp" not in x and "Band" not in x and "35" not in x]
	
	# Create the daily NRT dir 
	if os.path.exists(os.path.join(os.getcwd(),"daily_files_nrt")):
		pass
	else:
		os.mkdir(os.path.join(os.getcwd(),"daily_files_nrt"))
	daily_dir = os.path.join(os.getcwd(),"daily_files_nrt")

	# Copy over the final prods 
	for prod in final_prods:
		shutil.copy(prod,daily_dir)

	#Cleanup
	garbage = [x for x in os.listdir(os.getcwd()) if "MOD" in x or "MYD" in x]
	for i in garbage:
		os.remove(i)
	
if __name__ == '__main__':
	main()