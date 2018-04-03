'''
This file fetches the historic MOD09 and MYD09 (Q1 or GQ) surface reflectance data for the length of the MODIS record. 
Clouds are masked using the A1 bit fields for the Q1 product and the GA bit fields for the GQ product.
NDVI is calculated and files are stored in a directory called "historic"
To select specific products (Q1, GQ) and / or years (2000 - ), change the "product" and "doy" variables in the params section 

Developer: Aakash Ahamed (aakash.ahamed@nasa.gov)
NASA Goddard Space Flight Center
Applied Sciences Program
Hydrological Sciences Laboratory 
'''

import os
import ftplib
import datetime
from osgeo import gdal
import numpy as np
from numpy import *
import shutil
import itertools

###########################################
# Small function to grab the current date #
###########################################

def get_dates ():
        year = str(datetime.datetime.now().year)
        month = str(datetime.datetime.now().month)
        day = str(datetime.datetime.now().day)
        day_of_year = datetime.datetime.now().timetuple().tm_yday
        doy = str(day_of_year - 0) # set this to go back x number of days

        date= {
        "year" : year,
        "month": month,
        "day" : day,
        "doy" : doy}

        return date

dates = get_dates()

#############################
########## Params ###########
#############################

# Set these to the dates and products you wish to analyze 
product = "Q1"
doy = dates['doy']

years = [x for x in range(2000,2017)] # does not look for 2016 data

aqua_years = [x for x in range(2003,2017)] # does not look for 2002 historic data 

#############################
########## Classes ##########
#############################

class NDVI_8day():

	def __init__ (self, year, doy, sat, product):
		self.year = year
		self.doy = doy
		self.sat = sat
		self.product = product
		
	def get_current_directory(self):
		return(os.path.dirname(os.path.realpath(__file__)))	
			
	def get_dates (self): # Leave blank to select the current day or alter the dict to pick a date
		year = str(datetime.datetime.now().year)
		month = str(datetime.datetime.now().month)
		day = str(datetime.datetime.now().day)
		day_of_year = datetime.datetime.now().timetuple().tm_yday
		doy = str(day_of_year - 0) 
		
		date= {
		"year" : year,
		"month": month,
		"day" : day,
		"doy" : doy}
			
		return date

	def ftp_protocols(self,sat,product):
		print("Performing analysis for " + sat + " " + product)
		ftp = ftplib.FTP('ladsweb.nascom.nasa.gov') # https://ladsweb.nascom.nasa.gov
		ftp.login("anonymous","anonymous")
		if sat == "Terra" and product == "Q1":  
			ftp.cwd('allData/6/MOD09Q1')
		elif sat == "Terra" and product == "A1":
			ftp.cwd('allData/6/MOD09A1')
		elif sat == "Aqua" and product == "Q1":
			ftp.cwd('allData/6/MYD09Q1')
		elif sat == "Aqua" and product == "A1":
			ftp.cwd('allData/6/MYD09A1')
			
		else: 
			print("COULDN'T LOCATE MODIS PRODUCT, TRY AGAIN")
		ftp.set_pasv(True)
		
		return(ftp)

	def get_nearest_date(self, base, dates):
		nearness = { abs(base.timestamp() - date.timestamp()) : date for date in dates }
		return nearness[min(nearness.keys())]	
		
	def get_combo_list(self):
		h = ['h26','h27','h28'] # Change to go to a different area of the world
		v = ['v06','v07','v08'] # These too 
		combinations = list(itertools.product(*[h,v]))
		combo_strings = []
		for combo in combinations:
			combo_strings.append(''.join(combo))
		return combo_strings
	
	def hdf_to_gtiff_8day(self, file, band):
		if "Q1" in file and band == "b01" or "b02":
			command = '''gdal_translate -of Gtiff HDF4_EOS:EOS_GRID:"{}":MOD_Grid_250m_Surface_Reflectance:sur_refl_{} {}_{}.tif '''.format(file,band,file[0:23],band)
			print(command)
			os.system(command)
		# GQ command 
	
	def build_vrt_string(self, dir, band):
		files = os.listdir(dir)
		list = [x for x in files if x.endswith(band+".tif")]
		string = " ".join(list)
		return string

	def build_vrt_table(self, file_string):
		outfile = file_string[0:16] + "_" + file_string[-7:-4]
		command = '''gdalbuildvrt {}.vrt {}'''.format(outfile, file_string)
		os.system(command)
		return [x for x in os.listdir(os.getcwd()) if file_string[0:16] in x]

	def mosaic(self, vrt):
		mosaic_cmd = '''gdalwarp {} {}.tif -t_srs EPSG:4326 -dstnodata "-999" -overwrite'''.format(vrt,vrt[0:-4])
		os.system(mosaic_cmd)

	def cleanup_hdf(self, dir, combo):
		files = [x for x in os.listdir(os.getcwd()) if x.endswith(".hdf")]
		for i in files:
			os.remove (i)
		files = [x for x in os.listdir(os.getcwd()) if combo in x]
		for i in files:
			os.remove(i)
		files = [x for x in os.listdir(os.getcwd()) if ".vrt" in x]
		for i in files:
			os.remove(i)
	
	def read_bits(self, array):
		bitlist = []
		
		for i in np.nditer(array):
			bitlist.append(bin(i)[2:].zfill(16))
		
		cloudmask = []
		
		for item in bitlist:
			if item[-4:-2] == "01" or item[-4:-2] == "10" :# These are the MODIS bit fields 
				cloudmask.append(0)
			else:
				cloudmask.append(1)
		
		cloud_array = np.array(cloudmask).reshape(array.shape[0],array.shape[1])
		return(cloud_array)		
		
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
	
	def clip(self, extent, raster, outfiles):
		str_extent = ' '.join(map(str,extent))
		command = '''gdal_translate -projwin {} {} {}'''.format(str_extent, raster, outfiles)
		os.system(command)
		return(outfiles)
	
	def resample(self, resolution, raster, outfiles):
		str_resolution = ' '.join(map(str,resolution))
		command = '''gdalwarp -tr {} {} {}'''.format(str_resolution, raster, outfiles)
		print(command)
		os.system(command)

	def read_as_array(self, raster):
		ds = gdal.Open(raster)
		array = np.array(ds.GetRasterBand(1).ReadAsArray())
		return array
		
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
	

class CloudMask8day():

	def __init__ (self, year, doy, sat, product):
		self.year = year
		self.doy = doy
		self.sat = sat
		self.product = product
		
	def get_current_directory(self):
		return(os.path.dirname(os.path.realpath(__file__)))	
			
	def get_dates (self): # Leave blank to select the current day or alter the dict to pick a date
		year = str(datetime.datetime.now().year)
		month = str(datetime.datetime.now().month)
		day = str(datetime.datetime.now().day)
		day_of_year = datetime.datetime.now().timetuple().tm_yday
		doy = str(day_of_year - 0) #str(doy) #
		
		date= {
		"year" : year,
		"month": month,
		"day" : day,
		"doy" : doy}
			
		return date

	def ftp_protocols(self,sat,product):
		print("Performing analysis for " + sat + " " + product)
		ftp = ftplib.FTP('ladsweb.nascom.nasa.gov') # https://ladsweb.nascom.nasa.gov
		ftp.login("anonymous","anonymous")
		if sat == "Terra" and product == "Q1":  
			ftp.cwd('allData/6/MOD09Q1')
		elif sat == "Terra" and product == "A1":
			ftp.cwd('allData/6/MOD09A1')
		elif sat == "Aqua" and product == "Q1":
			ftp.cwd('allData/6/MYD09Q1')
		elif sat == "Aqua" and product == "A1":
			ftp.cwd('allData/6/MYD09A1')
			
		else: 
			print("COULDN'T LOCATE MODIS PRODUCT, TRY AGAIN")
		ftp.set_pasv(True)
		
		return(ftp)

	def get_combo_list(self):
		h = ['h26','h27','h28'] # Change to go to a different area of the world
		v = ['v06','v07','v08'] # These too 
		combinations = list(itertools.product(*[h,v]))
		combo_strings = []
		for combo in combinations:
			combo_strings.append(''.join(combo))
		return combo_strings
	
	def hdf_to_gtiff(self, file, band):
		if "A1" in file:
			command = '''gdal_translate -of Gtiff HDF4_EOS:EOS_GRID:"{}":MOD_Grid_500m_Surface_Reflectance:sur_refl_state_500m {}_{}.tif '''.format(file,file[0:23],"bqa")
			os.system(command)

		return([x for x in os.listdir(os.getcwd()) if x.endswith((band)+".tif")])
	
	def build_vrt_string(self, dir, band):
		files = os.listdir(dir)
		list = [x for x in files if x.endswith(band+".tif")]
		string = " ".join(list)
		return string

	def build_vrt_table(self, file_string):
		outfile = file_string[0:16] + "_" + file_string[-7:-4]
		command = '''gdalbuildvrt {}.vrt {}'''.format(outfile, file_string)
		os.system(command)
		return [x for x in os.listdir(os.getcwd()) if file_string[0:16] in x]

	def mosaic(self, vrt):
		mosaic_cmd = '''gdalwarp {} {}.tif -t_srs EPSG:4326 -dstnodata "-999" -overwrite'''.format(vrt,vrt[0:-4])
		os.system(mosaic_cmd)
		
	def read_bits(self, array):
		bitlist = []
		
		for i in np.nditer(array):
			bitlist.append(bin(i)[2:].zfill(16))
		
		cloudmask = []
		
		for item in bitlist:
			if item[-2:] == "01" or item[-2:] == "10" or item[-3] == "1":# These are the MODIS bit fields 
				cloudmask.append(0)
			else:
				cloudmask.append(1)
		
		cloud_array = np.array(cloudmask).reshape(array.shape[0],array.shape[1])
		return(cloud_array)	
		
	def cleanup_hdf(self, dir, combo):
		files = [x for x in os.listdir(os.getcwd()) if x.endswith(".hdf")]
		for i in files:
			os.remove (i)
		files = [x for x in os.listdir(os.getcwd()) if combo in x]
		for i in files:
			os.remove(i)
		files = [x for x in os.listdir(os.getcwd()) if ".vrt" in x]
		for i in files:
			os.remove(i)
			
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
		
	def clip(self, extent, raster, outfiles):
		str_extent = ' '.join(map(str,extent))
		command = '''gdal_translate -projwin {} {} {}'''.format(str_extent, raster, outfiles)
		os.system(command)
		return(outfiles)
		
	def resample(self, resolution, raster, outfiles):
		str_resolution = ' '.join(map(str,resolution))
		command = '''gdalwarp -tr {} {} {}'''.format(str_resolution, raster, outfiles)
		print(command)
		os.system(command)

	def read_as_array(self, raster):
		ds = gdal.Open(raster)
		array = np.array(ds.GetRasterBand(1).ReadAsArray())
		return array
		
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
		
#############################
########### Main ############
#############################

def ndvi_8day(product,year,doy,satellite):
	image = NDVI_8day(year, doy , satellite, product)
	year = image.year
	doy = image.doy
	sat = image.sat
	prod = image.product

	print("YEAR IS : " + str(year))
	print("DOY IS : " + str(doy))
	print("SATELLITE IS : " + str(sat))
	print("PRODUCT IS : " + image.product)

	ftp = image.ftp_protocols(image.sat,image.product)
	files = ftp.nlst()

	dir = ftp.cwd("".join(image.year)) # cd to year 
	files = ftp.nlst()
	
	dates = [int(f) for f in files] # Grab the nearest window to input date 
	
	nearest = min(dates, key=lambda x:abs(x-int(image.doy)))
	
	nearest = str(nearest)
	
	# Make sure "nearest" is 3 characters long 
	nearest = nearest.zfill(3)
		
	dir = ftp.cwd("".join(nearest)) # cd to day 
	
	files = ftp.nlst()

	combos = image.get_combo_list() # get list of tiles 
	
	files_of_interest = []

	for file in files:
		for combo in combos:
			if "".join(combo) in "".join(file):
				files_of_interest.append(file)
			
	if len(files_of_interest) !=9:
		print('''TILES OF INTEREST DO NOT EXIST FOR {}... MOVING TO NEXT YEAR'''.format(str(year)))
		return
				
	# Download
	for i in files_of_interest:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)

	files = [x for x in os.listdir(image.get_current_directory()) if x.endswith('.hdf')]
	print(files)

	for file in files:
		print(file[0:23])
		image.hdf_to_gtiff_8day("".join(file),"b01")
		image.hdf_to_gtiff_8day("".join(file),"b02")
		image.hdf_to_gtiff_8day("".join(file),"bqa")

	# Now you have the .tiff files. Mosaic them	
	b1 = image.build_vrt_string(image.get_current_directory(),"b01")
	b2 = image.build_vrt_string(image.get_current_directory(),"b02")

	image.build_vrt_table(b1)
	image.build_vrt_table(b2)

	vrt_files = [x for x in os.listdir(image.get_current_directory()) if x.endswith(".vrt")]

	for file in vrt_files:
		image.mosaic(file)
		
	for combo in combos:
		image.cleanup_hdf(os.listdir(image.get_current_directory()),"".join(combo))

	# Clip
	reference_file = [x for x in os.listdir(image.get_current_directory()) if "water_mask.tif" in x]
	reference_file = "".join(reference_file)

	extent = image.get_raster_extent(reference_file)
	# resolution = image.get_raster_resolution(reference_file)

	# Clip
	files = [x for x in os.listdir(image.get_current_directory()) if nearest in x]
	for i in files:
		image.clip (extent, "".join(i), "".join(i[:-4])+"_clip.tif") 

	# Cleanup 
	files = [x for x in os.listdir(image.get_current_directory()) if not os.path.isdir(x) and nearest in x and "clip" not in x]
	for i in files:
		os.remove(i)

	# Calculate NDVI with checks for funky data 
	b2 = "".join([x for x in os.listdir(image.get_current_directory()) if "b02" in x])
	b2 = os.path.join(image.get_current_directory(),b2)
	b1 = "".join([x for x in os.listdir(image.get_current_directory()) if "b01" in x])
	b1 = os.path.join(image.get_current_directory(),b1)

	g = gdal.Open(b1)
	red = g.ReadAsArray()	
	red.astype(float32)

	g = gdal.Open(b2)
	nir = g.ReadAsArray()
	nir.astype(float32)
			
	check = np.logical_and ( red > 0, nir > 0 )
	ndvi = np.where ( check,  (nir - red ) / ( nir + red ), -999)
			
	ndvi[ndvi>1] = 0
	ndvi[ndvi < -1] =0

	print("##################                NDVI                  ##############")

	print(ndvi)
	print(ndvi.shape)
	print(np.min(ndvi))
	print(np.max(ndvi))
	print(np.mean(ndvi))

	
	''' Cloud '''

	year = image.year
	doy = image.doy
	sat = image.sat
	if image.product == "GQ":
		product = "GA"
	if image.product == "Q1":
		product = "A1"
		
	print(year,doy,sat,product)
		
	mask = CloudMask8day(year,doy,sat,product)

	ftp = mask.ftp_protocols(sat,product)
	files = ftp.nlst()

	dir = ftp.cwd("".join(mask.year)) # Cd to a year
	files = ftp.nlst()

	dir = ftp.cwd("".join(nearest)) # Cd to a day
	files = ftp.nlst()

	combos = image.get_combo_list() # get relevant path/row combinations 

	files_of_interest = []

	for file in files:
		for combo in combos:
			if "".join(combo) in "".join(file):
				files_of_interest.append(file)
				
	# Download
	for i in files_of_interest:
		ftp.retrbinary('RETR %s' % i, open(i,"wb").write)

	files = [x for x in os.listdir(image.get_current_directory()) if x.endswith('.hdf')]

	for file in files:
		mask.hdf_to_gtiff("".join(file),"bqc")

	# Mosaic 
	qc = mask.build_vrt_string(mask.get_current_directory(),"bqa")

	mask.build_vrt_table(qc)

	vrt_files = [x for x in os.listdir(mask.get_current_directory()) if x.endswith(".vrt")]

	for file in vrt_files:
		mask.mosaic(file)

	for combo in combos:
		mask.cleanup_hdf(os.listdir(mask.get_current_directory()),"".join(combo))
		
	# Clip
	reference_file = [x for x in os.listdir(mask.get_current_directory()) if "water_mask.tif" in x] # "b01" in x] # dryBaseSEA.tif reference.tif 
	reference_file = "".join(reference_file)

	extent = mask.get_raster_extent(reference_file)

	# Clip
	files = [x for x in os.listdir(mask.get_current_directory()) if nearest and "bqa" in x]
	for i in files:
		mask.clip (extent, "".join(i), "".join(i[:-4])+"_clip.tif") 

	# Cleanup 
	files = [x for x in os.listdir(mask.get_current_directory()) if not os.path.isdir(x) and nearest in x and "clip" not in x]
	print(files)
	for i in files:
		os.remove(i)

	# Build cloud mask
	file = [x for x in os.listdir(mask.get_current_directory()) if "clip" in x and "bqa" in x]
	file = "".join(file)
	file = os.path.join(mask.get_current_directory(),file)
	print(file)

	cloud_array = mask.read_as_array(file)
	cloud_array = mask.read_bits(cloud_array)

	cloud_array = np.repeat(cloud_array,2, axis = 0)
	cloud_array = np.repeat(cloud_array,2, axis = 1)

	# Apply cloud mask
	almost_fin_array = np.multiply(cloud_array,ndvi)

	# Apply water mask
	water_array = [x for x in os.listdir(os.getcwd()) if "water" in x]
	water_array = "".join(water_array)
	water_array = image.read_as_array(water_array)

	fin_array = np.multiply(almost_fin_array,1)#water_array)

	fin_array[fin_array == 0] = np.nan

	# TODO: make some stats about image. 1) cloud fraction, 2) cloud free avg NDVI value 

	# Write NDVI raster 
	match_raster = [x for x in os.listdir(os.getcwd()) if "b01" in x and "clip" in x]
	match_raster = "".join(match_raster)

	mask.array2raster(match_raster,image.sat+"_"+image.year+ "_"+nearest+"_ndvi_8d.tif", fin_array)

	# Cleanup 
	files = [x for x in os.listdir(image.get_current_directory()) if not os.path.isdir(x) and "clip" in x]
	print(files)
	for f in files:
		os.remove(f)
		
	# Move to new directory 	
	file = [x for x in os.listdir(image.get_current_directory()) if x.endswith(".tif") and  "ndvi" in x]
	print(file)
	file = "".join(file)
	print(file)

	if os.path.exists(os.path.join(image.get_current_directory(),"historic")):
		pass
	else:
		os.mkdir(os.path.join(image.get_current_directory(),"historic"))
			
	new_dir = os.path.join(image.get_current_directory(),"historic")

	shutil.copy(file,new_dir)
	os.remove(file)
	
###########################################
# Small function to grab the current date #
###########################################


	
# dates = get_dates()

for year in aqua_years:
	ndvi_8day(product,str(year),str(doy),"Aqua")
	
for year in years:
	ndvi_8day(product,str(year),str(doy),"Terra")
