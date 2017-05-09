This readme describes the software routines for a regional flood detection and impact assessment module for Southeast Asia.

Developer: Aakash Ahamed (aakash.ahamed@nasa.gov; ahamednasa@gmail.com;  908 370 7738)
Contact and PI: John D. Bolten (john.bolten@nasa.gov)
Date: 1/2017

NASA Goddard Space Flight Center
Hydrological Sciences Laboratory (Code 617)
Applied Sciences Program 

DISCLAIMER: THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This system has been tested, ported, and run successfully on Ubuntu 14.04 OS with 16 GB ram. Initial storage requirements are ~ 30 GB, but accumulated results will soon exceed this. As such, either more storage or maintenance scripts are recommended. 

################################
##### Install Instructions #####
################################ 

1) Start with a fresh install of Ubuntu 14.04 (16 GB ram, 50+ GB storage, 1 core required)
2) Install the gdal library and its dependencies by running (from command line)

	sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update

	sudo apt-get install gdal-bin python3-gdal

3) Test the gdal install by entering the python shell and typing "from osgeo import gdal"

	python3

	> from osgeo import gdal

4) If there are no error messages, the software libraries and dependencies are correctly configured.

5) Install the crontab commands in the crontab.txt file. Relative file paths and execution times may have to be altered.

##############################
##### Usage Instructions #####
##############################

**************************************
******** Folder : NRT_Flood **********		Routines for near real-time and historic flooded areas products
**************************************

Description: This set of programs downloads, processes, and applies classification algorithm to MODIS data to determine flooded areas, non-flooded areas, and areas where there is not enough data (due to cloud obscuration) to make a determination.

Programs in this module and descriptions:
1) fetch_historic.py - determines the day of the year and downloads M*D09Q1 products for the length of the MODIS Aqua and Terra record for the nearest available date on each year. If products in the "historic" directory are current, the program does not run. For each image, clouds are masked and per-pixel NDVI is calculated (250m resolution).
2) calc_historic_average.py - averages the Aqua / Terra files in the "historic" dir to construct a combined MODIS Aqua / Terra average NDVI composite for that day of the year. Outputs are stored in the "reference_files" directory
3) fetch_NRT.py - fetches the latest MODIS Aqua / Terra data available on the NASA nrt servers. Clouds are filtered and NDVI is calculated at 250m resolution. Outputs are stored in the "daily_files_nrt" directory.
4) calc_NRT_flood.py - calculates the NDVI of permanent water bodies using the "permanent_water_mask.tif" file and the historic NDVI composite. Pixels below the mean H2O NDVI threshold are classified "surface water". Outputs are surface water extent products (0 = No water, 1 = water, NaN = clouds/permanent water) generated for the day of year (avg_SW) and the latest NRT 4-day composite (SW_doy_year).

The programs should be run in the following order: 

(1) --> (2) --> (3) --> (4)

Required static file: 
1) permanent_water_mask.tif - This file is a binary .tiff raster of the MOD44 Permanent water mask 

Each program contains a more detailed description within the docstring of the .py file as well as in-line comments for code clarity