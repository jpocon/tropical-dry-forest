#--------------------------------------------------
#
# CHELSA - Total Annual Precipitation
# (1) Sum annual precip
# (2) Compress resulting summed raster (uint16)
#
# Includes:
#   summing arrays
#   converting dtype/nodata
#
# J.P. Ocón - March 16, 2022
#
#--------------------------------------------------
# Set up the environment

print('Importing os, glob, rasterio, numpy, and gdal...')
print()
import os
import glob
import rasterio as r
import numpy as np
from osgeo import gdal, gdalnumeric
print('Modules imported.')
print()

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory to desired folder.')
cwd=os.getcwd()
print('Directory set to: ' + cwd)
print()

#--------------------------------------------------
# Overlap monthly 100mm into one tiff

print("Querying monthly rasters...")
crit100 = "*_V.2.1.tif" # criteria for glob files
q100 = os.path.join(cwd, crit100) # query filepaths
tile100 = glob.glob(q100) # glob tiles
t100 = [] # create empty array to append monthly data
print("Ready to append monthly data into one raster..." + "\n")

for fp in tile100:
    print("Appending " + fp + " ...")
    dset = r.open(fp, "r+")
    arr = dset.read()
    arr2d = arr[0, :, :] # take only the first band of the .tif
    del arr, dset, fp
    t100.append(arr2d)
    del arr2d
    print("Array appended." + "\n")

print("Summing all monthly rasters...")
t100_arr = np.array(t100)
del t100
dataOut = t100_arr.sum(axis = 0) # sum the monthly vals
dataOut = dataOut*0.1 # adjusting for CHELSA scaling
del t100_arr
print("Monthly rasters summed." + "\n")

print("Calling January .tif for metadata...")
# read in one of the monthly .tif's to reference spatial metadata later
ex = gdal.Open(cwd + "/CHELSA_pr_01_1981-2010_V.2.1.tif", gdal.GA_ReadOnly)
b = ex.GetRasterBand(1)
# filename for output file
out = cwd + "/CHELSA_pr_01-12_1981-2010_V.2.1.tif"
print("Spatial parameters set." + "\n")

print("Saving summed raster...")
driver = gdal.GetDriverByName("GTiff") # saving as GeoTiff
dsOut = driver.Create(out,
                      ex.RasterXSize,
                      ex.RasterYSize,
                      1,
                      b.DataType)
gdalnumeric.CopyDatasetInfo(ex, dsOut) # copy over all metadata from ref .tif to new file
bandOut = dsOut.GetRasterBand(1)
gdalnumeric.BandWriteArray(bandOut, dataOut) # write the appended array to new file band

del dsOut, bandOut, dataOut, ex
print("Monthly CHELSA summed.")

#--------------------------------------------------
# Function to compress and convert

def compress(fp):
    print('Calling in '+fp+' and converting nodataval to 0 and dtype to uint16...')
    print()
    with r.open(fp,'r+',compress='LZW') as dset:
        dset.nodata = 0
    src=r.open(fp,'r',compress='LZW')
    meta=src.meta.copy()
    meta.update({'dtype':'uint16'
                 }
                )
    arr=src.read()
    arr=arr.astype('uint16')
    print('Saving '+fp+' ...')
    print()
    with r.open(fp,'w',**meta,compress='LZW') as dest:
        dest.write(arr)

#--------------------------------------------------
# Append all full path filenames to array to iterate

print('Compressing all tifs in current directory...')
print()

for file in os.listdir(cwd):
    if 'pr_01-12' in file:
        fp=os.path.join(cwd,file)
        compress(fp)

print('The program in done running.')
