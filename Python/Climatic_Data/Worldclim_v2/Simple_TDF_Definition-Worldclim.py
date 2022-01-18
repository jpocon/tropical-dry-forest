#--------------------------------------------------
#
# Worldclim v2
# (1) Compress monthly raster files & extract monthly precip thresholds
# (2) Mask annual precip thresholds for each definition
# (3) Mask precip seasonality for each definition
# (4) Mask temp thresholds for each definition
# (5) Overlay per definition
#
# Includes:
#   converting nodata values
#   converting dtype
#
#--------------------------------------------------
# Set up the environment

print("Setting up the environment...")
# Base modules
import os
import gc
import fnmatch as fn
from shutil import copyfile as cp
import glob
import numpy as np
from numpy import *
import re

# Spatial modules
import rasterio as r
from osgeo import gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
print("Modules imported." + "\n")

#--------------------------------------------------
# Folder and file paths

print("Setting directories...")
os.chdir("")
prec = ""
bio = ""
mon = ""
out = ""
print("Directories set." + "\n")

#--------------------------------------------------
# Copy and rename the files for monthly 100mm

print("Copying raw data for processing...")
for file in os.listdir(prec):
    fname = os.fsdecode(file)
    if fname.endswith(".tif"):
        path = os.path.join(prec, fname)
        dst = os.path.join(mon,
                           os.path.basename(os.path.splitext(path)[0]) +
                           "_mon100" +
                           ".tif")
        cp(path, dst)
        
print("Data ready for processing." + "\n")

#--------------------------------------------------
# Functions to extract the months that are <= 100 and to compress the resulting .tif's

def clip(mon100, mon_clip):
    print("Clipping " + mon100 + " to 30N and 30S...")
    dset = gdal.Open(mon100, gdal.GA_ReadOnly)
    ds = gdal.Translate(mon_clip, dset,
                        projWin = [-180.0, 30.0, 180.0, -30.0],
                        outputType = gdalconst.GDT_Byte)
    print(mon_clip + " clipped." + "\n")

for file in os.listdir(mon):
     fname = os.fsdecode(file)
     if fname.endswith("_mon100.tif"):
        mon100 = os.path.join(mon, fname)
        mon_clip = os.path.join(mon,
                           os.path.basename(os.path.splitext(mon100)[0]) +
                           "_tropics" +
                           ".tif")
        print("Calling in " + mon100 + " ...")
        clip(mon100, mon_clip)

print("Clipped monthly rasters." + "\n")

def mon_100(mon100):
    print("Assigning pixels with <= 100 mm precip to vals of 1, and the rest to 0...")
    with r.open(mon100, "r+") as dset:
        arr = dset.read()
        arr[arr == -32768] = 101 # nodata val
        arr[(arr > 0) & (arr <= 100)] = 1
        arr[arr > 100] = -32768 # set all vals > 100 to nodata val
        dset.write(arr)

def compress(mon100):
    print("Compressing " +
          mon100 +
          " and converting nodata val to 0 and dtype to uint8...")
    with r.open(mon100, "r+", compress = "LZW") as dset:
        dset.nodata = 0
    src = r.open(mon100, "r", compress = "LZW")
    meta = src.meta.copy()
    meta.update({"dtype":"uint8"})
    arr = src.read()
    arr = arr.astype("uint8")
    print("Saving " +
          os.path.splitext(mon100)[0] +
          "_compressed" + ".tif" + " to current working directory...")
    with r.open(os.path.splitext(mon100)[0] +
                "_compressed" + ".tif", "w", **meta,compress = "LZW") as dest:
        dest.write(arr)
        print(mon100 + " successfully subsetted and compressed" + "\n")

for file in os.listdir(mon):
     fname = os.fsdecode(file)
     if fname.endswith("_tropics.tif"):
        mon100 = os.path.join(mon, fname)
        print("Calling in " + mon100 + " ...")
        mon_100(mon100)
        compress(mon100)

print("Monthly .tif's subsetted and compressed." + "\n")

#--------------------------------------------------
# Overlap monthly 100mm into one tiff

print("Querying monthly rasters...")
crit100 = "*_compressed.tif" # criteria for glob files
q100 = os.path.join(mon, crit100) # query filepaths
tile100 = glob.glob(q100) # glob tiles
t100 = [] # create empty array to append monthly data
print("Ready to append monthly data into one raster..." + "\n")

for fp in tile100:
    print("Appending " + fp + " ...")
    dset = r.open(fp, "r+")
    arr = dset.read()
    arr[arr == -32768] = 0 # change nodata vals to 0
    arr = arr.astype("uint8")
    arr2d = arr[0, :, :] # take only the first band of the .tif
    del arr, dset, fp
    t100.append(arr2d)
    gc.collect()
    del arr2d
    print("Array appended." + "\n")

print("Summing all monthly rasters...")
t100_arr = np.array(t100)
del t100
gc.collect()
dataOut = t100_arr.sum(axis = 0) # sum the monthly vals
gc.collect()
del t100_arr
print("Monthly rasters summed." + "\n")

print("Calling January .tif for metadata...")
# read in one of the monthly .tif's to reference spatial metadata later
ex100 = gdal.Open(mon + "/wc2_30s_prec_01_mon100_tropics_compressed.tif", gdal.GA_ReadOnly)
# we are only interested in the first band of the .tif (the precip data)
b100 = ex100.GetRasterBand(1)
# filename for output file
out100 = prec + "/wc2-monthly_dry.tif"
print("Spatial parameters set." + "\n")

print("Saving summed raster...")
driver = gdal.GetDriverByName("GTiff") # saving as GeoTiff
dsOut = driver.Create(out100, # filename
                      ex100.RasterXSize, # array dims from ref .tif
                      ex100.RasterYSize, 
                      1, # num of bands
                      b100.DataType) # ref .tif datatype
CopyDatasetInfo(ex100, dsOut) # copy over all metadata from ref .tif to new file
bandOut = dsOut.GetRasterBand(1)
BandWriteArray(bandOut, dataOut) # write the appended array to new file band

del dsOut, bandOut, dataOut, ex100, b100
print("Time to build the simple def!")

#--------------------------------------------------
# Define function to extract =< 2000mm
def Simple(simple):
    with r.open(simple,'r+') as dset:
        arr=dset.read()
        arr[arr<=2000]=1
        arr[arr!=1]=0
        dset.write(arr)

# copy raw data file into Simple def
for file in os.listdir(bio):
    fname=os.fsdecode(file)
    if fname.endswith('bio_30s_12.tif'):
        path = os.path.join(bio, fname)
        print("Subsetting annual precip =< 2000 mm...")
        simple = os.path.join(out, "annual_precip.tif")
        cp(path, simple)
        Simple(simple)
        print("AP ready!")
        continue
    else:
        continue

#--------------------------------------------------
# Define function to extract the dry seasons for each definition

def seas_simple(simple):
    with r.open(simple,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[arr>=4]=1
        arr[arr!=1]=0
        dset.write(arr)

# copy raw data file into Simple def
for file in os.listdir(prec):
    fname=os.fsdecode(file)
    if fname.endswith('wc2-monthly_dry.tif'):
        path = os.path.join(prec, fname)
        print("Subsetting seasonality >= 4 mos...")
        simple = os.path.join(out, "seasonal.tif")
        cp(path, simple)
        seas_simple(simple)
        print("Seasonality ready!")
        continue
    else:
        continue

#--------------------------------------------------
# Define function to extract temps that are >0C
def no_frz(n0):
    with r.open(n0,'r+') as dset:
        arr=dset.read()
        arr[arr>0]=1
        arr[arr!=1]=0
        dset.write(arr)

# Copy and rename the files for temp masks

for file in os.listdir(bio):
    fname = os.fsdecode(file)
    if fname.endswith('bio_30s_06.tif'):
        path = os.path.join(bio, fname)
        print("Subsetting no freeze...")
        n0 = os.path.join(out, "no_freeze.tif")
        cp(path, n0)
        no_frz(n0)
        print("No freeze ready!")
        continue
    else:
        continue

def clip2(bios, clip_simple):
    print("Clipping " + bios + " to 30N and 30S...")
    dset = gdal.Open(bios, gdal.GA_ReadOnly)
    ds = gdal.Translate(clip_simple, dset,
                        projWin = [-180.0, 30.0, 180.0, -30.0],
                        outputType = gdalconst.GDT_Byte)
    print(clip_simple + " clipped." + "\n")

for file in os.listdir(out):
     fname = os.fsdecode(file)
     if fname.endswith(".tif"):
        bios = os.path.join(out, fname)
        clip_simple = os.path.join(out,
                                   os.path.basename(os.path.splitext(bios)[0]) +
                                   "_simple" +
                                   ".tif")
        print("Calling in " + bios + " ...")
        clip2(bios, clip_simple)

print("Clipped simple rasters." + "\n")

#--------------------------------------------------
# function for thresholds of all definitions

def final_def(x):
    with r.open(x, 'r+') as dset:
        arr = dset.read() # read in the array
        arr[arr!=3] = 0 # if the pixel val does not have all three climatic parameters, set to 0
        arr[arr==3] = 1 # if all 3 params are met, set to 1
        dset.write(arr) # save the .tif

#--------------------------------------------------
# Out paths for defs
out_simple = 'Simple.tif'

crit_simple = '*_simple.tif'
q_simple = os.path.join(out, crit_simple)
tile_simple = glob.glob(q_simple)
t_simple_arr = []

print("Appending simple rasters to one array...")
for fp in tile_simple: # for each file name in the re generated list
    fname = r.open(fp, 'r+') # open the .tif in read format
    arr = fname.read()
    arr[arr==-1.7e+308] = 0 # read in the array and set no data vals to 0
    arr_2d = arr[0,:,:] # flatten the array by only calling one band
    t_simple_arr.append(arr_2d) # append to empty array
    fname.write(arr) # save .tif

print("Summing simple rasters...")
t_simple = np.array(t_simple_arr)
dataOut = t_simple.sum(axis = 0) # sum the stack of arrays into one .tif

print("Rasters summed!")
print("Saving summed raster to new file...")
ex = gdal.Open(out + "/seasonal_simple.tif", gdal.GA_ReadOnly) # get georeference metadata from existing file
b = ex.GetRasterBand(1) # retrive all band metadata from this file
driver = gdal.GetDriverByName('GTiff') # .tif file format
dsOut = driver.Create(out_simple, # using output name
                      ex.RasterXSize, # set the array dimensions to match existing file
                      ex.RasterYSize,
                      1, # number of bands in .tif
                      b.DataType) # match data types
CopyDatasetInfo(ex, dsOut) 
bandOut = dsOut.GetRasterBand(1)
BandWriteArray(bandOut, dataOut)

ex=None
bandOut=None
dsOut=None

print("Simple raster saved.")
print("Subsetting simple defintion...")
final_def(out_simple) # run the newly stacked/summed file and extract where pixels = 3
print("The program is done running!")
