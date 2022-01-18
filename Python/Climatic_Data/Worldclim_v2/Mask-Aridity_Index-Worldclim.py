#------------------------------------
#
# Mask temps for Aridity (1.0, and 0.65)
#
#------------------------------------
# Set up the environment

# Files
import os
import fnmatch as fn
from shutil import copyfile
import glob

# Maths
import numpy as np
from numpy import *

# GIS
import rasterio as r
from osgeo import gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

#------------------------------------
# Folder and file paths
ai='/ai_et0' # raw airidty file
#------------------------------------
# Copy and rename the files for tdf 1.0

for file in os.listdir(ai):
    fname=os.fsdecode(file)
    if fname.endswith('et0.tif'):
        path=os.path.join(ai,fname)
        dst=os.path.join(ai,os.path.basename(os.path.splitext(path)[0])+'_1.0_tdf'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Copy and rename the files for tdf 0.65

for file in os.listdir(ai):
    fname=os.fsdecode(file)
    if fname.endswith('et0.tif'):
        path=os.path.join(ai,fname)
        dst=os.path.join(ai,os.path.basename(os.path.splitext(path)[0])+'_0.65_tdf'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Define function to 1.0 AI

def tdf_01(tdf01):
    with r.open(tdf01,'r+') as dset:
        arr=dset.read()
        arr[arr<0]=11000
        arr[arr<=10000]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate tdf 1.0

for file in os.listdir(ai):
     fname=os.fsdecode(file)
     if fname.endswith('_1.0_tdf.tif'):
         tdf01=os.path.join(ai,fname)
         tdf_01(tdf01)
         continue
     else:
         continue

#------------------------------------
# Define function to 0.65 AI

def tdf_65(tdf65):
    with r.open(tdf65,'r+') as dset:
        arr=dset.read()
        arr[arr<0]=11000
        arr[arr<=6500]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate MAT file

for file in os.listdir(ai):
     fname=os.fsdecode(file)
     if fname.endswith('_0.65_tdf.tif'):
         tdf65=os.path.join(ai,fname)
         tdf_65(tdf65)
         continue
     else:
         continue

#------------------------------------
# Overlap/mask no freeze and MAT

# Read original tiff to get raster sizes and bands for writing new files
ex=gdal.Open('/ai_et0.tif',GA_ReadOnly)
b=ex.GetRasterBand(1)

# Set filepath for saving
out_ai='.tif'

# call in new files
crit_ai='*tdf.tif'
q=os.path.join(ai,crit_ai)
a=glob.glob(q)
a_arr=[]

for fp in a:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-2.14748e+09]=0
    arr_2d=arr[0,:,:]
    a_arr.append(arr_2d)
    fname.write(arr)

aarr=np.array(a_arr)
dataOut=aarr.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_ai,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None
