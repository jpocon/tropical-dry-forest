#------------------------------------
#
# Mask monthly 100mm thresholds from WorldClim v2
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
precip='' # folder with raw Worldclim precipitation files
dir100='' # destination folder

#------------------------------------
# Copy and rename the files for monthly 100mm

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(dir100,os.path.basename(os.path.splitext(path)[0])+'_mon100'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Define function to extract the months that are <= 100

def mon_100(mon100):
    with r.open(mon100,'r+') as dset:
        arr=dset.read()
        arr[arr==-32768]=101
        arr[arr<=100]=1
        arr[arr>100]=-32768
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate over each file

for file in os.listdir(dir100):
     fname=os.fsdecode(file)
     if fname.endswith('_mon100.tif'):
         mon100=os.path.join(dir100,fname)
         mon_100(mon100)
         continue
     else:
         continue

#------------------------------------
# Overlap monthly 100mm into one tiff

ex100=gdal.Open('/wc2_30s_prec_01_mon100.tif',GA_ReadOnly) # retrieve example tif for extent
b100=ex100.GetRasterBand(1)
out_100='.tif'

# 100mm file
crit_100='w*.tif'
q100=os.path.join(dir100,crit_100)
tile100=glob.glob(q100)
t100_arr=[]

for fp in tile100:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-32768]=0
    arr_2d=arr[0,:,:]
    t100_arr.append(arr_2d)
    fname.write(arr)

t100=np.array(t100_arr)
dataOut=t100.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_100,ex100.RasterXSize,ex100.RasterYSize,1,b100.DataType)
CopyDatasetInfo(ex100,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex100=None
bandOut=None
dsOut=None
