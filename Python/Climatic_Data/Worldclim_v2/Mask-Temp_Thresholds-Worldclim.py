#------------------------------------
#
# Mask temps for TDF definitions
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
temp=''
mask=''

#------------------------------------
# Copy and rename the files for no freeze

for file in os.listdir(temp):
    fname=os.fsdecode(file)
    if fname.endswith('06.tif'):
        path=os.path.join(temp,fname)
        dst=os.path.join(mask,os.path.basename(os.path.splitext(path)[0])+'_n0'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Copy and rename the files for >17C

for file in os.listdir(temp):
    fname=os.fsdecode(file)
    if fname.endswith('01.tif'):
        path=os.path.join(temp,fname)
        dst=os.path.join(mask,os.path.basename(os.path.splitext(path)[0])+'_17'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Define function to extract temps that are >0C

def no_frz(n0):
    with r.open(n0,'r+') as dset:
        arr=dset.read()
        arr[arr>0]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate no freeze file

for file in os.listdir(mask):
     fname=os.fsdecode(file)
     if fname.endswith('_n0.tif'):
         n0=os.path.join(mask,fname)
         no_frz(n0)
         continue
     else:
         continue

#------------------------------------
# Define function to extract temps that are >17C

def MAT(mat):
    with r.open(mat,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[arr>=17]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate MAT file

for file in os.listdir(mask):
     fname=os.fsdecode(file)
     if fname.endswith('_17.tif'):
         mat=os.path.join(mask,fname)
         MAT(mat)
         continue
     else:
         continue

#------------------------------------
# Overlap/mask no freeze and MAT

# Read original tiff to get raster sizes and bands for writing new files
ex=gdal.Open('.tif',GA_ReadOnly)
b=ex.GetRasterBand(1)

# Set filepath for saving
out_temp='.tif'

# call in new files
crit_temp='w*.tif'
q=os.path.join(mask,crit_temp)
t=glob.glob(q)
t_arr=[]

for fp in t:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-1.7e+308]=0
    arr_2d=arr[0,:,:]
    t_arr.append(arr_2d)
    fname.write(arr)

tarr=np.array(t_arr)
dataOut=tarr.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_temp,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None

#------------------------------------
# Mask the overlaps from value 2 to 1 for overlaying into final climate maps
temps='.tif'
with r.open(temps,'r+') as dset:
    arr=dset.read()
    arr[arr==1]=0
    arr[arr==2]=1
    dset.write(arr)
    
