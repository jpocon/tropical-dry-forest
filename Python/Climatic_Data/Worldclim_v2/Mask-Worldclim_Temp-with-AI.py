#--------------------------------------------------
#
# Masking compressed temp raster with each AI map
#
#--------------------------------------------------
# Set up the environment

print('Setting the environment...')
print()
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
print('Modules imported.')
print()

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory.')
print()
os.chdir('')
cwd=os.getcwd()

#--------------------------------------------------
# Copy and rename compressed temp raster for each AI threshold
print('Copying temp rasters to AI destinations...')
for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if fname.endswith('compress_clip.tif'):
        path=os.path.join(cwd,fname)
        dst1=os.path.join(cwd,'temp_65.tif')
        dst2=os.path.join(cwd,'temp_1.tif')
        copyfile(path,dst1)
        copyfile(path,dst2)
        continue
    if fname.endswith('65_tdf-clip.tif'):
        path=os.path.join(cwd,fname)
        dst=os.path.join(cwd,'AI_65.tif')
        copyfile(path,dst)
        continue
    if fname.endswith('0_tdf-clip.tif'):
        path=os.path.join(cwd,fname)
        dst=os.path.join(cwd,'AI_1.tif')
        copyfile(path,dst)
        continue
    else:
        continue
print('Files copied.')
print()

#--------------------------------------------------
# Paths for each raster
print('Setting the environment to merge temp and AI rasters...')
ex=gdal.Open('.tif',GA_ReadOnly) # choose compressed raster for its extent
b=ex.GetRasterBand(1)
AI_65='.tif'
AI_1='.tif'

# AI 0.65
crit_65='*65.tif'
q_65=os.path.join(cwd,crit_65)
tile_65=glob.glob(q_65)
arr65=[]

for fp in tile_65:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr_2d=arr[0,:,:]
    arr65.append(arr_2d)
    fname.write(arr)

ai65=np.array(arr65)
dataOut=ai65.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(AI_65,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

arr65=None

# AI 1.0
crit_1='*1.tif'
q_1=os.path.join(cwd,crit_1)
tile_1=glob.glob(q_1)
arr1=[]

for fp in tile_1:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr_2d=arr[0,:,:]
    arr1.append(arr_2d)
    fname.write(arr)

ai1=np.array(arr1)
dataOut=ai1.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(AI_1,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None
arr1=None

print('Rasters merged.')
print()

#--------------------------------------------------
# Mask compressed temp raster with both AI extents
print('Masking AI 0.65...')
def ai65(AI_65):
    with r.open(AI_65,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0        
        arr[arr==2]=1
        dset.write(arr)
        arr=None
        dset=None

ai65(AI_65)
print('AI 0.65 masked.')
print()

print('Masking AI 1.0...')
def ai1(AI_1):
    with r.open(AI_1,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[arr==2]=1
        dset.write(arr)
        arr=None
        dset=None

ai1(AI_1)
print('AI 1.0 masked.')
print()

print('The program is done running!')
