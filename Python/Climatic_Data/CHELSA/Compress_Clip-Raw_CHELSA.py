#--------------------------------------------------
#
# Compress GB rasters to MB and clip to 30N/s
# Includes:
#   converting nodata values
#   converting dtype
#   clipping extent
#
#--------------------------------------------------
# Set up the environment

print('Importing os, rasterio, numpy, and gdal...')
print()
import os
import rasterio
import numpy as np
from osgeo import gdal
print('Modules imported.')
print()

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory to desired folder.')
print()
os.chdir('') #directory pointing to raw CHELSA rasters
cwd=os.getcwd()

#--------------------------------------------------
# Function to compress and convert

def compress(fp):
    print('Calling in '+fp+' and converting nodataval to -32768 and dtype to int16...')
    print()
    with rasterio.open(fp,'r+',compress='LZW') as dset:
        dset.nodata=-32768
    src=rasterio.open(fp,'r',compress='LZW')
    meta=src.meta.copy()
    meta.update({'dtype':'int16'
                 }
                )
    arr=src.read()
    arr=arr.astype('int16')
    print('Saving '+fp+' ...')
    print()
    with rasterio.open(fp,'w',**meta,compress='LZW') as dest:
        dest.write(arr)
        dest=None
        arr=None
        src=None
        meta=None

#--------------------------------------------------
# Append all full path filenames to array to iterate

print('Compressing all tifs in current directory...')
print()

for r, d, f in os.walk(cwd):
    for file in f:
        if 'solrad' in file:
            fp=os.path.join(r,file).replace('\\', '/')
            compress(fp)

print('CHELSA rasters compressed.')
print()

#--------------------------------------------------
# Clip to 30N/S

def latlon(area):
    dset=gdal.Open(area,gdal.GA_ReadOnly)
    arr=dset.ReadAsArray()
    arr=arr.astype('int16')
    ds=gdal.Translate(area,dset,projWin=[-180.0001388888500173,
    30.0,179.9998596711500340,-30.0])
    ds=None

print('Clipping CHELSA rasters to 30N/S...')
for file in os.listdir(cwd):
     fname=os.fsdecode(file)
     if fname.endswith('.tif'):
         area=os.path.join(cwd,fname)
         print('Clipping '+area+' ...')
         latlon(area)
         print(area+' clipped.')
         print()
         continue
     else:
         continue
print('CHELSA rasters clipped.')
print()

print('The program in done running.')
