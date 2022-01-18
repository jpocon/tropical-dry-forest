#--------------------------------------------------
#
# Compress GB rasters to MB
# Includes:
#   converting nodata values
#   converting dtype
#
#--------------------------------------------------
# Set up the environment

print('Importing os, rasterio, and numpy...')
print()
import os
import rasterio
import numpy as np
print('Modules imported.')
print()

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory to desired folder.')
print()
os.chdir('')
cwd=os.getcwd()

#--------------------------------------------------
# Function to compress and convert

def compress(fp):
    print('Calling in '+fp+' and converting nodataval to 0 and dtype to uint8...')
    print()
    with rasterio.open(fp,'r+',compress='LZW') as dset:
        dset.nodata=0
    src=rasterio.open(fp,'r',compress='LZW')
    meta=src.meta.copy()
    meta.update({'dtype':'uint8'
                 }
                )
    arr=src.read()
    arr=arr.astype('uint8')
    print('Saving '+fp+' ...')
    print()
    with rasterio.open(fp,'w',**meta,compress='LZW') as dest:
        dest.write(arr)

#--------------------------------------------------
# Append all full path filenames to array to iterate

print('Compressing CHELSA AI in current directory...')
print()

# Call in the large raster
fp='.tif'
with rasterio.open(fp,'r+') as data:
    arr=data.read()
    arr[arr>0]=1
    arr[arr<=0]=0
    data.write(arr)
compress(fp)

print('The program in done running.')
