#--------------------------------------------------
#
# Compress rasters
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
import numpy
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

def binary(fp):
    with rasterio.open(fp,'r+') as data:
        arr=data.read()
        arr[arr>0]=1
        arr[arr<=0]=0
        data.write(arr)

#--------------------------------------------------
# Append all full path filenames to array to iterate

print('Compressing treecover tiles...')
print()

folders=["fao_10_00", "fao_10_18", "fao_40_00", "fao_40_18",
         "ml_10_00", "ml_10_18", "ml_40_00", "ml_40_18"] # Google Earth Engine export folders

for i in range(0,8):
    for file in os.listdir(folders[i]):
        if file.endswith('.tif'):
            fp=os.path.join(folders[i],file)
            print('Calling in '+fp+' and converting to binary...')
            binary(fp)
            print(fp+' converted to binary...')
            compress(fp)
            print(fp+' compressed.')
            print()
            continue
        else:
            continue

print('The program in done running.')
