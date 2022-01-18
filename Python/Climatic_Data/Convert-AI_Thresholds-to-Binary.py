#--------------------------------------------------
#
# Binary AI Maps:
#   band math
#   convert nodata + dtype
#
#--------------------------------------------------
# Set up the environment

print('Importing os, rasterio, and numpy...')
print()
import os
import rasterio as r
import numpy as np
import glob
print('Modules imported.')
print()

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory to desired folder.')
print()
os.chdir('')
output = ''
cwd = os.getcwd()

#--------------------------------------------------
# Function to compress .tif and convert nodata / dtype

def compress(fp):
    print('Calling in ' + fp + ' and converting nodataval to 0 and dtype to uint8...')
    print()
    with r.open(fp, 'r+', compress = 'LZW') as dset:
        dset.nodata = 0
    src = r.open(fp, 'r', compress = 'LZW')
    meta = src.meta.copy()
    meta.update({'dtype':'uint8'
                 }
                )
    arr = src.read()
    arr[arr <= 0.2] = 0
    arr[arr >= 0.65] = 0
    arr[arr >= 0.2] = 1
    arr = arr.astype('uint8')
    print('Saving ' + fp + ' ...')
    print()
    with r.open(fp, 'w', **meta, compress = 'LZW') as dest:
        dest.write(arr)

#--------------------------------------------------
# 0.2 < AI < 0.65
crit_ai = '*humid.tif'
q_ai = os.path.join(cwd, crit_ai)
tile_ai = glob.glob(q_ai)

for fp in tile_ai:
    compress(fp)

print('The program in done running.')
