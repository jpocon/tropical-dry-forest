#------------------------------------
#
# Download Hansen Datamask files and save VRT
#
#------------------------------------
# Set up the environment

import urllib.request
import os
import rasterio
import numpy as np
from rasterio.merge import merge
from rasterio.plot import show
from osgeo import gdal

#------------------------------------
# Define the variables

os.chdir('/datamask')
vrt_dmask='datamask.vrt'

#------------------------------------
# Download the files and iterate over each
# to mask land values, and append to VRT and mosaic arrays

print('Beginning file download with urllib2...')
print()

src_dmask=[]

# Call in each URL in txt file and download from server
with open('datamask.txt','r') as download:
    for line in download:
        
        print('Downloading...',str(line))
        
        tiff=os.path.basename(os.path.splitext(line)[0])+'.tif'
        urllib.request.urlretrieve(line,tiff)

        print('Appending tif to VRT array...')
        print()
        
        src_dmask.append(tiff)
        continue

print('All files downloaded.')
print()

print('Creating the VRT...')
print()

# Create VRT and save.
my_vrt=gdal.BuildVRT(vrt_dmask,src_dmask)
my_vrt=None

print('VRT saved.')
print()

print('The program is finished running.')
