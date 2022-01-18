#------------------------------------
#
# Mask each definition
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
climate=''

# Out paths for defs
ex=gdal.Open('/monthly_ML.tif',GA_ReadOnly)
b=ex.GetRasterBand(1)
out_dry='/DryFlor.tif'

# DryFlor
crit_dry='*Dry.tif'
q_dry=os.path.join(climate,crit_dry)
tile_dry=glob.glob(q_dry)
t_dry_arr=[]

for fp in tile_dry:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-32768]=0
    arr_2d=arr[0,:,:]
    t_dry_arr.append(arr_2d)
    fname.write(arr)

t_dry=np.array(t_dry_arr)
dataOut=t_dry.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_dry,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None

#Run function to mask 3 settings to 1
def final_dry(out_dry):
    with r.open(out_dry,'r+') as dset:
        arr=dset.read()
        arr[arr!=3]=0
        arr[arr==3]=1
        dset.write(arr)

final_dry(out_dry)
