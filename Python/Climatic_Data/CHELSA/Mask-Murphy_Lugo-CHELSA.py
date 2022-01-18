#------------------------------------
#
# Mask each definition (ML)
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
climate=''
precip=''

#------------------------------------
# Copy and rename the files for Murphy and Lugo, FAO, DryFlor

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('ML.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(climate,'AP_chelsa_ML.tif')
        copyfile(path,dst)
        continue
    if fname.endswith('Dry.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(climate,'AP_chelsa_Dry.tif')
        copyfile(path,dst)
        continue
    if fname.endswith('FAO.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(climate,'AP_chelsa_FAO.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Out paths for defs
ex=gdal.Open('/monthly_ML.tif',GA_ReadOnly)
b=ex.GetRasterBand(1)
out_murph='/Murphy_Lugo.tif'

# Murphy and Lugo
crit_ml='*ML.tif'
q_ml=os.path.join(climate,crit_ml)
tile_ml=glob.glob(q_ml)
t_ml_arr=[]

for fp in tile_ml:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-32768]=0
    arr_2d=arr[0,:,:]
    t_ml_arr.append(arr_2d)
    fname.write(arr)

t_ml=np.array(t_ml_arr)
dataOut=t_ml.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_murph,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None

#Run function to mask 3 settings to 1
def final_murph(out_murph):
    with r.open(out_murph,'r+') as dset:
        arr=dset.read()
        arr[arr!=3]=0
        arr[arr==3]=1
        dset.write(arr)

final_murph(out_murph)
