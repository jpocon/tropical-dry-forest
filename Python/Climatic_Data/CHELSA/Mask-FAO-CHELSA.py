#------------------------------------
#
# Mask each definition (FAO)
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
out_fao='/FAO.tif'

# FAO
crit_fao='*FAO.tif'
q_fao=os.path.join(climate,crit_fao)
tile_fao=glob.glob(q_fao)
t_fao_arr=[]

for fp in tile_fao:
    fname=r.open(fp,'r+')
    arr=fname.read()
    arr[arr==-32768]=0
    arr_2d=arr[0,:,:]
    t_fao_arr.append(arr_2d)
    fname.write(arr)

t_fao=np.array(t_fao_arr)
dataOut=t_fao.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_fao,ex.RasterXSize,ex.RasterYSize,1,b.DataType)
CopyDatasetInfo(ex,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex=None
bandOut=None
dsOut=None

#Run function to mask 3 settings to 1
def final_fao(out_fao):
    with r.open(out_fao,'r+') as dset:
        arr=dset.read()
        arr[arr!=3]=0
        arr[arr==3]=1
        dset.write(arr)

final_fao(out_fao)
