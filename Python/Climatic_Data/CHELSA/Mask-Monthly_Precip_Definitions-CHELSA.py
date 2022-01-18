#------------------------------------
#
# Mask for each climatic definition
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
precip='/monthly_100'
temp=''

#------------------------------------
# Copy and rename the files for Murphy and Lugo, FAO, DryFlor

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('100mm-chelsa.tif'):
        path=os.path.join(precip,fname)
        dst1=os.path.join(climate,'monthly_ML.tif')
        dst2=os.path.join(climate,'monthly_FAO.tif')
        dst3=os.path.join(climate,'monthly_Dry.tif')
        copyfile(path,dst1)
        copyfile(path,dst2)
        copyfile(path,dst3)
        continue
    else:
        continue

for file in os.listdir(temp):
    fname=os.fsdecode(file)
    if fname.endswith('temp-chelsa.tif'):
        path=os.path.join(temp,fname)
        dst1=os.path.join(climate,'temp_ML.tif')
        dst2=os.path.join(climate,'temp_FAO.tif')
        dst3=os.path.join(climate,'temp_Dry.tif')
        copyfile(path,dst1)
        copyfile(path,dst2)
        copyfile(path,dst3)
        continue
    else:
        continue

#------------------------------------
# Climatic maps file paths
murph='/monthly_ML.tif'
dry='/monthly_Dry.tif'
fao='/monthly_FAO.tif'

#------------------------------------
# Define function to extract the dry seasons for each definition

def ml(murph):
    with r.open(murph,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[(arr>=4)&(arr<=7)]=1
        arr[arr!=1]=0
        dset.write(arr)

def df(dry):
    with r.open(dry,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[(arr>=3)&(arr<=6)]=1
        arr[arr!=1]=0
        dset.write(arr)

def f(fao):
    with r.open(fao,'r+') as dset:
        arr=dset.read()
        arr[arr==1]=0
        arr[(arr>=5)&(arr<=8)]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run each function to mask the dry seasons in each deifnition
ml(murph)
df(dry)
f(fao)
