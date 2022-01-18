#------------------------------------
#
# Mask AP for TDF definitions
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
precip='' # folder with raw Worldclim precipitation files
mask='' # destination folder

#------------------------------------
# Copy and rename the files for each AP

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('bio_30s_12.tif'):
        path=os.path.join(precip,fname)
        dst1=os.path.join(mask,os.path.basename(os.path.splitext(path)[0])+'_ML'+'.tif')
        dst2=os.path.join(mask,os.path.basename(os.path.splitext(path)[0])+'_Dry'+'.tif')
        dst3=os.path.join(mask,os.path.basename(os.path.splitext(path)[0])+'_FAO'+'.tif')
        copyfile(path,dst1)
        copyfile(path,dst2)
        copyfile(path,dst3)
        continue
    else:
        continue

#------------------------------------
# Define function to extract AP 250-2000mm (Murphy and Lugo)

def ML(ml):
    with r.open(ml,'r+') as dset:
        arr=dset.read()
        arr[(arr>=250)&(arr<=2000)]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate for Murphy and Lugo AP

for file in os.listdir(mask):
     fname=os.fsdecode(file)
     if fname.endswith('_ML.tif'):
         ml=os.path.join(mask,fname)
         ML(ml)
         continue
     else:
         continue

#------------------------------------
# Define function to extract AP 500-1500mm (FAO)

def FAO(fao):
    with r.open(fao,'r+') as dset:
        arr=dset.read()
        arr[(arr>=500)&(arr<=1500)]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate for FAO AP

for file in os.listdir(mask):
     fname=os.fsdecode(file)
     if fname.endswith('_FAO.tif'):
         fao=os.path.join(mask,fname)
         FAO(fao)
         continue
     else:
         continue

#------------------------------------
# Define function to extract AP <1800mm (DryFlor)

def Dry(dry):
    with r.open(dry,'r+') as dset:
        arr=dset.read()
        arr[(arr<=1800)&(arr>=0)]=1
        arr[arr!=1]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate for DryFlor AP

for file in os.listdir(mask):
     fname=os.fsdecode(file)
     if fname.endswith('_Dry.tif'):
         dry=os.path.join(mask,fname)
         Dry(dry)
         continue
     else:
         continue
