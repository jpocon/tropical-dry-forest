#------------------------------------
#
# (1) Reading and writing (GPM-IMERG) HDF5 Files
# (2) Extracting Monthly Avgs from IMERG
# (3) Copy and save as 60mm/month, 100mm/month
# (4) Overlap all 51 months / 4.25 for precip final
# (5) Mask the 60/100mm thresholds from WorldClim v2
#
#------------------------------------
# Set up the environment

import os
import fnmatch as fn
from shutil import copyfile
import numpy as np
from osgeo import gdal, osr
import rasterio as r
from affine import Affine
from rasterio.merge import merge
from shutil import copyfile
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import glob
from numpy import *

#------------------------------------
# Define the variables

pull='/IMERG'
save='/IMERG/gpm_daily'
month='/IMERG/gpm_monthly'
# List of coords doesn't come in a,b,c,d,e,f order - source world file (tfw) has the correct order
b,c,e,f,a,d=np.loadtxt('/IMERG/3B-MO-GIS.MS.MRG.3IMERG.20140401-S000000-E235959.04.V05B.tfw')
geotrans=[a,b,c,d,e,f]
bnd=1 # Index of the precip band

# (3)
seas='/IMERG/seasonality'

# (4)/(5) folders/files are dependent on functions from #1-3 and are defined below

#------------------------------------
# (1) Function to extract single band from HDF5 and save as GTiff

def hdf_band_extract(hdf,save,bnd):
    # Open HDF5
    dset=gdal.Open(hdf,gdal.GA_ReadOnly)
    band=gdal.Open(dset.GetSubDatasets()[bnd][0],gdal.GA_ReadOnly)
    # Read the band into an array
    arr=band.ReadAsArray()
    arr=arr[:,::-1]
    arr=arr.T
    # Convert nan values
    arr[arr==-9999.9]=np.nan
    # Convert to mm/day
    arr=arr*24.
    # Build save path
    fpath=os.path.join(save,os.path.basename(os.path.splitext(hdf)[0])+'_B'+str(bnd+1)+'_mmday'+'.tiff')
    # Write GTiff
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        band.RasterYSize,
        band.RasterXSize,
        1,# Number of bands
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform(geotrans)
    tiff.GetRasterBand(1).WriteArray(arr)
    tiff.GetRasterBand(1).SetNoDataValue(-9999.9)
    tiff=None # Close file and write to 'save'

# Run a for loop to iterate over each file
for file in os.listdir(pull):
     fname=os.fsdecode(file)
     if fname.endswith('.HDF5'):
         hdf=os.path.join(pull,fname)
         hdf_band_extract(hdf,save,bnd)
         continue
     else:
         continue

#------------------------------------
#
# (2) Convert from mm/day to mm/month
#
#------------------------------------
# Function for months with 31 days

def monthly31_avg(tif31,month):
    ds31=gdal.Open(tif31,gdal.GA_ReadOnly)
    arr=ds31.ReadAsArray()
    arr=arr*31.
    fpath=os.path.join(month,os.path.basename(os.path.splitext(tif31)[0])+'_mmM'+'.tiff')
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        3600, # Y-size
        1800, # X-size
        1,
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform(geotrans)
    tiff.GetRasterBand(1).WriteArray(arr)
    tiff.GetRasterBand(1).SetNoDataValue(-9999.9)
    tiff=None

# Function for months with 30 days
def monthly30_avg(tif30,month):
    ds30=gdal.Open(tif30,gdal.GA_ReadOnly)
    arr=ds30.ReadAsArray()
    arr=arr*30.
    fpath=os.path.join(month,os.path.basename(os.path.splitext(tif30)[0])+'_mmM'+'.tiff')
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        3600,
        1800,
        1,
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform(geotrans)
    tiff.GetRasterBand(1).WriteArray(arr)
    tiff.GetRasterBand(1).SetNoDataValue(-9999.9)
    tiff=None

# Function for non-leap year Feb    
def monthly28_avg(tif28,month):
    ds28=gdal.Open(tif28,gdal.GA_ReadOnly)
    arr=ds28.ReadAsArray()
    arr=arr*28.
    fpath=os.path.join(month,os.path.basename(os.path.splitext(tif28)[0])+'_mmM'+'.tiff')
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        3600,
        1800,
        1,
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform(geotrans)
    tiff.GetRasterBand(1).WriteArray(arr)
    tiff.GetRasterBand(1).SetNoDataValue(-9999.9)
    tiff=None

# Function for leap year
def monthly29_avg(tif29,month):
    ds29=gdal.Open(tif29,gdal.GA_ReadOnly)
    arr=ds29.ReadAsArray()
    arr=arr*29.
    fpath=os.path.join(month,os.path.basename(os.path.splitext(tif29)[0])+'_mmM'+'.tiff')
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        3600,
        1800,
        1,
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform(geotrans)
    tiff.GetRasterBand(1).WriteArray(arr)
    tiff.GetRasterBand(1).SetNoDataValue(-9999.9)
    tiff=None

# For loop to define the files of certain months
for file in os.listdir(save):
    fname=os.fsdecode(file)
    if fname.endswith('.tiff'):
        for i in range(10):
            for j in range(10):
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,1):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,2):
                    tif28=os.path.join(save,fname)
                    monthly28_avg(tif28,month)
                    continue
                if fname.startswith('3B-MO.MS.MRG.3IMERG.20160201-S000000-E235959.02.V05B_B2_mmday.'):
                    tif29=os.path.join(save,fname)
                    monthly29_avg(tif29,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,3):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,4):
                    tif30=os.path.join(save,fname)
                    monthly30_avg(tif30,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,5):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,6):
                    tif30=os.path.join(save,fname)
                    monthly30_avg(tif30,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,7):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,8):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(0,9):
                    tif30=os.path.join(save,fname)
                    monthly30_avg(tif30,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(1,0):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(1,1):
                    tif30=os.path.join(save,fname)
                    monthly30_avg(tif30,month)
                    continue
                if fn.fnmatch(fname,'*.%s%s.V05B_B2_mmday.tiff'%(i,j)) is True and (i,j)==(1,2):
                    tif31=os.path.join(save,fname)
                    monthly31_avg(tif31,month)
                    continue
                else:
                    continue

#/////////////////////////////////////
# Test to ensure the leap year values are unique compared to other Feb values, and that Feb values did not get processed as leap

src='/IMERG/gpm_monthly/3B-MO.MS.MRG.3IMERG.20160201-S000000-E235959.02.V05B_B2_mmday_mmM.tiff'
src2='/IMERG/gpm_monthly/3B-MO.MS.MRG.3IMERG.20180201-S000000-E235959.02.V05B_B2_mmday_mmM.tiff'
srcc='/IMERG/gpm_daily/3B-MO.MS.MRG.3IMERG.20160201-S000000-E235959.02.V05B_B2_mmday.tiff'
srcc2='/IMERG/gpm_daily/3B-MO.MS.MRG.3IMERG.20180201-S000000-E235959.02.V05B_B2_mmday.tiff'

print('Leap Year - Feb 2016, mm/month')
def test(src):
    tif=gdal.Open(src,gdal.GA_ReadOnly)
    arr=tif.ReadAsArray()
    return arr[1400,1400]
print('#1 - Leap from for loop',test(src))
print()
print('Non-leap Year - Feb 2018, mm/month')
def test(src2):
    tif=gdal.Open(src2,gdal.GA_ReadOnly)
    arr=tif.ReadAsArray()
    return arr [1400,1400]
print('#2 - Non-leap from for loop',test(src2))
print()
def test(srcc):
    tif=gdal.Open(srcc,gdal.GA_ReadOnly)
    arr=tif.ReadAsArray()
    arr=arr*29.
    return arr[1400,1400]
print('#3 - Leap from direct calc',test(srcc))
print()
def test(srcc2):
    tif=gdal.Open(srcc2,gdal.GA_ReadOnly)
    arr=tif.ReadAsArray()
    arr=arr*29.
    return arr[1400,1400]
print('#4 - Non-leap from direct calc',test(srcc2),'multiplied Feb 2018 by 29, to show diff from #2')

#------------------------------------
# (3) Rename the files for seas60/100

for file in os.listdir(month):
    fname=os.fsdecode(file)
    if fname.endswith('.tiff'):
        path=os.path.join(month,fname)
        dst=os.path.join(seas,os.path.basename(os.path.splitext(path)[0])+'_seas60'+'.tiff')
        copyfile(path,dst)
        continue
    else:
        continue

for file in os.listdir(month):
    fname=os.fsdecode(file)
    if fname.endswith('.tiff'):
        path=os.path.join(month,fname)
        dst=os.path.join(seas,os.path.basename(os.path.splitext(path)[0])+'_seas100'+'.tiff')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# (4) Function to extract the months from respective years that are <= 60&100

def seas_60(seas60):
    with r.open(seas60,'r+') as dset:
        arr=dset.read()
        arr[arr==np.nan]=0
        arr[arr<=60.]=1
        arr[arr>60.]=0
        dset.write(arr)

def seas_100(seas100):
    with r.open(seas100,'r+') as dset:
        arr=dset.read()
        arr[arr==np.nan]=0
        arr[arr<=100.]=1
        arr[arr>100.]=0
        dset.write(arr)

#------------------------------------
# For loop to iterate over each file

for file in os.listdir(seas):
     fname=os.fsdecode(file)
     if fname.endswith('_seas60.tiff'):
         seas60=os.path.join(seas,fname)
         seas_60(seas60)
         continue
     else:
         continue

for file in os.listdir(seas):
     fname=os.fsdecode(file)
     if fname.endswith('_seas100.tiff'):
         seas100=os.path.join(seas,fname)
         seas_100(seas100)
         continue
     else:
         continue

#------------------------------------
# Copy each seasonality file type to its own folder and overlap

dir_60='/IMERG/seasonality/60'
dir_100='/IMERG/seasonality/100'
fpath='/IMERG/seasonality/3B-MO.MS.MRG.3IMERG.20140401-S000000-E235959.04.V05B_B2_mmday_mmM_seas60.tiff'
extent60=gdal.Open('/IMERG/seasonality/3B-MO.MS.MRG.3IMERG.20140401-S000000-E235959.04.V05B_B2_mmday_mmM_seas60.tiff',GA_ReadOnly)
ex60_band=extent60.GetRasterBand(1)
extent100=gdal.Open('/IMERG/seasonality/3B-MO.MS.MRG.3IMERG.20140401-S000000-E235959.04.V05B_B2_mmday_mmM_seas100.tiff',GA_ReadOnly)
ex100_band=extent100.GetRasterBand(1)
out_60='/IMERG/seasonality/Seasonality_60mm.tiff'
out_100='/IMERG/seasonality/Seasonality_100mm.tiff'

# 60mm file
for file in os.listdir(seas):
    fname=os.fsdecode(file)
    if fname.endswith('seas60.tiff'):
        path=os.path.join(seas,fname)
        dst=os.path.join(dir_60,fname)
        copyfile(path,dst)
        continue
    else:
        continue

crit_60='3*.tiff'
q60=os.path.join(dir_60,crit_60)
tile60=glob.glob(q60)
t60_arr=[]

for fp in tile60:
    fname=r.open(fp,'r+')
    arr=fname.read()
    NaNs=isnan(arr)
    arr[NaNs]=0
    arr=arr/4.25
    arr_2d=arr[0,:,:]
    t60_arr.append(arr_2d)
    fname.write(arr)

t60=np.array(t60_arr)
dataOut=t60.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_60,extent60.RasterXSize,extent60.RasterYSize,1,ex60_band.DataType)
CopyDatasetInfo(extent60,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

extent60=None
bandOut=None
dsOut=None

# 100mm file
for file in os.listdir(seas):
    fname=os.fsdecode(file)
    if fname.endswith('seas100.tiff'):
        path=os.path.join(seas,fname)
        dst=os.path.join(dir_100,fname)
        copyfile(path,dst)
        continue
    else:
        continue

crit_100='3*.tiff'
q100=os.path.join(dir_100,crit_100)
tile100=glob.glob(q100)
t100_arr=[]

for fp in tile100:
    fname=r.open(fp,'r+')
    arr=fname.read()
    NaNs=isnan(arr)
    arr[NaNs]=0
    arr=arr/4.25
    arr_2d=arr[0,:,:]
    t100_arr.append(arr_2d)
    fname.write(arr)

t100=np.array(t100_arr)
dataOut=t100.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_100,extent100.RasterXSize,extent100.RasterYSize,1,ex100_band.DataType)
CopyDatasetInfo(extent100,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

extent100=None
bandOut=None
dsOut=None

#------------------------------------
# (5) Folder and file paths
precip='/WorldClim-v2/precip_30s'
dir60='/WorldClim-v2/precip_30s/seas/seas60'
fpath='/WorldClim-v2/precip_30s/wc2_30s_prec_01.tif'

#------------------------------------
# Rename the files for seas60/100

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(dir60,os.path.basename(os.path.splitext(path)[0])+'_seas60'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Define function to extract the months from respective years that are <= 60&100

def seas_60(seas60):
    with r.open(seas60,'r+') as dset:
        arr=dset.read()
        #arr[arr==np.nan]=0 The IMERG datasets needed this, but WC appears to invert the desired
        arr[arr<=60.]=1
        arr[arr>60.]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate over each file

for file in os.listdir(dir60):
     fname=os.fsdecode(file)
     if fname.endswith('_seas60.tif'):
         seas60=os.path.join(dir60,fname)
         seas_60(seas60)
         continue
     else:
         continue

#------------------------------------
# Copy each seasonality file into its dir and overlap

ex60=gdal.Open('/WorldClim-v2/precip_30s/seas/seas60/wc2_30s_prec_01_seas60.tif',GA_ReadOnly)
b60=ex60.GetRasterBand(1)
out_60='/WorldClim-v2/precip_30s/seas/Seasonality_WorldClim_60mm.tif'

# 60mm file
crit_60='w*.tif'
q60=os.path.join(dir60,crit_60)
tile60=glob.glob(q60)
t60_arr=[]

for fp in tile60:
    fname=r.open(fp,'r+')
    arr=fname.read()
    NaNs=isnan(arr)
    arr[NaNs]=0
    arr_2d=arr[0,:,:]
    t60_arr.append(arr_2d)
    fname.write(arr)

t60=np.array(t60_arr)
dataOut=t60.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_60,ex60.RasterXSize,ex60.RasterYSize,1,b60.DataType)
CopyDatasetInfo(ex60,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex60=None
bandOut=None
dsOut=None

#------------------------------------
# Folder and file paths
precip='/WorldClim-v2/precip_30s'
dir100='/WorldClim-v2/precip_30s/seas/seas100'
fpath='/WorldClim-v2/precip_30s/wc2_30s_prec_01.tif'

#------------------------------------
# Rename the files for seas60/100

for file in os.listdir(precip):
    fname=os.fsdecode(file)
    if fname.endswith('.tif'):
        path=os.path.join(precip,fname)
        dst=os.path.join(dir100,os.path.basename(os.path.splitext(path)[0])+'_seas100'+'.tif')
        copyfile(path,dst)
        continue
    else:
        continue

#------------------------------------
# Define function to extract the months from respective years that are <= 60&100

def seas_100(seas100):
    with r.open(seas100,'r+') as dset:
        arr=dset.read()
        #arr[arr==np.nan]=0 The IMERG datasets needed this, but WC appears to invert the desired
        arr[arr<=100.]=1
        arr[arr>100.]=0
        dset.write(arr)

#------------------------------------
# Run a for loop to iterate over each file

for file in os.listdir(dir100):
     fname=os.fsdecode(file)
     if fname.endswith('_seas100.tif'):
         seas100=os.path.join(dir100,fname)
         seas_100(seas100)
         continue
     else:
         continue

#------------------------------------
# Copy each seasonality file into its dir and overlap

ex100=gdal.Open('/WorldClim-v2/precip_30s/seas/seas100/wc2_30s_prec_01_seas100.tif',GA_ReadOnly)
b100=ex100.GetRasterBand(1)
out_100='/WorldClim-v2/precip_30s/seas/Seasonality_WorldClim_100mm.tif'

# 100mm file
crit_100='w*.tif'
q100=os.path.join(dir100,crit_100)
tile100=glob.glob(q100)
t100_arr=[]

for fp in tile100:
    fname=r.open(fp,'r+')
    arr=fname.read()
    NaNs=isnan(arr)
    arr[NaNs]=0
    arr_2d=arr[0,:,:]
    t100_arr.append(arr_2d)
    fname.write(arr)

t100=np.array(t100_arr)
dataOut=t100.sum(axis=0)

driver=gdal.GetDriverByName('GTiff')
dsOut=driver.Create(out_100,ex100.RasterXSize,ex100.RasterYSize,1,b100.DataType)
CopyDatasetInfo(ex100,dsOut)
bandOut=dsOut.GetRasterBand(1)
BandWriteArray(bandOut,dataOut)

ex100=None
bandOut=None
dsOut=None
