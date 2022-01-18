#--------------------
#
#(1)Create 7200x43200 array filled with desired areas for deg to km2
#   Every row, values increase from 742249591 to 85734000 in 3191 increments then back to 742249591
#(2)Sum both arrays and all values that == 742249591 ... 85734000 = 0, convert to int32
#(3)Clip global rasters to each shapefile
#(4)Total all remaining values in each file and print to xlsx
#
#--------------------
#Set up the environment
print('Importing necessary modules...')
#Files
import os,sys
import fnmatch as fn
from shutil import copyfile
import glob
import xlsxwriter

#Maths
import numpy as np
from numpy import *

#GIS/Raster
from PIL import ImageDraw,Image
import fiona
import rasterio
import rasterio.mask
from rasterstats import zonal_stats
from osgeo import gdal,osr,ogr,gdalnumeric
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
gdal.UseExceptions()
print('Modules successfully imported.')
print()

#--------------------
#(1) Create an array matching 30N/S clip full of accurate area values for each row
print('Creating 30N/S array for deg to km2 area conversion...')
#Fill an array with desired area values
a=np.arange(74249591,85737191,3191)
b=np.arange(85734000,74246400,-3191)
c=np.append(a,b)
d=np.array([c,c])
e=np.transpose(d)

#Duplicate the second column 43198 times for 7200x43200 array
def dup_cols(e,indx,num_dups=1):
    return np.insert(e,[indx+1]*num_dups,e[:,[indx]],axis=1)

f=dup_cols(e,indx=1,num_dups=43198)
print('Array created.')
print()

#--------------------
#Folder paths 'C:/Users/Jonathan Ocon/Desktop/' if on laptop, 'D:/' if using SSD
print('Defining all relevant folder paths...')
#Climate/land/AI raster files
rast=''
land=''

#Regional shapefiles
htspt=''
cntry=''

#Final outputs
FINAL=''
print('Folders defined.')
print()
#--------------------
#(1) Call in each raster and sum with the custom array above, any values == original custom vals = 0
'''
def km(area):
    with rasterio.open(area,'r+',compress='LZW') as dset:
        arr=dset.read()
        arr=arr.astype(rasterio.float64)
        arr=arr+f
        arr[arr==f]=0
        with rasterio.Env():
            profile=dset.profile
            profile.update(
                dtype=rasterio.float64)
            with rasterio.open(area,'w',**profile) as dset:
                dset.write(arr)

print('Converting deg areas to km2...')
for file in os.listdir(rast):
    fname=os.fsdecode(file)
    if fname.endswith('-area.tif'):
        area=os.path.join(rast,fname)
        print('Converting '+area+' ...')
        km(area)
        print(area+' converted.')
        print()
        continue
    else:
        continue
print('Deg areas converted to km2.')
print()
'''
#--------------------
#(3) Call in each raster, clip to each shapefile
#Raster criteria
crit_rast='*-area.tif'
q_rast=os.path.join(rast,crit_rast)
tile_rast=glob.glob(q_rast)

#COUNTRIES
crit_cntry='*.shp'
q_cntry=os.path.join(cntry,crit_cntry)
tile_cntry=glob.glob(q_cntry)

def clip_cntry(raster):
    for fp in tile_cntry:
        if fp.endswith('.shp'):
            shp=fp
            print('Clipping '+raster+' to '+shp+' ...')
            with fiona.open(shp,'r') as shapefile:
                features=[feature['geometry'] for feature in shapefile]
            with rasterio.open(raster) as src:
                out_image,out_transform=rasterio.mask.mask(src,features,
                                                    crop=True)
                out_meta=src.meta.copy()
                out_meta.update({'driver':'GTiff',
                                 'height':out_image.shape[1],
                                 'width':out_image.shape[2],
                                 'transform':out_transform})
                with rasterio.open(os.path.join(FINAL,os.path.basename(os.path.splitext(raster)[0])+
                                                '_'+os.path.basename(os.path.splitext(shp)[0])+
                                                '_FINAL_cntry'+'.tif'),'w',**out_meta) as dest:
                    dest.write(out_image)
            print(raster+' clipped.')
            print()
            continue
        else:
            continue

#HOTSPOTS
crit_htspt='*.shp'
q_htspt=os.path.join(htspt,crit_htspt)
tile_htspt=glob.glob(q_htspt)

def clip_htspt(raster):
    for fp in tile_htspt:
        if fp.endswith('.shp'):
            shp=fp
            print('Clipping '+raster+' to '+shp+' ...')
            with fiona.open(shp,'r') as shapefile:
                features=[feature['geometry'] for feature in shapefile]
            with rasterio.open(raster) as src:
                out_image,out_transform=rasterio.mask.mask(src,features,
                                                    crop=True)
                out_meta=src.meta.copy()
                out_meta.update({'driver':'GTiff',
                                 'height':out_image.shape[1],
                                 'width':out_image.shape[2],
                                 'transform':out_transform})
                with rasterio.open(os.path.join(FINAL,os.path.basename(os.path.splitext(raster)[0])+
                                                '_'+os.path.basename(os.path.splitext(shp)[0])+
                                                '_FINAL_htspt'+'.tif'),'w',**out_meta) as dest:
                    dest.write(out_image)
            print(raster+' clipped.')
            print()
            continue
        else:
            continue

#--------------------
print('Clipping climate rasters to shapefiles...')
for fp in tile_rast:
    if fp.endswith('.tif'):
        raster=fp
        print('Clipping '+raster+' to countries...')
        clip_cntry(raster)
        print(raster+' countries clipped.')
        print()
        print('Clipping '+raster+' to hotspots...')
        clip_htspt(raster)
        print(raster+' hotspots clipped.')
        print()
        continue
    else:
        continue
print('Clipping climate rasters to shapefiles complete.')
print()

#--------------------
#(4) Run raster stats on each file to determine the area and output each result to excel file

crit_final='*.tif'
q_final=os.path.join(FINAL,crit_final)
tile_final=glob.glob(q_final)
tfinal_arr=[]

def stats(fp):
    with rasterio.open(fp,'r+',compress='LZW') as src:
        arr=src.read(1,masked=True) # covers up nan vals
        arr=arr/1e8
        if np.sum(arr)>=0:
            tfinal_arr.append([fp,np.sum(arr)])
            src=None

for fp in tile_final:
    if fp.endswith('.tif'):
        print('Summing the areas for '+fp+' ...')
        stats(fp)
        print('Areas summed.')
        print()
        continue
    else:
        continue

#--------------------------------------------------
# Define the variables and cwd

print('Setting working directory to desired folder.')
os.chdir('')
cwd=os.getcwd()

print('Compiling Excel sheet...')
print()
wkbk=xlsxwriter.Workbook('tdf_cntry_htspt_areas.xlsx')
wksht=wkbk.add_worksheet()

bold=wkbk.add_format({'bold':True})
num=wkbk.add_format({'num_format':'#,##0.00'})

wksht.write('A1','Region',bold)
wksht.write('B1','Area',bold)

row=1
col=0

for region, area in (tfinal_arr):
    wksht.write(row,col,region)
    wksht.write(row,col+1,area,num)
    row+=1
    
wkbk.close()
print('Excel sheet saved.')
print()

#--------------------------------------------------
# Function to compress and convert

def compress(fp):
    print('Calling in '+fp+' and converting nodataval to 65535 and dtype to uint8...')
    print()
    src=rasterio.open(fp,'r+')
    arr=src.read()
    arr[arr>0]=1
    arr[arr==-3.4e+38]=2
    with rasterio.open(fp,'r+',compress='LZW') as dset:
        dset.nodata=2
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
        dest=None
        arr=None
        src=None
        meta=None

#--------------------------------------------------
# Append all full path filenames to array to iterate

print('Compressing all tifs in current directory...')
print()

for file in os.listdir(FINAL):
     fname=os.fsdecode(file)
     if fname.endswith('.tif'):
        fp=os.path.join(FINAL,file)
        compress(fp)

print('The program in done running!')
