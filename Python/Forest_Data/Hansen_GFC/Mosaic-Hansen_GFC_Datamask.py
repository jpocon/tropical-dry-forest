#------------------------------------
#
# Mosaic Hansen Datamask files
#
#------------------------------------
# Set up the environment

import os
import rasterio
import numpy as np
from rasterio.merge import merge
from osgeo import gdal

#------------------------------------
# Define the variables
print('"mosaic.py" created by Jon Ocon, Aug 12, 2019.')
print('This program mosaics GFC Datamask tiles (Hansen et al 2013) for the pantropics.')
print()

# Directory variables
os.chdir('/datamask')
cwd=os.getcwd()

# Copy source tif metadata
tif=rasterio.open('Hansen_GFC-2018-v1.6_datamask_00N_000E.tif')

# Empty arrays for each section of tiles
src00N_E=[]
src10N_E=[]
src10S_E=[]
src20N_E=[]
src20S_E=[]
src30N_E=[]
src00N_W=[]
src10N_W=[]
src10S_W=[]
src20N_W=[]
src20S_W=[]
src30N_W=[] # Once code finishes, each array variable = None

#------------------------------------
# Merge rows of Lat (20S to 30N) by Lon (E or W)
print('Mosaic tifs by latitudes:')
print()

#------------------------------------
#00N - E (all following this section are identical minus lat/lon)

# For loop to pull the tiles matching the corresponding empty array
for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '00N' in fname and fname.endswith('E.tif'):
        src00N_E.append(rasterio.open(fname,'r',compress='LZW'))
        continue
    else:
        continue

print('Merging 00N_E tiles...')

# Merge the tiles
arr00N_E,out_trans=merge(src00N_E)
src00N_E=None

# Copy the metadata
print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 00N_E tif...')

# Write the merged tif to the working dir
with rasterio.open('merge00N_E.tif','w',**out_meta,compress='LZW') as dest:
    dest.write(arr00N_E)
    dest=None
    arr00N_E=None
    out_trans=None
    out_meta=None

print('merge00N_E.tif saved.')
print()

#------------------------------------
#00N - W

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '00N' in fname and fname.endswith('W.tif'):#
        src00N_W.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 00N_W tiles...')#

arr00N_W,out_trans=merge(src00N_W)#
src00N_W=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 00N_W tif...')#

with rasterio.open('merge00N_W.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr00N_W)#
    dest=None
    arr00N_W=None#
    out_trans=None
    out_meta=None

print('merge00N_W.tif saved.')#
print()

#------------------------------------
#10N - E (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '10N' in fname and fname.endswith('E.tif'):
        src10N_E.append(rasterio.open(fname,'r',compress='LZW'))
        continue
    else:
        continue

print('Merging 10N_E tiles...')

arr10N_E,out_trans=merge(src10N_E)
src10N_E=None

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 10N_E tif...')

with rasterio.open('merge10N_E.tif','w',**out_meta,compress='LZW') as dest:
    dest.write(arr10N_E)
    dest=None
    arr10N_E=None
    out_trans=None
    out_meta=None

print('merge10N_E.tif saved.')
print()

#------------------------------------
#10N - W (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '10N' in fname and fname.endswith('W.tif'):
        src10N_W.append(rasterio.open(fname,'r',compress='LZW'))
        continue
    else:
        continue

print('Merging 10N_W tiles...')

arr10N_W,out_trans=merge(src10N_W)
src10N_W=None

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 10N_W tif...')

with rasterio.open('merge10N_W.tif','w',**out_meta,compress='LZW') as dest:
    dest.write(arr10N_W)
    dest=None
    arr10N_W=None
    out_trans=None
    out_meta=None

print('merge10N_W.tif saved.')
print()

#------------------------------------
#10S - E (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '10S' in fname and fname.endswith('E.tif'):#
        src10S_E.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 10S_E tiles...')#

arr10S_E,out_trans=merge(src10S_E)#
src10S_E=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 10S_E tif...')#

with rasterio.open('merge10S_E.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr10S_E)#
    dest=None
    arr10S_E=None#
    out_trans=None
    out_meta=None

print('merge10S_E.tif saved.')#
print()

#------------------------------------
#10S - W (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '10S' in fname and fname.endswith('W.tif'):#
        src10S_W.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 10S_W tiles...')#

arr10S_W,out_trans=merge(src10S_W)#
src10S_W=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 10S_W tif...')#

with rasterio.open('merge10S_W.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr10S_W)#
    dest=None
    arr10S_W=None#
    out_trans=None
    out_meta=None

print('merge10S_W.tif saved.')#
print()

#------------------------------------
#20N - E (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '20N' in fname and fname.endswith('E.tif'):#
        src20N_E.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 20N_E tiles...')#

arr20N_E,out_trans=merge(src20N_E)#
src20N_E=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 20N_E tif...')#

with rasterio.open('merge20N_E.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr20N_E)#
    dest=None
    arr20N_E=None#
    out_trans=None
    out_meta=None

print('merge20N_E.tif saved.')#
print()

#------------------------------------
#20N - W (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '20N' in fname and fname.endswith('W.tif'):#
        src20N_W.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 20N_W tiles...')#

arr20N_W,out_trans=merge(src20N_W)#
src20N_W=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 20N_W tif...')#

with rasterio.open('merge20N_W.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr20N_W)#
    dest=None
    arr20N_W=None#
    out_trans=None
    out_meta=None

print('merge20N_W.tif saved.')#
print()

#------------------------------------
#20S - E (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '20S' in fname and fname.endswith('E.tif'):#
        src20S_E.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 20S_E tiles...')#

arr20S_E,out_trans=merge(src20S_E)#
src20S_E=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 20S_E tif...')#

with rasterio.open('merge20S_E.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr20S_E)#
    dest=None
    arr20S_E=None#
    out_trans=None
    out_meta=None

print('merge20S_E.tif saved.')#
print()

#------------------------------------
#20S - W (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '20S' in fname and fname.endswith('W.tif'):#
        src20S_W.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 20S_W tiles...')#

arr20S_W,out_trans=merge(src20S_W)#
src20S_W=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 20S_W tif...')#

with rasterio.open('merge20S_W.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr20S_W)#
    dest=None
    arr20S_W=None#
    out_trans=None
    out_meta=None

print('merge20S_W.tif saved.')#
print()

#------------------------------------
#30N - E (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '30N' in fname and fname.endswith('E.tif'):#
        src30N_E.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 30N_E tiles...')#

arr30N_E,out_trans=merge(src30N_E)#
src30N_E=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 30N_E tif...')#

with rasterio.open('merge30N_E.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr30N_E)#
    dest=None
    arr30N_E=None#
    out_trans=None
    out_meta=None

print('merge30N_E.tif saved.')#
print()

#------------------------------------
#30N - W (all following this section are identical minus lat/lon)

for file in os.listdir(cwd):
    fname=os.fsdecode(file)
    if '30N' in fname and fname.endswith('W.tif'):#
        src30N_W.append(rasterio.open(fname,'r',compress='LZW'))#
        continue
    else:
        continue

print('Merging 30N_W tiles...')#

arr30N_W,out_trans=merge(src30N_W)#
src30N_W=None#

print('Copying metadata...')

out_meta=tif.meta.copy()
out_meta.update({'width':720000,
                 'height':40000,
                 'transform':out_trans,
                 }
                )

print('Writing 30N_W tif...')#

with rasterio.open('merge30N_W.tif','w',**out_meta,compress='LZW') as dest:#
    dest.write(arr30N_W)#
    dest=None
    arr30N_W=None#
    out_trans=None
    out_meta=None

print('merge30N_W.tif saved.')#
print()
