#-------------------------------------------
#
# (1) Convert csv file to shp (points)(Need internet for EPSG WKT)
# (2) Buffer the TDF plots to .1 and .01-degrees
# (3) Rasterize the point buffers for analysis
#
#-------------------------------------------
# Set up the environment

import os
import osgeo
from osgeo import gdal,osr,ogr
import shapefile as shp
import csv
import numpy as np

# Set the working directory
os.chdir('')

#-------------------------------------------
# (1) Funtion to generate a .prj file
def getWKT_PRJ(epsg):
    from urllib.request import urlopen
    wkt = urlopen('http://spatialreference.org/ref/epsg/{0}/prettywkt/'.format(epsg))
    output = wkt.read()
    return output

# Create a point shapefile
tdf_shp = shp.Writer('tdf_buffer.shp')

# For every record there must be a corresponding geometry.
tdf_shp.autoBalance = 1

# Create the field names and data type for each.
tdf_shp.field('long','F')
tdf_shp.field('lat','F')
tdf_shp.field('MAP', 'F')
tdf_shp.field('presence', 'N')

# Count the features
counter = 1

# Access the CSV file
with open('tdf_pts.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    # Skip the header
    next(reader, None)
    # Loop through each of the rows and assign the attributes to variables
    for row in reader:
        long = row[0]
        lat = row[1]
        MAP = row[2]
        presence = row[3]
        # Create the point geometry
        tdf_shp.point(float(long),float(lat))
        # Add attribute data
        tdf_shp.record(long,lat,MAP,presence)
        print('Feature' + str(counter) + " added to Shapefile.")
        counter = counter + 1

# Save the Shapefile
tdf_shp.close()

# Create a projection file
prj = open('tdf.prj', 'w')
epsg = getWKT_PRJ('4326')
print(epsg)
prj.write('epsg')
prj.close()

#-------------------------------------------
# (2) Function to buffer points and retain attributes

def buffer_from_pts(inShpPath,bufferSize):
    buffShpPath=inShpPath.replace('.shp','_buffer-{}deg.shp'.format(bufferSize))
    print("\nCreating buffered shapefile...")
    inputds=ogr.Open(inShpPath)
    inputlyr=inputds.GetLayer()
    shpdriver=ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(buffShpPath):
        shpdriver.DeleteDataSource(buffShpPath)
    outputBufferds=shpdriver.CreateDataSource(buffShpPath)
    bufferlyr=outputBufferds.CreateLayer(buffShpPath,geom_type=ogr.wkbPolygon)
    featureDefn=bufferlyr.GetLayerDefn()
    # Create new fields in the output shp and get a list of field names for feature creation
    fieldNames=[]
    for i in range(inputlyr.GetLayerDefn().GetFieldCount()):
        fieldDefn=inputlyr.GetLayerDefn().GetFieldDefn(i)
        bufferlyr.CreateField(fieldDefn)
        fieldNames.append(fieldDefn.name)
    for feature in inputlyr:
        ingeom=feature.GetGeometryRef()
        fieldVals=[] # Make list of field values for feature
        for f in fieldNames: fieldVals.append(feature.GetField(f))
        outFeature=ogr.Feature(featureDefn)
        geomBuffer=ingeom.Buffer(bufferSize)
        outFeature.SetGeometry(geomBuffer)
        for v, val in enumerate(fieldVals): # Set output feature attributes
            outFeature.SetField(fieldNames[v],val)
        bufferlyr.CreateFeature(outFeature)
        outFeature=None
    # Copy the input .prj file
    from shutil import copyfile
    copyfile(inShpPath.replace('.shp','.prj'),buffShpPath.replace('.shp','.prj'))
    return buffShpPath

#-------------------------------------------
# Buffer .1-deg
buffer_from_pts('tdf_buffer.shp',.1)

# Buffer .01-deg
buffer_from_pts('tdf_buffer.shp',.01)

#-------------------------------------------
# (3) Pull and save folders
pull=''
save=''

# Function to rasterize shp and save as GTiff
def rasterize(shp,save):
    # Open the shp
    v_ds=gdal.OpenEx(shp,gdal.OF_VECTOR)
    lyr=v_ds.GetLayer()
    # Set pixel size
    pixW=pixH=.01
    # Set extent
    xmin,xmax,ymin,ymax=lyr.GetExtent()
    cols=int((xmax-xmin)/pixH)
    rows=int((ymax-ymin)/pixW)
    # Build save path
    fpath=os.path.join(save,os.path.basename(os.path.splitext(shp)[0])+'.tiff')
    # Write GTiff
    tiff=gdal.GetDriverByName('GTiff').Create(fpath,
        cols,
        rows,
        1,# Number of bands
        gdal.GDT_Float32,
        ['TILED=YES'])
    srs=osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    tiff.SetProjection(srs.ExportToWkt())
    tiff.SetGeoTransform((xmin,pixW,0,ymax,0,-pixH))
    tiff.GetRasterBand(1).SetNoDataValue(np.nan)
    gdal.Rasterize(tiff,v_ds)
    tiff=None # Close file and write to 'save'

# For loop to iterate over each buffer
for file in os.listdir(pull):
     fname=os.fsdecode(file)
     if fname.endswith('deg.shp'):
         shp=os.path.join(pull,fname)
         rasterize(shp,save)
         continue
     else:
         continue
