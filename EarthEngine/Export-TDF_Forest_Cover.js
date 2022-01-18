/////////////////////////////////////////////////////////////////////////////////////
//
// Pantropical Forest Cover Estimates
// Using Tropical Dry Forest Definitions
// Jon Ocón, April 22, 2021 (UPDATED)
//
// Update line 40 to change canopy cover threshold, and update export descriptions.
//
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
//Areas of interest (.shp)//
/////////////////////////////////////////////////////////////////////////////////////

var wwf = ee.FeatureCollection("users/.../wwf_terr_ecos-TDF");
var bio = ee.FeatureCollection("users/.../bio_hotspots");

var tropics = ee.Geometry.Polygon(
        [[[-180.0, 30.0],
          [-180.0, -30.0],
          [180.0, -30.0],
          [180.0, 30.0]]],null,false);


/////////////////////////////////////////////////////////////////////////////////////
//Hansen Global Forest Cover//
/////////////////////////////////////////////////////////////////////////////////////


//Hansen GFC
var gfc = ee.Image("UMD/hansen/global_forest_change_2020_v1_8");

// Select the treecover, lossyear bands
var treecover = gfc.select(['treecover2000']);
var lossyear = gfc.select(['lossyear']);
var loss = lossyear.gte(0).and(lossyear.lte(20));

//Mask desired canopy threshold
var gfc2000 = treecover.gte(60).and(treecover.lte(100)); // >=40%
var loss2020 = gfc2000.updateMask(loss); //mask the forest covers


/////////////////////////////////////////////////////////////////////////////////////
//Clim Defs (.tif)//
/////////////////////////////////////////////////////////////////////////////////////


// Definitions
var FAO = ee.Image("users/.../FAO_CH");
var fao2000 = FAO.updateMask(gfc2000);
var fao2020 = FAO.updateMask(loss2020);

var ML = ee.Image("users/.../ML_WC");
var ml2000 = ML.updateMask(gfc2000);
var ml2020 = ML.updateMask(loss2020);

// Simple Definitions
var SW = ee.Image("users/.../Simple_WC");
var sw2000 = SW.updateMask(gfc2000);
var sw2020 = SW.updateMask(loss2020);

var SC = ee.Image("users/.../Simple_CH");
var sc2000 = SC.updateMask(gfc2000);
var sc2020 = SC.updateMask(loss2020);

// Consensus
var OW = ee.Image("users/.../Consensus-WC");
var ow2000 = OW.updateMask(gfc2000);
var ow2020 = OW.updateMask(loss2020);

var OC = ee.Image("users/.../Consensus-CH");
var oc2000 = OC.updateMask(gfc2000);
var oc2020 = OC.updateMask(loss2020);


/////////////////////////////////////////////////////////////////////////////////////
//Export Image//
/////////////////////////////////////////////////////////////////////////////////////


// FAO
var FAO_00 = fao2000.visualize({
  bands: ['b1'],
  max: 1
});
Export.image.toDrive({
  image: FAO_00,
  description: 'FAO_00',
  scale: 30,
  maxPixels: 1e13,
  region: tropics
});

var FAO_20 = fao2020.visualize({
  bands: ['b1'],
  max: 1
});
Export.image.toDrive({
  image: FAO_20,
  description: 'FAO_20',
  scale: 30,
  maxPixels: 1e13,
  region: tropics
});


/////////////////////////////////////////////////////////////////////////////////////
//Quantifying Forest Change 2000-2020//
/////////////////////////////////////////////////////////////////////////////////////


// World Wildlife Fund Tropical/Subtropical Broadleaf Deciduous Forest
var forestArea = gfc2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegions({
    reducer: ee.Reducer.sum(),
    collection: wwf,
    scale: 30
});

Export.table.toDrive({
  collection: forestSize,
  description: 'WWF_00',
  fileFormat: 'CSV'
});

var forestArea_loss = loss2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegions({
    reducer: ee.Reducer.sum(),
    collection: wwf,
    scale: 30
});

Export.table.toDrive({
  collection: forestSize_loss,
  description: 'WWF_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Biodiversity Hotspots
var forestArea = gfc2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegions({
    reducer: ee.Reducer.sum(),
    collection: bio,
    scale: 30
});

Export.table.toDrive({
  collection: forestSize,
  description: 'BIO_00',
  fileFormat: 'CSV'
});

var forestArea_loss = loss2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegions({
    reducer: ee.Reducer.sum(),
    collection: bio,
    scale: 30
});

Export.table.toDrive({
  collection: forestSize_loss,
  description: 'BIO_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// FAO CHELSA
var forestArea = fao2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'FAO_00',
  fileFormat: 'CSV'
});

var forestArea_loss = fao2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'FAO_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Murphy and Lugo Worldclim
var forestArea = ml2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'ML_00',  
  fileFormat: 'CSV'
});

var forestArea_loss = ml2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'ML_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Simple Worldlcim
var forestArea = sw2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'SW_00',
  fileFormat: 'CSV'
});

var forestArea_loss = sw2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'SW_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Simple CHELSA
var forestArea = sc2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'SC_00',
  fileFormat: 'CSV'
});

var forestArea_loss = sc2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'SC_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Worldlcim Consensus
var forestArea = ow2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'OW_00',
  fileFormat: 'CSV'
});

var forestArea_loss = ow2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'OW_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// CHELSA Consensus
var forestArea = oc2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'OC_00',
  fileFormat: 'CSV'
});

var forestArea_loss = oc2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'OC_20',
  fileFormat: 'CSV'
});


//////////////////////////////////////////////////////////////////////////////////////


// Pantropics
var forestArea = gfc2000.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize = forestArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize)
  ]),
  description: 'Pan_00',
  fileFormat: 'CSV'
});


var forestArea_loss = loss2020.multiply(ee.Image.pixelArea()).divide(1e6);
var forestSize_loss = forestArea_loss.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: tropics,
    scale: 30,
    maxPixels: 1e13
});

Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, forestSize_loss)
  ]),
  description: 'Pan_20',
  fileFormat: 'CSV'
});


