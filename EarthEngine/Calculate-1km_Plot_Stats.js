///////////////////////////////////////////////////////////////////////////////
//
// Extracting forest cover stats for 1 km plots
// J. P. Ocón, April 22, 2021
//
///////////////////////////////////////////////////////////////////////////////
//Variables & Functions//
///////////////////////////////////////////////////////////////////////////////

//Tropics bounds
var tropics = ee.Geometry.Polygon(
        [[[-180.0, 30.0],
          [-180.0, -30.0],
          [180.0, -30.0],
          [180.0, 30.0]]],null,false);

//540 TDF plots
var plots = ee.FeatureCollection('users/.../plots');
//print('Plots metadata:', plots);

//Buffer plots
function bufferPoints(radius, bounds) { //bounds: true or false (rect or circle)
  return function(pt) {
    pt = ee.Feature(pt);
    return bounds ? pt.buffer(radius).bounds() : pt.buffer(radius);
  };
}
var plots1km = plots.map(bufferPoints(1000, false)); //radius (m)
print('Buffer metadata:', plots1km);

//Climate definitions
var FAO = ee.Image('users/.../FAO_CH');
var FAO_wc = ee.Image('users/.../FAO_WC');
var ML = ee.Image('users/.../ML_CH');
var ML_wc = ee.Image('users/.../ML_WC');
var Dry = ee.Image('users/.../Dry_CH');
var Dry_wc = ee.Image('users/.../Dry_WC');
var AI = ee.Image('users/.../AI_CH');
var AI_wc = ee.Image('users/.../AI_WC');

//Global forest cover (Hansen et al. 2013)
var hansen = ee.Image('UMD/hansen/global_forest_change_2018_v1_6');
var gfc = hansen.reduceResolution({ //resample to 1 km
  reducer: ee.Reducer.mean(),
  maxPixels: 4096 //make big to avoid error
}).reproject(FAO.projection()); //reproject using FAO crs

///////////////////////////////////////////////////////////////////////////////
//Count Plots//
///////////////////////////////////////////////////////////////////////////////

//Which plots fall within each definition
var FAOplots = FAO.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: FAOplots,
  description: 'FAO_plots',
  fileFormat: 'CSV'
});

var fwc = FAO_wc.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: fwc,
  description: 'FAO_wc',
  fileFormat: 'CSV'
});

var ml = ML.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: ml,
  description: 'ML',
  fileFormat: 'CSV'
});

var mlwc = ML_wc.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: mlwc,
  description: 'ML_wc',
  fileFormat: 'CSV'
});

var dry = Dry.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: dry,
  description: 'Dry',
  fileFormat: 'CSV'
});

var dwc = Dry_wc.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: dwc,
  description: 'Dry_wc',
  fileFormat: 'CSV'
});

var ai = AI.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: ai,
  description: 'AI',
  fileFormat: 'CSV'
});

var awc = AI_wc.reduceRegions({
  reducer: ee.Reducer.count(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: awc,
  description: 'AI_wc',
  fileFormat: 'CSV'
});


///////////////////////////////////////////////////////////////////////////////
//Mean Canopy Cover//
///////////////////////////////////////////////////////////////////////////////

//Scale: 30 + Buffers
var treecover = hansen.select('treecover2000').clip(tropics);
print('treecover metadata:', treecover); //check metadata 30m scale
var mean30buffs = treecover.reduceRegions({
  reducer: ee.Reducer.mean(),
  collection: plots1km,
  scale: 30 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: mean30buffs,
  description: 'mean30buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Buffers
var treecover = gfc.select('treecover2000').clip(tropics);
print('treecover metadata:', treecover); //check metadata 1000m scale
var mean1000buffs = treecover.reduceRegions({
  reducer: ee.Reducer.mean(),
  collection: plots1km,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: mean1000buffs,
  description: 'mean1000buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Plots
var treecover = gfc.select('treecover2000').clip(tropics);
var mean1000plots = treecover.reduceRegions({
  reducer: ee.Reducer.mean(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: mean1000plots,
  description: 'mean1000plots',
  fileFormat: 'CSV'
});

///////////////////////////////////////////////////////////////////////////////
//Min Canopy Cover//
///////////////////////////////////////////////////////////////////////////////

//Scale: 30 + Buffers
var treecover = hansen.select('treecover2000').clip(tropics);
var min30buffs = treecover.reduceRegions({
  reducer: ee.Reducer.min(),
  collection: plots1km,
  scale: 30 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: min30buffs,
  description: 'min30buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Buffers
var treecover = gfc.select('treecover2000').clip(tropics);
var min1000buffs = treecover.reduceRegions({
  reducer: ee.Reducer.min(),
  collection: plots1km,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: min1000buffs,
  description: 'min1000buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Plots
var treecover = gfc.select('treecover2000').clip(tropics);
var min1000plots = treecover.reduceRegions({
  reducer: ee.Reducer.min(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: min1000plots,
  description: 'min1000plots',
  fileFormat: 'CSV'
});

///////////////////////////////////////////////////////////////////////////////
//Max Canopy Cover//
///////////////////////////////////////////////////////////////////////////////

//Scale: 30 + Buffers
var treecover = hansen.select('treecover2000').clip(tropics);
var max30buffs = treecover.reduceRegions({
  reducer: ee.Reducer.max(),
  collection: plots1km,
  scale: 30 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: max30buffs,
  description: 'max30buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Buffers
var treecover = gfc.select('treecover2000').clip(tropics);
var max1000buffs = treecover.reduceRegions({
  reducer: ee.Reducer.max(),
  collection: plots1km,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: max1000buffs,
  description: 'max1000buffs',
  fileFormat: 'CSV'
});

//Scale: 1000 + Plots
var treecover = gfc.select('treecover2000').clip(tropics);
var max1000plots = treecover.reduceRegions({
  reducer: ee.Reducer.max(),
  collection: plots,
  scale: 1000 //match img scale
});

Export.table.toDrive({ //workaround session timeout by exporting
  collection: max1000plots,
  description: 'max1000plots',
  fileFormat: 'CSV'
});

