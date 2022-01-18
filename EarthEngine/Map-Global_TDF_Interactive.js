
//////////////////////////////////////////////////////////////////////////////////////
//Tropics Clip//
/////////////////////////////////////////////////////////////////////////////////////


var tropics = ee.Geometry.Polygon(
        [[[-180.0, 30.0],
          [-180.0, -30.0],
          [180.0, -30.0],
          [180.0, 30.0]]],null,false);


//////////////////////////////////////////////////////////////////////////////////////
//WWF Ecoregions (Dry Forest)//
/////////////////////////////////////////////////////////////////////////////////////


//Dry Forest Ecoregions
var ecos = ee.FeatureCollection("users/.../wwf_terr_ecos-fix_geo-tdf");
Map.addLayer(ecos.draw({color:'FDF7CA'}),{},'Tropical Dry Forest Ecoregions');


//////////////////////////////////////////////////////////////////////////////////////
//Hansen Global Forest Cover//
/////////////////////////////////////////////////////////////////////////////////////


//Hansen GFC >5%
var gfc = ee.Image("UMD/hansen/global_forest_change_2017_v1_5");
var viz_gfc = {bands: 'treecover2000',
                min:5,
                max:100,
                palette:['0AFF37','08b228']};
// Select the land/water mask
var datamask = gfc.select('datamask');
var treecover = gfc.select('treecover2000');
// Create a binary mask
var mask = datamask.eq(1);
//Create an extent mask
var gfc_mask = treecover.gte(5);
var maskedGFC = gfc.updateMask(gfc_mask);
var mask_clip_GFC = maskedGFC.clip(tropics);
//Map.addLayer(mask_clip_GFC,viz_gfc,'Hansen GFC 2000 >5%');


//////////////////////////////////////////////////////////////////////////////////////
//Simard Canopy Height//
/////////////////////////////////////////////////////////////////////////////////////


//Tree canopy height 4.1-27.9m
var canopy_og = ee.Image("NASA/JPL/global_forest_canopy_height_2005");
var b1_canopy = canopy_og.select('1');
var viz_canopy = {bands: '1',
                  min:4.1,
                  max:27.9,
                  palette: ['f2ba7b','e08a28']};
var canopy_Lmask = b1_canopy.gte(4.1);
var canopy_Umask = b1_canopy.lte(27.9);
var Lmask_canopy = b1_canopy.updateMask(canopy_Lmask);
var ULmask_canopy = Lmask_canopy.updateMask(canopy_Umask);
var clip_canopy = ULmask_canopy.clip(tropics);
//Map.addLayer(clip_canopy,viz_canopy,'Simard Canopy Height 4.1-27.9m');


//////////////////////////////////////////////////////////////////////////////////////
//WorldClim v2 Temperature//
/////////////////////////////////////////////////////////////////////////////////////


//Temperature Variables
var MAT = ee.Image("users/.../wc2_bio/wc2_bio_30s_01"),
    MinT = ee.Image("users/.../wc2_bio/wc2_bio_30s_06"),
    MAP = ee.Image("users/.../wc2_bio/wc2_bio_30s_12");

//Visualization varibale
var viz_temp = {bands: 'b1',
                min:17,
                max:35,
                opacity:0.5,
                palette:['68E8C4','69D297']}; //Mint-green palette

//Mean Annual Temp >=17C
var b1_MAT = MAT.select('b1');
var MAT_mask = b1_MAT.gte(17);
var mask_MAT = b1_MAT.updateMask(MAT_mask);
var clip_MAT = mask_MAT.clip(tropics);
//Map.addLayer(clip_MAT,viz_temp,'WorldClim MAT >=17C');

//No freeze >0C
var b1_MinT = MinT.select('b1');
var viz_MinT = {bands:'b1',
              min:1,
              max:30,
              palette:['d58aea','c01eed']};
var MinT_mask = b1_MinT.gte(1);
var mask_MinT = b1_MinT.updateMask(MinT_mask);
var clip_MinT = mask_MinT.clip(tropics);
//Map.addLayer(clip_MinT,viz_temp,'WorldClim >0C');

//Temperature Mask
var T = clip_MinT.updateMask(clip_MAT);
Map.addLayer(T,viz_temp,'Temperature Mask >17C and >0C');


//////////////////////////////////////////////////////////////////////////////////////
//WorldClim v2 Precipitation//
/////////////////////////////////////////////////////////////////////////////////////


//MAP <1800mm (Dryflor)
var b1_MAP1800 = MAP.select('b1');
var viz_1800 = {bands: 'b1',
                min:0,
                max:1800,
                opacity:0.75,
                palette:['6DB6DE','6D7ADE','87B3FA','6B8EC7','6B8EC7']};
var MAP1800_mask = b1_MAP1800.lte(1800);
var mask_MAP1800 = b1_MAP1800.updateMask(MAP1800_mask);
var clip_MAP1800 = mask_MAP1800.clip(tropics);
Map.addLayer(clip_MAP1800,viz_1800,'WorldClim AP <1800 mm/yr');

//MAP 500-1500mm (FAO)
var b1_MAP1500 = MAP.select('b1');
var viz_1500 = {bands: 'b1',
                min:500,
                max:1500,
                opacity:0.75,
                palette:['6DB6DE','6D7ADE','87B3FA','6B8EC7','6B8EC7']};
var MAP1500_Umask = b1_MAP1500.lte(1500);
var MAP1500_Lmask = b1_MAP1500.gte(500);
var Lmask_MAP1500 = b1_MAP1500.updateMask(MAP1500_Lmask);
var ULmask_MAP1500 = Lmask_MAP1500.updateMask(MAP1500_Umask);
var clip_MAP1500 = ULmask_MAP1500.clip(tropics);
Map.addLayer(clip_MAP1500,viz_1500,'WorldClim AP 500-1500 mm/yr');

//MAP 250-2000mm (Murphy & Lugo)
var b1_MAP2000 = MAP.select('b1');
var viz_2000 = {bands: 'b1',
                min:250,
                max:2000,
                opacity:0.75,
                palette:['6DB6DE','6D7ADE','87B3FA','6B8EC7','6B8EC7']};
var MAP2000_Umask = b1_MAP2000.lte(2000);
var MAP2000_Lmask = b1_MAP2000.gte(250);
var Lmask_MAP2000 = b1_MAP2000.updateMask(MAP2000_Lmask);
var ULmask_MAP2000 = Lmask_MAP2000.updateMask(MAP2000_Umask);
var clip_MAP2000 = ULmask_MAP2000.clip(tropics);
Map.addLayer(clip_MAP2000,viz_2000,'WorldClim AP 250-2000 mm/yr');


//////////////////////////////////////////////////////////////////////////////////////
//GPM IMERG Seasonality//
/////////////////////////////////////////////////////////////////////////////////////


//Global Precipitation Measurement: IMERG
var i100 = ee.Image("users/.../Seasonality_100mm"),
    i60 = ee.Image("users/.../Seasonality_60mm");
var i100_b1 = i100.select('b1');
var viz_i100 = {bands: 'b1',
                min:0,
                max:100,
                opacity:0.75,
                palette:['B56BC2']};
var i60_b1 = i60.select('b1');
var viz_i60 = {bands: 'b1',
                min:0,
                max:60,
                opacity:0.75,
                palette:['B56BC2']};

//GPM <100mm 3-6 months/year
var i100_Lmask = i100_b1.gte(3);
var i100_Umask = i100_b1.lte(6);
var Lmask_i100 = i100_b1.updateMask(i100_Lmask);
var Umask_i100 = Lmask_i100.updateMask(i100_Umask);
var i100_tropics = Umask_i100.clip(tropics);
var i100_land = i100_tropics.updateMask(mask);
Map.addLayer(i100_land,viz_i100,'IMERG Seasonality <100mm for 3-6 mos/yr');

//GPM <60mm 4-7 months/year
var i60_Lmask = i60_b1.gte(4);
var i60_Umask = i60_b1.lte(7);
var Lmask_i60 = i60_b1.updateMask(i60_Lmask);
var Umask_i60 = Lmask_i60.updateMask(i60_Umask);
var i60_tropics = Umask_i60.clip(tropics);
var i60_land = i60_tropics.updateMask(mask);
Map.addLayer(i60_land,viz_i60,'IMERG Seasonality <60mm for 4-7 mos/yr');

//GPM <60mm 5-8 months/year
var i60_Lmask2 = i60_b1.gte(5);
var i60_Umask2 = i60_b1.lte(8);
var Lmask_i602 = i60_b1.updateMask(i60_Lmask2);
var Umask_i602 = Lmask_i602.updateMask(i60_Umask2);
var i60_tropics2 = Umask_i602.clip(tropics);
var i60_land2 = i60_tropics2.updateMask(mask);
Map.addLayer(i60_land2,viz_i60,'IMERG Seasonality <60mm for 5-8 mos/yr');


//////////////////////////////////////////////////////////////////////////////////////
//WorldClim v2 Seasonality//
/////////////////////////////////////////////////////////////////////////////////////


//WorldClim 2 - Dry Periods
var wc100 = ee.Image("users/.../Seasonality_WorldClim_100mm"),
    wc60 = ee.Image("users/.../Seasonality_WorldClim_60mm");
var wc100_b1 = wc100.select('b1');
var viz_wc100 = {bands: 'b1',
                min:0,
                max:100,
                opacity:0.75,
                palette:['B56BC2']};
var wc60_b1 = wc60.select('b1');
var viz_wc60 = {bands: 'b1',
                min:0,
                max:100,
                opacity:0.75,
                palette:['B56BC2']};

//WC <100mm 3-6 months/year
var wc100_Lmask = wc100_b1.gte(3);
var wc100_Umask = wc100_b1.lte(6);
var Lmask_wc100 = wc100_b1.updateMask(wc100_Lmask);
var Umask_wc100 = Lmask_wc100.updateMask(wc100_Umask);
var wc100_tropics = Umask_wc100.clip(tropics);
Map.addLayer(wc100_tropics,viz_wc100,'WorldClim Seasonality <100mm for 3-6 mos/yr');


//WC <60mm 4-7 months/year
var wc60_Lmask = wc60_b1.gte(4);
var wc60_Umask = wc60_b1.lte(7);
var Lmask_wc60 = wc60_b1.updateMask(wc60_Lmask);
var Umask_wc60 = Lmask_wc60.updateMask(wc60_Umask);
var wc60_tropics = Umask_wc60.clip(tropics);
Map.addLayer(wc60_tropics,viz_wc60,'WorldClim Seasonality <60mm for 4-7 mos/yr');


//WC <60mm 5-8 months/year
var wc60_Lmask2 = wc60_b1.gte(5);
var wc60_Umask2 = wc60_b1.lte(8);
var Lmask_wc602 = wc60_b1.updateMask(wc60_Lmask2);
var Umask_wc602 = Lmask_wc602.updateMask(wc60_Umask2);
var wc60_tropics2 = Umask_wc602.clip(tropics);
Map.addLayer(wc60_tropics2,viz_wc60,'WorldClim Seasonality <60mm for 5-8 mos/yr');


//////////////////////////////////////////////////////////////////////////////////////
//Defintions//
/////////////////////////////////////////////////////////////////////////////////////


//Murphy & Lugo (IMERG)//////////////////////////////////////////////////////////////
//Visualization
var viz_lugo = {bands: 'b1',
              min:250,
              max:2000,
              opacity:0.85,
              palette:['08B228','1BC902','5BBF02']};
//Temperature
var temp1 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip1 = clip_MAP2000.updateMask(i60_land);
//Maps
var lugo_clim = precip1.updateMask(temp1);
Map.addLayer(lugo_clim,viz_lugo,'Murphy & Lugo 250-2000mm, 4-7mos (IMERG)');

//Murphy & Lugo (WorldClim)//////////////////////////////////////////////////////////
//Temperature
var temp1 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip1 = clip_MAP2000.updateMask(wc60_tropics);
//Maps
var lugo_clim = precip1.updateMask(temp1);
Map.addLayer(lugo_clim,viz_lugo,'Murphy & Lugo 250-2000mm, 4-7mos (WorldClim)');


//DRYFLOR (IMERG)///////////////////////////////////////////////////////////////////
//Visualization
var viz_dflor = {bands: 'b1',
              min:0,
              max:1800,
              opacity:0.85,
              palette:['08B228','1BC902','5BBF02']};
//Temperature
var temp2 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip2 = clip_MAP1800.updateMask(i100_land);

//Maps
var dryflor_clim = precip2.updateMask(temp2);
Map.addLayer(dryflor_clim,viz_dflor,'DRYFLOR <1800mm, 3-6mos (IMERG)');

//DRYFLOR (WorldClim)///////////////////////////////////////////////////////////////
//Temperature
var temp2 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip2 = clip_MAP1800.updateMask(wc100_tropics);
//Maps
var dryflor_clim = precip2.updateMask(temp2);
Map.addLayer(dryflor_clim,viz_dflor,'DRYFLOR <1800mm, 3-6mos (WorldClim)');


//FAO (IMERG)///////////////////////////////////////////////////////////////////////
//Visualization
var viz_fao = {bands: 'b1',
              min:500,
              max:1500,
              opacity:0.85,
              palette:['08B228','1BC902','5BBF02']};
//Temperature
var temp3 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip3 = clip_MAP1500.updateMask(i60_land2);
//Maps
var fao_clim = precip3.updateMask(temp3);
Map.addLayer(fao_clim,viz_fao,'FAO 500-1500, 5-8mos (IMERG)');

//FAO (WorldClim)///////////////////////////////////////////////////////////////////
//Temperature
var temp3 = clip_MinT.updateMask(clip_MAT);
//Precipitation
var precip3 = clip_MAP1500.updateMask(wc60_tropics2);
//Maps
var fao_clim = precip3.updateMask(temp3);
Map.addLayer(fao_clim,viz_fao,'FAO 500-1500, 5-8mos (WorldClim)');


//////////////////////////////////////////////////////////////////////////////////////
//Known TDF Plots (Ibanez et al. 2018)//
/////////////////////////////////////////////////////////////////////////////////////


//TDF Plots
var plots = ee.FeatureCollection("users/.../tdf_pts");
Map.addLayer(plots.draw({color:'FF450A',pointRadius:5}),{},'Known TDF Plots');


//////////////////////////////////////////////////////////////////////////////////////
//Legend//
/////////////////////////////////////////////////////////////////////////////////////


//Position the panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
var legendTitle = ui.Label({
  value: 'Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
legend.add(legendTitle);
 
//Create and style 1 row of the legend
var makeRow = function(color, name) {
      //Create the label that is actually the colored box
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          //Use padding to give the box height and width
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
      //Create the label filled with the description text
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
      //Return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};

//Palette with the colors
var palette =['FDF7CA','FF450A','68E8C4','87B3FA','B56BC2','1BC902'];
 
//Name of the legend items
var names = ['TDF Ecoregions','TDF Plots','WorldClim Temp Mask','WorldClim Precip Masks','WorldClim/IMERG Seasonality Masks','TDF Extents (by definition)'];
 
//Add color and and names
for (var i = 0; i < 6; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
//Add legend to map
Map.add(legend);
