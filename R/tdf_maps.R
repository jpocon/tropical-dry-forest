#set up the env####
library(easypackages)

#data processing
libraries("tidyverse", "dplyr", "tidyr", "magrittr", "png")
#spatial data processing
libraries("raster", "rgdal", "sf", "lwgeom", "spData", "foreign", "maptools", "gdalUtils", "rgeos")
#geovisualization
libraries("tmap", "tmaptools", "ggplot2", "gridExtra", "rasterVis", "grid") #static
libraries("mapview", "leaflet", "shiny", "dismo") #interactive

setwd("")
cwd <- getwd()
windowsFonts(Times = windowsFont("TT Times")) 

#polygonize rasters into shapefile - ONLY RUN ONCE####
gdal_polygonizeR <- function(x, outshape = NULL, gdalformat = "ESRI Shapefile",
                             pypath = NULL, readpoly = T, quiet = T) {
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which("gdal_polygonize.py")
  }
  if (!file.exists(pypath)) stop("gdal_polygonize.py not on system!")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub("\\.shp$", "", outshape)
    f.exists <- file.exists(paste(outshape, c("shp", "shx", "dbf"), sep = "."))
    if (any(f.exists))
      stop(sprintf("File already exists: %s",
                   toString(paste(outshape, c("shp", "shx", "dbf"),
                                  sep = '.')[f.exists])), call. = F)
  } else outshape <- tempfile()
  if (is(x, "Raster")) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext = ".asc")})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop("'x' must be a file path (char str), or a raster object.")
  system2("python", args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                  pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}

land <- gdal_polygonizeR(file.path(cwd, "data/baselayers/land_tropics.tif"), readpoly = T) %>%
  gBuffer(byid = T, width = 0) #fix geoms
#tropics
poly <- readWKT("POLYGON((
                -180.0 -30.0,
                180.0 -30.0,
                180.0 30.0,
                -180.0 30.0,
                -180.0 -30.0
))", p4s = CRS("+proj=longlat +datum=WGS84"))
land <- crop_shape(land, poly)
#rasters
climate <- list.files("data/climate/proc_clim/", pattern = ".tif", full.names = T)
#AI
wc_65 <- gdal_polygonizeR(file.path(cwd, climate[11]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
ch_65 <- gdal_polygonizeR(file.path(cwd, climate[1]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
#Worldclim
ml_wc <- gdal_polygonizeR(file.path(cwd, climate[7]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
fao_wc <- gdal_polygonizeR(file.path(cwd, climate[5]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
dry_wc <- gdal_polygonizeR(file.path(cwd, climate[3]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
#CHELSA
ml_ch <- gdal_polygonizeR(file.path(cwd, climate[6]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
fao_ch <- gdal_polygonizeR(file.path(cwd, climate[4]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)
dry_ch <- gdal_polygonizeR(file.path(cwd, climate[2]), readpoly = T) %>%
  gBuffer(byid = T, width = 0)




#call in shapefiles####
subcons <- list.files("data/baselayers/proc_regions/macro/mapping/",
                      pattern = ".shp", full.names = T)
hotspots_tdf <- list.files("data/baselayers/proc_regions/meso/mapping",
                           pattern = ".shp", full.names = T)
archs <- list.files("data/baselayers/proc_regions/micro/mapping",
                    pattern = ".shp", full.names = T)
#baselayers
countries <- readOGR("data/baselayers/raw_regions/features/ne_10m_admin_0_countries/ne_10m_admin_0_countries-tropics-plots.shp") %>%
  gBuffer(byid = T, width = 0)
ecoregions <- readOGR("data/ecoregions/proc_ecos/wwf_terr_ecos-TDF.shp") %>%
  gBuffer(byid = T, width = 0)
plots <- readOGR("data/plots/TDF_plots-MASTER-091919.shp")




#mapping####
tmap_mode("plot")
abc <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.")




#PLOTS####
tdf_plots <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(ecoregions) + #ecoregions
  tm_fill(col = "#758000", 
          alpha = 0.8) +
  
  tm_shape(plots) + #TDF plots
  tm_dots(size = 0.1,
          col = "#CC0000",
          legend.z = 0) +
  
  tm_add_legend("symbol", col = "#758000", 
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "WWF TSBF Ecoregions") +
  tm_add_legend("symbol", col = "#CC0000", 
                shape = 16, size = 0.75,
                labels = "Tropical Dry Forest Plots") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")
  

tmap_save(tdf_plots, "maps/tdf_plots.png", 
          width = 7, units = "in", asp = 0)




#AI####
tdf_wc_65 <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(wc_65) + #Worldclim AI
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "Worldclim v2",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "AI <= 0.65") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "a.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_ch_65 <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(ch_65) + #CHLESA AI
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "CHELSA",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "AI <= 0.65") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "b.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_overlap_ai <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(raster(climate[8])) + #Consensus AI
  tm_raster(palette = c("#3DBC74", "#33638D"),
            legend.show = F) +
  
  tm_add_legend("symbol", col = c("#3DBC74", "#33638D"), title = "AI Comparison",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus")) +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "c.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


list_ai <- list(tdf_wc_65, tdf_ch_65, tdf_overlap_ai) #list ai maps
tdf_ai <- tmap_arrange(list_ai, ncol = 1, nrow = 3, 
                       asp = NA, outer.margins = NULL) #stack/arrange maps vertically


tmap_save(tdf_ai, "maps/tdf_ai.png", width = 7, units = "in")




#WORLDCLIM####
tdf_wc_ml <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(ml_wc) + #Worldclim Murphy and Lugo
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "Worldclim v2",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "Murphy and Lugo") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "a.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_wc_fao <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(fao_wc) + #Worldclim FAO
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "Worldclim v2",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "FAO") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "b.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_wc_dry <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(dry_wc) + #Worldclim Dryflor
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "Worldclim v2",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "Dryflor") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "c.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_overlap_wc <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(raster(climate[10])) + #Consensus WC
  tm_raster(palette = c("#33638D", "#3DBC74", "#CAE11F"),
            legend.show = F) +
  
  tm_add_legend("symbol", col = c("#33638D", "#3DBC74", "#CAE11F"), title = "Worldclim Comparison",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus Among 2", "Consensus Among All")) +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "d.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


list_wc <- list(tdf_wc_ml, tdf_wc_fao, tdf_wc_dry, tdf_overlap_wc) #list wc maps
tdf_wc <- tmap_arrange(list_wc, ncol = 1, nrow = 4, 
                       asp = NA, outer.margins = NULL) #stack/arrange maps vertically


tmap_save(tdf_wc, "maps/tdf_wc.png", width = 7, units = "in")




#CHELSA####
tdf_ch_ml <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(ml_ch) + #CHELSA Murphy and Lugo
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "CHELSA",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "Murphy and Lugo") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "a.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_ch_fao <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(fao_ch) + #CHELSA FAO
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "CHELSA",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "FAO") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "b.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_ch_dry <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(dry_ch) + #CHELSA Dryflor
  tm_fill(col = "#33638D") +
  
  tm_add_legend("symbol", col = "#33638D", title = "CHELSA",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = "Dryflor") +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "c.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


tdf_overlap_ch <- tm_shape(land) + #basemap
  tm_fill(col = "#FFFFFF") + 
  
  tm_shape(countries) + #countries
  tm_borders(col = "#CCCCCC") +
  
  tm_shape(raster(climate[9])) + #Consensus CH
  tm_raster(palette = c("#33638D", "#3DBC74", "#CAE11F"),
            legend.show = F) +
  
  tm_add_legend("symbol", col = c("#33638D", "#3DBC74", "#CAE11F"), title = "CHELSA Comparison",
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus Among 2", "Consensus Among All")) +
  
  tm_layout(bg.color = "#A6CEE3", #ocean
            main.title = "d.",
            main.title.position = 0,
            outer.margins = c(0, 0, 0, 0),
            legend.position = c(0.01, 0.05), #coords (left, bottom)
            legend.title.size = 0.8,
            legend.text.size = 0.5, #sizing helps stack text
            legend.width = -1, #prevents any auto re-sizing
            fontfamily = "Times")


list_ch <- list(tdf_ch_ml, tdf_ch_fao, tdf_ch_dry, tdf_overlap_ch) #list ch maps
tdf_ch <- tmap_arrange(list_ch, ncol = 1, nrow = 4, 
                       asp = NA, outer.margins = NULL) #stack/arrange maps vertically


tmap_save(tdf_ch, "maps/tdf_ch.png", width = 7, units = "in")


#Subcontinents####
for (i in 1:5){
  subcon <- readOGR(subcons[i]) %>% gBuffer(byid = T, width = 0)
  
  if (i == 3){
    
    
    PACsubcon <- spTransform(subcon, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    PACcountries <- spTransform(countries, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    PACplots <- spTransform(plots, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
    
    PACml_ch <- spTransform(ml_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    PACfao_ch <- spTransform(fao_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    PACdry_ch <- spTransform(dry_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    
    
    PACsubcon_ml_ch <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACml_ch, PACsubcon)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    PACsubcon_fao_ch <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACfao_ch, PACsubcon)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    PACsubcon_dry_ch <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACdry_ch, PACsubcon)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #inset of Hawaii
    hawaii <- readOGR(archs[3]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    PACH_ml_ch <-  tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACml_ch, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    PACH_fao_ch <- tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACfao_ch, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    PACH_dry_ch <- tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACdry_ch, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    #save maps with insets
    tmap_save(PACsubcon_ml_ch, paste("maps/pac_ch_1.png", sep = ""), 
              insets_tm = PACH_ml_ch, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(PACsubcon_fao_ch, paste("maps/pac_ch_2.png", sep = ""), 
              insets_tm = PACH_fao_ch, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(PACsubcon_dry_ch, paste("maps/pac_ch_3.png", sep = ""), 
              insets_tm = PACH_dry_ch, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "pac_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    pac_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig6_", i, ".png", sep = ""), 
           pac_ch, width = 7, height = 2, units = "in")
    
    
  } else {
    
    
    subcon_ml_ch <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, subcon)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    subcon_fao_ch <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, subcon)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    subcon_dry_ch <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, subcon)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(subcon_ml_ch, paste("maps/sub_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(subcon_fao_ch, paste("maps/sub_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(subcon_dry_ch, paste("maps/sub_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "sub_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    sub_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig6_", i, ".png", sep = ""), 
           sub_ch, width = 7, height = 2, units = "in")
    
    
  }
} #CHELSA


for (i in 1:5){
  subcon <- readOGR(subcons[i]) %>% gBuffer(byid = T, width = 0)
  
  if (i == 3){
    
    
    PACsubcon <- spTransform(subcon, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    PACcountries <- spTransform(countries, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    PACplots <- spTransform(plots, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
    
    PACml_wc <- spTransform(ml_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    PACfao_wc <- spTransform(fao_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    PACdry_wc <- spTransform(dry_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>% gBuffer(byid = T, width = 0)
    
    
    PACsubcon_ml_wc <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACml_wc, PACsubcon)) + #ML Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    PACsubcon_fao_wc <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACfao_wc, PACsubcon)) + #FAO Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    PACsubcon_dry_wc <- tm_shape(PACsubcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(PACcountries, PACsubcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACdry_wc, PACsubcon)) + #Dryflor Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, PACsubcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Wordclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #inset of Hawaii
    hawaii <- readOGR(archs[3]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    PACH_ml_wc <-  tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACml_wc, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    PACH_fao_wc <- tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACfao_wc, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    PACH_dry_wc <- tm_shape(hawaii) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(PACcountries, hawaii)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(PACdry_wc, hawaii)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(PACplots, hawaii)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    #save maps with insets
    tmap_save(PACsubcon_ml_wc, paste("maps/pac_wc_1.png", sep = ""), 
              insets_tm = PACH_ml_wc, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(PACsubcon_fao_wc, paste("maps/pac_wc_2.png", sep = ""), 
              insets_tm = PACH_fao_wc, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(PACsubcon_dry_wc, paste("maps/pac_wc_3.png", sep = ""), 
              insets_tm = PACH_dry_wc, insets_vp = viewport(x = 0.75, y = 0.65, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "pac_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    pac_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig7_", i, ".png", sep = ""), 
           pac_wc, width = 7, height = 2, units = "in")
    
    
  } else {
    
    
    subcon_ml_wc <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, subcon)) + #ML Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    subcon_fao_wc <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_wc, subcon)) + #FAO Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    subcon_dry_wc <- tm_shape(subcon) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, subcon)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, subcon)) + #Dryflor Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, subcon)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(subcon_ml_wc, paste("maps/sub_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(subcon_fao_wc, paste("maps/sub_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(subcon_dry_wc, paste("maps/sub_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "sub_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    sub_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig7_", i, ".png", sep = ""), 
           sub_wc, width = 7, height = 2, units = "in")
    
    
  }
} #Worldclim




#Biodiversity Hotspots####
for (i in 1:7){
  
  
  if (i == 4){
    
    
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, hot)) + #ML CHELSA
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, hot)) + #FAO CHELSA
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, hot)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(hot_ml_ch, paste("maps/hot_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_ch, paste("maps/hot_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_ch, paste("maps/hot_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "hot_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    hot_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig8_", i, ".png", sep = ""), 
           hot_ch, width = 7, height = 2, units = "in")
    
    
  } else if (i == 5){
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, hot)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, hot)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, hot)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #inset of Hawaii
    gala <- readOGR(archs[2]) %>% gBuffer(byid = T, width = 0)
    
    
    gala_ml_ch <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, gala)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    gala_fao_ch <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, gala)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    gala_dry_ch <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") + #no Dryflor for Galapagos
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    #save maps with insets
    tmap_save(hot_ml_ch, paste("maps/gala_ch_1.png", sep = ""), 
              insets_tm = gala_ml_ch, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_ch, paste("maps/gala_ch_2.png", sep = ""), 
              insets_tm = gala_fao_ch, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_ch, paste("maps/gala_ch_3.png", sep = ""), 
              insets_tm = gala_dry_ch, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "gala_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    gala_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig8_", i, ".png", sep = ""), 
           gala_ch, width = 7, height = 2, units = "in")
    
    
  } else if (i == 7){
    
    
    poly_plots <- spTransform(plots, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    poly_countries <- spTransform(countries, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_ml_ch <- spTransform(ml_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_fao_ch <- spTransform(fao_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_dry_ch <- spTransform(dry_ch, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    
    
    #hawaii
    haw <- readOGR(archs[3]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    haw_ml_ch <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_ch, haw)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "Murphy and Lugo"), 
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_ml_ch, paste("maps/poly/a_ch_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    haw_fao_ch <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_fao_ch, haw)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "FAO"),
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_fao_ch, paste("maps/poly/a_ch_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    haw_dry_ch <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_ch, haw)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "Dryflor"),
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_dry_ch, paste("maps/poly/a_ch_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #marquesas
    marq <- readOGR("data/baselayers/proc_regions/meso/mapping/poly/marquesas.shp") %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    marq_ml_ch <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_ch, marq)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_ml_ch, paste("maps/poly/b_ch_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    marq_fao_ch <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_fao_ch, marq)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_fao_ch, paste("maps/poly/b_ch_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    marq_dry_ch <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_ch, marq)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_dry_ch, paste("maps/poly/b_ch_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #fiji
    fiji <- readOGR(archs[1]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    fiji_ml_ch <-  tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_ch, fiji)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_ml_ch, paste("maps/poly/c_ch_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    fiji_fao_ch <- tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") + #no FAO in Fiji
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_fao_ch, paste("maps/poly/c_ch_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    fiji_dry_ch <- tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_ch, fiji)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_dry_ch, paste("maps/poly/c_ch_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #northern mariana
    nmi <- readOGR("data/baselayers/proc_regions/meso/mapping/poly/northern_mariana.shp") %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    nmi_ml_ch <-  tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_ml_ch, nmi)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_ml_ch, paste("maps/poly/d_ch_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    nmi_fao_ch <- tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") + #no FAO overlap NMI
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_fao_ch, paste("maps/poly/d_ch_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    nmi_dry_ch <- tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_dry_ch, nmi)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_dry_ch, paste("maps/poly/d_ch_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #ML
    pnglist = lapply(list.files("maps/poly", pattern = "*_ch_ml.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_ml <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_1.png", poly_ml, width = 7/3*2, height = 7/3*2, units = "in")
    
    #FAO
    pnglist = lapply(list.files("maps/poly", pattern = "*_ch_fao.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_fao <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_2.png", poly_fao, width = 7/3*2, height = 7/3*2, units = "in")
    
    #Dryflor
    pnglist = lapply(list.files("maps/poly", pattern = "*_ch_dry.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_dry <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_3.png", poly_dry, width = 7/3*2, height = 7/3*2, units = "in")
    
    
  } else {
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, hot)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, hot)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_ch <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, hot)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(hot_ml_ch, paste("maps/hot_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_ch, paste("maps/hot_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_ch, paste("maps/hot_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "hot_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    hot_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig8_", i, ".png", sep = ""), 
           hot_ch, width = 7, height = 2, units = "in")
    
    
  }
} #CHELSA


for (i in 1:7){
  
  
  if (i == 4){
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, hot)) + #ML WC
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_wc, hot)) + #FAO WC
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, hot)) + #Dryflor WC
      tm_fill(col = "#33638D") + #no plots in Sundaland
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(hot_ml_wc, paste("maps/hot_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_wc, paste("maps/hot_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_wc, paste("maps/hot_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "hot_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    hot_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig9_", i, ".png", sep = ""), 
           hot_wc, width = 7, height = 2, units = "in")
    
    
  } else if (i == 5){
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, hot)) + #ML WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_wc, hot)) + #FAO WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, hot)) + #Dryflor WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.05, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #inset of Hawaii
    gala <- readOGR(archs[2]) %>% gBuffer(byid = T, width = 0)
    
    
    gala_ml_wc <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC for Galapagos
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    gala_fao_wc <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC for Galapagos
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    gala_dry_wc <-  tm_shape(gala) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(countries, gala)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC for Galapagos
      tm_shape(raster::crop(plots, gala)) +
      tm_symbols(size = 0.01, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3")
    
    
    #save maps with insets
    tmap_save(hot_ml_wc, paste("maps/gala_wc_1.png", sep = ""), 
              insets_tm = gala_ml_wc, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_wc, paste("maps/gala_wc_2.png", sep = ""), 
              insets_tm = gala_fao_wc, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_wc, paste("maps/gala_wc_3.png", sep = ""), 
              insets_tm = gala_dry_wc, insets_vp = viewport(x = 0.25, y = 0.5, width = 0.4, height = 0.4), 
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "gala_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    gala_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig9_", i, ".png", sep = ""), 
           gala_wc, width = 7, height = 2, units = "in")
    
    
  } else if (i == 7){
    
    
    poly_plots <- spTransform(plots, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    poly_countries <- spTransform(countries, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_ml_wc <- spTransform(ml_wc, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_fao_wc <- spTransform(fao_wc, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    poly_dry_wc <- spTransform(dry_wc, CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>% gBuffer(byid = T, width = 0)
    
    
    #hawaii
    haw <- readOGR(archs[3]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    haw_ml_wc <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_wc, haw)) + #ML WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "Murphy and Lugo"), 
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_ml_wc, paste("maps/poly/a_wc_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    haw_fao_wc <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_fao_wc, haw)) + #FAO WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "FAO"),
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_fao_wc, paste("maps/poly/a_wc_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    haw_dry_wc <- tm_shape(haw) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, haw)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_wc, haw)) + #Dryflor WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, haw)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = paste(abc[i], "Dryflor"),
                main.title.position = 0,
                inner.margins = c(0, 0, .2, 0),
                title = "Hawai'i",
                title.position = c("LEFT", "TOP"),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    tmap_save(haw_dry_wc, paste("maps/poly/a_wc_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #marquesas
    marq <- readOGR("data/baselayers/proc_regions/meso/mapping/poly/marquesas.shp") %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    marq_ml_wc <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_wc, marq)) + #ML WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_ml_wc, paste("maps/poly/b_wc_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    marq_fao_wc <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_fao_wc, marq)) + #FAO WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_fao_wc, paste("maps/poly/b_wc_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    marq_dry_wc <- tm_shape(marq) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(poly_countries, marq)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_wc, marq)) + #Dryflor WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, marq)) + #TDF plots
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                title = "Marquesas Islands",
                title.position = c("LEFT", "TOP"),
                fontfamily = "Times")
    
    tmap_save(marq_dry_wc, paste("maps/poly/b_wc_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #fiji
    fiji <- readOGR(archs[1]) %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    fiji_ml_wc <-  tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_ml_wc, fiji)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_ml_wc, paste("maps/poly/c_wc_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    fiji_fao_wc <- tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") + #no FAO in Fiji
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_fao_wc, paste("maps/poly/c_wc_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    fiji_dry_wc <- tm_shape(fiji) +
      tm_fill(col = "#FFFFFF") +
      tm_shape(raster::crop(poly_countries, fiji)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(poly_dry_wc, fiji)) +
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(poly_plots, fiji)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Fiji",
                title.position = c("LEFT", "TOP"),
                bg.color = "#A6CEE3")
    
    tmap_save(fiji_dry_wc, paste("maps/poly/c_wc_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    #northern mariana
    nmi <- readOGR("data/baselayers/proc_regions/meso/mapping/poly/northern_mariana.shp") %>% 
      spTransform(CRSobj = "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      gBuffer(byid = T, width = 0)
    
    
    nmi_ml_wc <-  tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") + #no ML overlap NMI
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_ml_wc, paste("maps/poly/d_wc_ml.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    nmi_fao_wc <- tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") + #no FAO overlap NMI
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_fao_wc, paste("maps/poly/d_wc_fao.png", sep = ""),
              width = 7/3, height = 2, units = "in",  asp = 0)
    
    
    nmi_dry_wc <- tm_shape(nmi) +
      tm_fill(col = "#FFFFFF") + #no dryflor overlap NMI
      tm_shape(raster::crop(poly_plots, nmi)) +
      tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_layout(title = "Northern Mariana Islands",
                title.position = c("LEFT", "TOP"),
                inner.margins = c(0, 0, 0.2, 0),
                bg.color = "#A6CEE3")
    
    tmap_save(nmi_dry_wc, paste("maps/poly/d_wc_dry.png", sep = ""),
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #ML
    pnglist = lapply(list.files("maps/poly", pattern = "*_wc_ml.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_ml <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_wc_1.png", poly_ml, width = 7/3*2, height = 7/3*2, units = "in")
    
    #FAO
    pnglist = lapply(list.files("maps/poly", pattern = "*_wc_fao.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_fao <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_wc_2.png", poly_fao, width = 7/3*2, height = 7/3*2, units = "in")
    
    #Dryflor
    pnglist = lapply(list.files("maps/poly", pattern = "*_wc_dry.png", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    poly_dry <- gridExtra::grid.arrange(grobs = groblist, ncol = 2)
    ggsave("maps/poly/poly_wc_3.png", poly_dry, width = 7/3*2, height = 7/3*2, units = "in")
    
    
  } else {
    hot <- readOGR(hotspots_tdf[i]) %>% gBuffer(byid = T, width = 0)
    
    
    hot_ml_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, hot)) + #ML WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_fao_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_wc, hot)) + #FAO WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    hot_dry_wc <- tm_shape(hot) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, hot)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, hot)) + #Dryflor WC
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, hot)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(hot_ml_wc, paste("maps/hot_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_fao_wc, paste("maps/hot_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(hot_dry_wc, paste("maps/hot_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "hot_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    hot_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig9_", i, ".png", sep = ""), 
           hot_wc, width = 7, height = 2, units = "in")
    
    
  }
} #Worldclim




#Archipelagos####
for (i in 1:4){
  arch <- readOGR(archs[i]) %>% gBuffer(byid = T, width = 0)
  
  
  if (i == 1){
    
    
    arch_ml_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, arch)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no FAO in Fiji
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, arch)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_ch, paste("maps/arch_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_ch, paste("maps/arch_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_ch, paste("maps/arch_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig10_", i, ".png", sep = ""), 
           arch_ch, width = 7, height = 2, units = "in") 
    
    
  } else if (i == 2){
    
    
    arch_ml_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, arch)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, arch)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no Dryflor in the Galapagos
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_ch, paste("maps/arch_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_ch, paste("maps/arch_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_ch, paste("maps/arch_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig10_", i, ".png", sep = ""), 
           arch_ch, width = 7, height = 2, units = "in") 
    
    
  } else {
    
    
    arch_ml_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_ch, arch)) + #ML CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_ch, arch)) + #FAO CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_ch <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_ch, arch)) + #Dryflor CHELSA
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, CHELSA") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_ch, paste("maps/arch_ch_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_ch, paste("maps/arch_ch_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_ch, paste("maps/arch_ch_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_ch_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_ch <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig10_", i, ".png", sep = ""), 
           arch_ch, width = 7, height = 2, units = "in") 
    
    
  }
  
  
} #CHELSA


for (i in 1:4){
  arch <- readOGR(archs[i]) %>% gBuffer(byid = T, width = 0)
  
  
  if (i == 1){
    
    
    arch_ml_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, arch)) + #ML Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no FAO in Fiji
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, arch)) + #Dryflor Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_wc, paste("maps/arch_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_wc, paste("maps/arch_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_wc, paste("maps/arch_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig11_", i, ".png", sep = ""), 
           arch_wc, width = 7, height = 2, units = "in")
    
    
  } else if (i == 2){
    
    
    arch_ml_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC in Galapagos
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC in Galapagos
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") + #no WC in Galapagos
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_wc, paste("maps/arch_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_wc, paste("maps/arch_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_wc, paste("maps/arch_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig11_", i, ".png", sep = ""), 
           arch_wc, width = 7, height = 2, units = "in")
    
    
  } else {
    
    
    arch_ml_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(ml_wc, arch)) + #ML Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Murphy and Lugo, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = abc[i],
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_fao_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(fao_wc, arch)) + #FAO Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "FAO, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    arch_dry_wc <- tm_shape(arch) + #basemap
      tm_fill(col = "#FFFFFF") + 
      tm_shape(raster::crop(countries, arch)) + #countries
      tm_borders(col = "#CCCCCC") +
      tm_shape(raster::crop(dry_wc, arch)) + #Dryflor Worldclim
      tm_fill(col = "#33638D") +
      tm_shape(raster::crop(plots, arch)) + #TDF plots
      tm_symbols(size = 0.025, shape = 21, col = "#CC0000",
                 border.lwd = 0.5, border.col = "#CCCCCC") +
      tm_add_legend("symbol", col = "#CC0000", 
                    shape = 16, size = 0.75,
                    labels = "Tropical Dry Forest Plots") +
      tm_add_legend("symbol", col = "#33638D",
                    shape = 15, size = 0.75, #sizing helps stack symbols
                    alpha = 0.8, labels = "Dryflor, Worldclim") +
      tm_layout(bg.color = "#A6CEE3", #ocean
                main.title = " ",
                main.title.position = 0,
                inner.margins = c(0.3, 0, 0.1, 0),
                legend.position = c(0.01, 0.05), #coords (left, bottom)
                legend.text.size = 0.5, #sizing helps stack text
                legend.width = -1, #prevents any auto re-sizing
                fontfamily = "Times")
    
    
    #save maps with insets
    tmap_save(arch_ml_wc, paste("maps/arch_wc_1.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_fao_wc, paste("maps/arch_wc_2.png", sep = ""), 
              width = 7/3, height = 2, units = "in", asp = 0)
    tmap_save(arch_dry_wc, paste("maps/arch_wc_3.png", sep = ""),  
              width = 7/3, height = 2, units = "in", asp = 0)
    
    
    #call maps back in and save
    pnglist = lapply(list.files("maps", pattern = "arch_wc_*", full.names = T), png::readPNG)
    groblist = lapply(pnglist, grid::rasterGrob)
    arch_wc <- gridExtra::grid.arrange(grobs = groblist, nrow = 1)
    ggsave(paste("maps/tdf_fig11_", i, ".png", sep = ""), 
           arch_wc, width = 7, height = 2, units = "in")
    
    
  }
  
  
} #Worldclim




#Hawai'i####
hawaii <- list.files("data/baselayers/proc_regions/micro/mapping/hawaii/",
                     pattern = ".shp", full.names = T)

#all islands + plots
hi_all <- readOGR(archs[3]) %>% gBuffer(byid = T, width = 0)


hi_plots <- tm_shape(hi_all) +
  tm_fill(col = "#FFFFFF") + 
  tm_shape(raster::crop(countries, hi_all)) + 
  tm_borders(col = "#CCCCCC") +
  tm_shape(raster::crop(ecoregions, hi_all)) + 
  tm_fill(col = "#758000", 
          alpha = 0.8) +
  tm_shape(raster::crop(plots, hi_all)) + 
  tm_symbols(size = 0.05, shape = 21, col = "#CC0000",
             border.lwd = 0.5, border.col = "#CCCCCC") +
  tm_add_legend("symbol", col = "#CC0000", 
                shape = 16, size = 0.75,
                labels = "Tropical Dry Forest Plots") +
  tm_add_legend("symbol", col = "#758000",
                shape = 15, size = 0.75, 
                alpha = 0.8, labels = "WWF TSBF Ecoregions") +
  tm_layout(bg.color = "#A6CEE3", 
            title = "Hawai'i",
            title.position = c("LEFT", "TOP"),
            inner.margins = c(0.2, 0.01, 0.2, 0.01),
            legend.position = c(0.025, 0.05), 
            legend.text.size = 0.5, 
            legend.width = -1, 
            fontfamily = "Times")

tmap_save(hi_plots, paste("maps/hawaii/hawaii.png", sep = ""),  
          width = 7/4, height = 2, units = "in", asp = 0)


#consensus maps
hi_AI <- tm_shape(hi_all) +
  tm_fill(col = "#FFFFFF") + 
  tm_shape(raster::crop(countries, hi_all)) + 
  tm_borders(col = "#CCCCCC") +
  tm_shape(raster::crop(raster(climate[8]), hi_all)) +
  tm_raster(palette = c("#3DBC74", "#33638D"),
            legend.show = F) +
  tm_add_legend("symbol", col = c("#3DBC74", "#33638D"),
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus")) +
  tm_layout(bg.color = "#A6CEE3", 
            title = "AI Consensus",
            title.position = c("LEFT", "TOP"),
            inner.margins = c(0.2, 0.01, 0.2, 0.01),
            legend.position = c(0.025, 0.05), 
            legend.text.size = 0.5, 
            legend.width = -1, 
            fontfamily = "Times")

tmap_save(hi_AI, paste("maps/hawaii/hawaii_AI.png", sep = ""),  
          width = 7/4, height = 2, units = "in", asp = 0)


hi_CH <- tm_shape(hi_all) +
  tm_fill(col = "#FFFFFF") + 
  tm_shape(raster::crop(countries, hi_all)) + 
  tm_borders(col = "#CCCCCC") +
  tm_shape(raster::crop(raster(climate[9]), hi_all)) +
  tm_raster(palette = c("#33638D", "#3DBC74", "#CAE11F"),
            legend.show = F) +
  tm_add_legend("symbol", col = c("#33638D", "#3DBC74", "#CAE11F"),
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus Among 2", "Consensus Among All")) +
  tm_layout(bg.color = "#A6CEE3", 
            title = "CHELSA Consensus",
            title.position = c("LEFT", "TOP"),
            inner.margins = c(0.2, 0.01, 0.2, 0.01),
            legend.position = c(0.025, 0.05), 
            legend.text.size = 0.5, 
            legend.width = -1, 
            fontfamily = "Times")

tmap_save(hi_CH, paste("maps/hawaii/hawaii_CH.png", sep = ""),  
          width = 7/4, height = 2, units = "in", asp = 0)


hi_WC <- tm_shape(hi_all) +
  tm_fill(col = "#FFFFFF") + 
  tm_shape(raster::crop(countries, hi_all)) + 
  tm_borders(col = "#CCCCCC") +
  tm_shape(raster::crop(raster(climate[10]), hi_all)) +
  tm_raster(palette = c("#33638D", "#3DBC74", "#CAE11F"),
            legend.show = F) +
  tm_add_legend("symbol", col = c("#33638D", "#3DBC74", "#CAE11F"),
                shape = 15, size = 0.75, #sizing helps stack symbols
                alpha = 0.8, labels = c("No Consensus", "Consensus Among 2", "Consensus Among All")) +
  tm_layout(bg.color = "#A6CEE3", 
            title = "Worldclim Consensus",
            title.position = c("LEFT", "TOP"),
            inner.margins = c(0.2, 0.01, 0.2, 0.01),
            legend.position = c(0.025, 0.05), 
            legend.text.size = 0.5, 
            legend.width = -1, 
            fontfamily = "Times")

tmap_save(hi_WC, paste("maps/hawaii/hawaii_WC.png", sep = ""),  
          width = 7/4, height = 2, units = "in", asp = 0)
