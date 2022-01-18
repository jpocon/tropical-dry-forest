# load raster packages
require(raster)
require(rgdal)

# define working directory
cwd = ''
setwd(cwd)

# create a list of tifs -- CHELSA -- stack, sum, plot, write file to cwd
list_ch <- list.files(cwd , pattern = "*ch-clip.tif$")
stack_ch <- raster::stack(list_ch)
sum_ch <- sum(stack_ch , na.rm = TRUE)
plot(sum_ch)
writeRaster(sum_ch , 'tdf_overlap_ch' , format = 'GTiff')

# repeat for Worldclim
list_wc <- list.files(cwd , pattern = "*wc-clip.tif$")
stack_wc <- raster::stack(list_wc)
sum_wc <- sum(stack_wc , na.rm = TRUE)
plot(sum_wc)
writeRaster(sum_wc , 'tdf_overlap_wc', format = 'GTiff')

# repeat for Aridity Index
list_ai <- list.files(cwd , pattern = "AI")
stack_ai <- raster::stack(list_ai)
sum_ai <- sum(stack_ai , na.rm = TRUE)
plot(sum_ai)
writeRaster(sum_ai , 'tdf_overlap_ai', format = 'GTiff')