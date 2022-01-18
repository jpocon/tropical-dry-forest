library(tidyverse)
library(tm)
library(formattable)

######################################################
#### DATA PREP ####
setwd('')

## read in plot and area csv's produced from python and QGIS
cntry <- read_csv('tdf_cntry_areas.csv')
htspt <- read_csv('tdf_htspt_areas.csv')
cntry_plots <- read_csv('countries-plots.csv')
htspt_plots <- read_csv('hotspots-plots.csv')

## remove unwanted words from filenaming
stpwrds=c('tif',
          '.',
          'cntry',
          'htspt',
          'FINAL',
          'SOVEREIGNT',
          'NAME',
          'area',
          'clip')

## break apart filenaming - 'climate dataset' & 'definition' + 'region'
x <- cntry$Region %>%
  str_replace_all("-","- ") %>%
  str_replace_all("_","- ") %>%
  removeWords(stpwrds) %>%
  str_remove_all("-")

y <- htspt$Region %>%
  str_replace_all("-","- ") %>%
  str_replace_all("_","- ") %>%
  removeWords(stpwrds) %>%
  str_remove_all("-")

## add as a new col
cntry$names <- x
htspt$names <- y

## remove previous col
cntry <- cntry %>% 
  select(-Region) %>%
  group_by(names)

htspt <- htspt %>%
  select(-Region) %>%
  group_by(names)

######################################################
#### DEFINITIONS ####
## all climate extents excluding ecoregions
cntry_def <- cntry %>%
  filter_all(any_vars(str_detect(. , 'wc') | str_detect(. , 'ch')))

htspt_def <- htspt %>%
  filter_all(any_vars(str_detect(. , 'wc') | str_detect(. , 'ch')))

# Regions
defs <- c("wc" , "ch" , "65" , "ml" , "fao" , "dry")
a_no <- cntry_def$names %>%
  removeWords(defs)
cntry_def$Region <- a_no

b_no <- htspt_def$names %>%
  removeWords(defs)
htspt_def$Region <- b_no

# ML - Countries
cntry_def$Murphy_and_Lugo <- ifelse(grepl("ml ch" , cntry_def$names) , cntry_def$Area ,
                                  ifelse(grepl("ml wc" , cntry_def$names) , cntry_def$Area , NA))
# FAO
cntry_def$FAO <- ifelse(grepl("fao ch" , cntry_def$names) , cntry_def$Area ,
                                    ifelse(grepl("fao wc" , cntry_def$names) , cntry_def$Area , NA))
# Dryflor
cntry_def$Dryflor <- ifelse(grepl("dry ch" , cntry_def$names) , cntry_def$Area ,
                                    ifelse(grepl("dry wc" , cntry_def$names) , cntry_def$Area , NA))
# AI
cntry_def$Aridity <- ifelse(grepl("ch 65" , cntry_def$names) , cntry_def$Area ,
                                    ifelse(grepl("wc 65" , cntry_def$names) , cntry_def$Area , NA))
# Squish
cntry_def$Region <- str_squish(cntry_def$Region)

# ML - Hotspots
htspt_def$Murphy_and_Lugo <- ifelse(grepl("ml ch" , htspt_def$names) , htspt_def$Area ,
                                  ifelse(grepl("ml wc" , htspt_def$names) , htspt_def$Area , NA))
# FAO
htspt_def$FAO <- ifelse(grepl("fao ch" , htspt_def$names) , htspt_def$Area ,
                                    ifelse(grepl("fao wc" , htspt_def$names) , htspt_def$Area , NA))
# Dryflor
htspt_def$Dryflor <- ifelse(grepl("dry ch" , htspt_def$names) , htspt_def$Area ,
                                    ifelse(grepl("dry wc" , htspt_def$names) , htspt_def$Area , NA))
# AI
htspt_def$Aridity <- ifelse(grepl("ch 65" , htspt_def$names) , htspt_def$Area ,
                                    ifelse(grepl("wc 65" , htspt_def$names) , htspt_def$Area , NA))
# Squish
htspt_def$Region <- str_squish(htspt_def$Region)

## create final dataframe and print to csv
cntry_def <- cntry_def %>%
  select(-Area) %>%
  ungroup() %>%
  select(-names) %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarize_all(~paste(unique(na.omit(.)) , collapse = ',')) %>%
  separate_rows(Murphy_and_Lugo , FAO , Dryflor , Aridity , convert = TRUE) %>%
  write_csv('cntry-extent-km2.csv')

htspt_def <- htspt_def %>%
  select(-Area) %>%
  ungroup() %>%
  select(-names) %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarize_all(~paste(unique(na.omit(.)) , collapse = ',')) %>%
  separate_rows(Murphy_and_Lugo , FAO , Dryflor , Aridity , convert = TRUE) %>%
  write_csv('htspt-extent-km2.csv')
