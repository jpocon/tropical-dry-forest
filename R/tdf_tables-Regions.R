library(tidyverse)
library(tm)
library(formattable)

######################################################
#### DATA PREP ####
setwd('')

## read in plot and area csv's produced from python and QGIS
areas <- read_csv('tdf_areas.csv')
clim_counts <- read_csv('clim_counts.csv')
ecos_plots <- read_csv('ecos_plots.csv')

## remove unwanted words from filenaming
stpwrds=c('no' , 
          'fixed' , 
          'macro' , 
          'micro' , 
          'meso' , 
          'FINAL' , 
          'tif' , 
          '.' , 
          'countries' , 
          'islands' , 
          'ecos' , 
          'TDF' , 
          'tdf' , 
          'terr' , 
          'single' , 
          'parts' , 
          'the' , 
          'Islands' ,
          'clip')

## break apart filenaming - 'climate dataset' & 'definition' + 'region'
x <- areas$Region %>%
  str_replace_all("-","- ") %>%
  str_replace_all("_","- ") %>%
  removeWords(stpwrds) %>%
  str_remove_all("-")

## add as a new col
areas$names <- x

## remove previous col
areas <- areas %>% 
  select(-Region) %>%
  group_by(names)

######################################################
#### (TABLE 2) CLIMATE DEFS ####
## all climate extents excluding ecoregions
table_2 <- areas %>%
  filter_all(any_vars(str_detect(. , 'wc') | str_detect(. , 'ch'))) %>%
  filter_all(all_vars(!grepl('wwf' , .) & !grepl('tropics' , .)))

# Regions
defs_table_2 <- c("wc" , "ch" , "65" , "ml" , "fao" , "dry")
d_no <- table_2$names %>%
  removeWords(defs_table_2)
table_2$Region <- d_no

# ML
table_2$Murphy_and_Lugo <- ifelse(grepl("ml ch" , table_2$names) , table_2$Area ,
                                  ifelse(grepl("ml wc" , table_2$names) , table_2$Area , NA))
# FAO
table_2$FAO <- ifelse(grepl("fao ch" , table_2$names) , table_2$Area , 
                      ifelse(grepl("fao wc" , table_2$names) , table_2$Area , NA))
# Dryflor
table_2$Dryflor <- ifelse(grepl("dry ch" , table_2$names) , table_2$Area , 
                          ifelse(grepl("dry wc" , table_2$names) , table_2$Area , NA))
# AI
table_2$Aridity <- ifelse(grepl("ch 65" , table_2$names) , table_2$Area ,
                          ifelse(grepl("wc 65" , table_2$names) , table_2$Area , NA))

## relevel the order of Region col
table_2$Region <- str_squish(table_2$Region)
Regions <- c("Global" , 
             "Africa" ,
             "North and Central America" ,
             "South America" ,
             "South Asia" ,
             "South East Asia and Asia Pacific" ,
             "Caribbean" ,
             "East Melanesian" ,
             "Madagascar and Indian Ocean" ,
             "New Caledonia" ,
             "Polynesia and Micornesia" ,
             "Sundaland and Nicobar of India" ,
             "Tumbes Choco Magdalena" ,
             "Wallacea" ,
             "Fiji" ,
             "Hawaii" ,
             "Galapagos" ,
             "Puerto Rico")

## create final dataframe and print to csv
table_2 <- table_2 %>%
  select(-Area) %>%
  ungroup() %>%
  select(-names) %>%
  mutate(Region = fct_relevel(Region, Regions)) %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarize_all(~paste(unique(na.omit(.)) , collapse = ',')) %>%
  separate_rows(Murphy_and_Lugo , FAO , Dryflor , Aridity , convert = TRUE) %>%
  write_csv('table_2_extent-km2.csv')

######################################################
#### ECOREGIONS ####
## all climate extents in ecoregions
ecos <- areas %>%
  filter_all(any_vars(str_detect(. , 'ch') & str_detect(. , 'wwf') | str_detect(. , 'wc') & str_detect(. , 'wwf'))) %>%
  filter_all(all_vars(!grepl('tropics' , .)))
# Regions
defs_ecos <- c("ch" , "wc" , "65" , "ml" , "fao" , "dry" , "wwf")
d <- ecos$names %>%
  removeWords(defs_ecos)
ecos$Region <- d
# ML
ecos$Murphy_and_Lugo <- ifelse(grepl("ml ch" , ecos$names) , ecos$Area , 
                                  ifelse(grepl("ml wc" , ecos$names) , ecos$Area , NA))
# FAO
ecos$FAO <- ifelse(grepl("fao ch" , ecos$names) , ecos$Area , 
                      ifelse(grepl("fao wc" , ecos$names) , ecos$Area , NA))
# Dryflor
ecos$Dryflor <- ifelse(grepl("dry ch" , ecos$names) , ecos$Area , 
                          ifelse(grepl("dry wc" , ecos$names) , ecos$Area , NA))
# AI
ecos$Aridity <- ifelse(grepl("ch 65" , ecos$names) , ecos$Area ,
                          ifelse(grepl("wc 65" , ecos$names) , ecos$Area , NA))

## relevel the order of Region col
ecos$Region <- str_squish(ecos$Region)
Regions_ecos <- c("Global" ,
                  "Africa" ,
                  "North and Central America" ,
                  "South America" ,
                  "South Asia" ,
                  "South East Asia and Asia Pacific" ,
                  "Caribbean" ,
                  "East Melanesia" ,
                  "Madagascar and Indian Ocean" ,
                  "New Caledonia" ,
                  "Polynesia and Micronesia" ,
                  "Sundaland and Nicobar of India" ,
                  "Tumbes Choco Magdalena" ,
                  "Wallacea" ,
                  "Fiji" ,
                  "Galapagos" ,
                  "Hawaii" ,
                  "Puerto Rico")

## create final dataframe and print to csv
ecos <- ecos %>%
  select(-Area) %>%
  ungroup() %>%
  select(-names) %>%
  mutate(Region = fct_relevel(Region, Regions_ecos)) %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarize_all(~paste(na.omit(.), collapse = ',')) %>%
  separate_rows(Murphy_and_Lugo , FAO , Dryflor , Aridity , convert = TRUE) %>% 
  write_csv('ecoregion_extent-km2.csv')

######################################################
#### GLOBAL LAND AREAS ####
## no climate extents, just land areas
land <- areas %>%
  filter_all(all_vars(!grepl('ch' , .) & !grepl('wc' , .)))

## relevel the order of names col - needed to peek 'l' for spaces
land$Region <- land$names %>%
  removeWords("land") %>%
  removeWords("tropics") %>%
  str_squish()

Regions_land <- c("Global" ,
                  "Africa" , 
                  "North and Central America" ,
                  "South America" ,
                  "South Asia" ,
                  "South East Asia and Asia Pacific" ,
                  "Caribbean" , 
                  "East Melanesian" ,
                  "Madagascar and Indian Ocean" ,
                  "New Caledonia" ,
                  "Polynesia and Micornesia" ,
                  "Sundaland and Nicobar of India" ,
                  "Tumbes Choco Magdalena" ,
                  "Wallacea" ,
                  "Fiji" ,
                  "Galapagos" ,
                  "Hawaii" ,
                  "Puerto Rico" ,
                  "wwf Global" ,
                  "wwf Africa" ,
                  "wwf North and Central America" ,
                  "wwf South America" ,
                  "wwf South Asia" ,
                  "wwf South East Asia and Asia Pacific" ,
                  "wwf Caribbean" ,
                  "wwf East Melanesia" ,
                  "wwf Madagascar and Indian Ocean" ,
                  "wwf New Caledonia" ,
                  "wwf Polynesia and Micronesia" ,
                  "wwf Sundaland and Nicobar of India" ,
                  "wwf Tumbes Choco Magdalena" ,
                  "wwf Wallacea" ,
                  "wwf Fiji" ,
                  "wwf Galapagos" ,
                  "wwf Hawaii" ,
                  "wwf Puerto Rico")

## create final dataframe and print to csv
land <- land %>%
  ungroup() %>%
  mutate(Region = fct_relevel(Region, Regions_land)) %>%
  arrange(Region) %>%
  select(-names) %>%
  write_csv('total_land_areas.csv')

######################################################
#### (TABLE 3) ECOREGION AGREEMENT ####
## calculate the overlap of total ecoregion area and definition
table_3 <- ecos
table_3$total_ecos_area <- rep(land$Area[19:36], each = 2)

table_3 <- table_3 %>%
  mutate(ML_ecos_agree = (Murphy_and_Lugo / total_ecos_area * 100) ,
            FAO_ecos_agree = (FAO / total_ecos_area * 100) ,
            Dry_ecos_agree = (Dryflor / total_ecos_area * 100) ,
            AI_ecos_agree = (Aridity / total_ecos_area * 100)) %>%
  mutate_if(is.numeric, round, digits = 0) %>%
  select(-Murphy_and_Lugo) %>%
  select(-FAO) %>%
  select(-Dryflor) %>%
  select(-Aridity)

table_3[is.na(table_3)] <- 0
write_csv(table_3 , 'table_3_ecos_agree.csv')

######################################################
#### (TABLE 4) TDF PLOT AGREEMENT ####
## in order to match rows/cols, had to manually add 'Galapagos' + 'Sundaland' 
## and 'Global' (which included replacing missing areas - 
## i.e. Fiji, Galapagos, and East Melanesian) to 'areas' csv
clim_plots <- clim_counts %>% 
  group_by(Region) %>%
  select(-c(LAT , LON , 
            wc65 , ch65 , mlwc , mlch , faowc , faoch , drywc , drych , 
            COUNTRY , MACRO , MESO , MICRO , REF , DATABASE)) %>%
  gather(defs , value , ALL:DRY_CH) %>%
  unite("Regions" , Region:defs , remove = TRUE) %>%
  filter_all(all_vars(!grepl('NA_' , .))) %>%
  ungroup()

## break apart filenaming - 'climate dataset' & 'definition' + 'region'
defs_clim_plots <- c("WC" , "CH" , "65" , "ML" , "FAO" , "DRY" , "ALL")
clim_plots$Region <- clim_plots$Regions %>%
  str_replace_all("_","- ") %>%
  removeWords(defs_clim_plots) %>%
  str_remove_all("-")

# Global
clim_plots$All <- ifelse(grepl("ALL" , clim_plots$Regions) , clim_plots$value ,
                            ifelse(grepl("ALL" , clim_plots$Regions) , clim_plots$value , NA))
#ML
clim_plots$Murphy_and_Lugo <- ifelse(grepl("_ML_CH" , clim_plots$Regions) , clim_plots$value ,
                                  ifelse(grepl("_ML_WC" , clim_plots$Regions) , clim_plots$value , NA))
# FAO
clim_plots$FAO <- ifelse(grepl("_FAO_CH" , clim_plots$Regions) , clim_plots$value ,
                                     ifelse(grepl("_FAO_WC" , clim_plots$Regions) , clim_plots$value , NA))
# Dryflor
clim_plots$Dryflor <- ifelse(grepl("_DRY_CH" , clim_plots$Regions) , clim_plots$value ,
                                     ifelse(grepl("_DRY_WC" , clim_plots$Regions) , clim_plots$value , NA))
# AI
clim_plots$Aridity <- ifelse(grepl("_CH_65" , clim_plots$Regions) , clim_plots$value ,
                                     ifelse(grepl("_WC_65" , clim_plots$Regions) , clim_plots$value , NA))

## relevel the order of Region col
clim_plots$Region <- str_squish(clim_plots$Region)
Regions <- c("Global" ,              
             "Africa" ,
             "North and Central America" ,
             "South America" ,
             "South Asia" ,
             "South East Asia and Asia Pacific" ,
             "Caribbean" ,
             "East Melanesia" ,
             "Madagascar and Indian Ocean" ,
             "New Caledonia" ,
             "Polynesia and Micronesia" ,
             "Sundaland and Nicobar" ,
             "Tumbes Choco Magdalena" ,
             "Wallacea" ,
             "Fiji" ,
             "Galapagos" ,
             "Hawaii" ,
             "Puerto Rico")

clim_plots <- clim_plots %>%
  ungroup() %>%
  select(-Regions) %>%
  select(-value) %>%
  mutate(Region = fct_relevel(Region, Regions)) %>%
  arrange(Region) %>%
  group_by(Region) %>%
  summarize_all(~paste(na.omit(.) , collapse = ',')) %>%
  separate_rows(All , Murphy_and_Lugo , FAO , Dryflor , Aridity , convert = TRUE) %>%
  write_csv('clim_plots.csv')

## create final dataframe and print to csv
table_4 <- table_3

table_4$Clim_plots <- rep(clim_plots$All[1:36])
table_4$Ecos_plots <- rep(ecos_plots$count[1:18] , each = 2)
table_4$ML_plots <- rep(clim_plots$Murphy_and_Lugo[1:36])
table_4$FAO_plots <- rep(clim_plots$FAO[1:36])
table_4$Dry_plots <- rep(clim_plots$Dryflor[1:36])
table_4$AI_plots <- rep(clim_plots$Aridity[1:36])

table_4 <- table_4 %>% 
  mutate(Eco_plot_agree = (Ecos_plots / Clim_plots * 100) ,
         ML_plot_agree = (ML_plots / Clim_plots * 100) ,
         FAO_plot_agree = (FAO_plots / Clim_plots * 100) ,
         Dry_plot_agree = (Dry_plots / Clim_plots * 100) ,
         AI_plot_agree = (AI_plots / Clim_plots * 100)) %>%
  mutate_if(is.numeric, round, digits = 0) %>%
  select(-total_ecos_area) %>%
  select(-ML_ecos_agree) %>%
  select(-FAO_ecos_agree) %>%
  select(-Dry_ecos_agree) %>%
  select(-AI_ecos_agree) %>%
  select(-Ecos_plots) %>%
  select(-ML_plots) %>%
  select(-FAO_plots) %>%
  select(-Dry_plots) %>%
  select(-AI_plots)

table_4[is.na(table_4)] <- 0
write_csv(table_4 , 'table_4_plots_agree.csv')

######################################################
#### (TABLES 5-7) BIODIVERSITY HOTSPOTS ####
## HOTSPOT COUNTS
bio_2 <- table_2 %>%
  slice(13:28) %>%
  write_csv('bio_2_extent-km2.csv')

## HOTSPOT ECOS AGREEMENT
bio_3 <- table_3 %>%
  slice(13:28) %>%
  write_csv('bio_3_ecos_agree.csv')

## HOTSPOT PLOT AGREEMENT
bio_4 <- table_4 %>%
  slice(13:28) %>%
  write_csv('bio_4_plots_agree.csv')

######################################################
#### (TABLE 2) FIGURES ####
## create new dframes for each region and climate set
## Global
tab2_figs_global <- table_2[1:2,]
tab2_figs_global$Data <- "NA"
tab2_figs_global$Data[1] <- "CHELSA"
tab2_figs_global$Data[2] <- "Worldclim"

## index the rows for each region and climate set
idx_sub_ch <- seq(1,9,2)
idx_sub_wc <- seq(2,10,2)
## Subcontinents
tab2_figs_sub <- table_2[3:12,]
tab2_figs_sub$Data <- "NA"
tab2_figs_sub$Data[idx_sub_ch] <- "CHELSA"
tab2_figs_sub$Data[idx_sub_wc] <- "Worldclim"

idx_bio_ch <- seq(1,15,2)
idx_bio_wc <- seq(2,16,2)
## Biodiversity Hotspots
tab2_figs_bio <- table_2[13:28,]
tab2_figs_bio$Data <- "NA"
tab2_figs_bio$Data[idx_bio_ch] <- "CHELSA"
tab2_figs_bio$Data[idx_bio_wc] <- "Worldclim"

idx_arch_ch <- seq(1,7,2)
idx_arch_wc <- seq(2,8,2)
## Archipelagos
tab2_figs_arch <- table_2[29:36,]
tab2_figs_arch$Data <- "NA"
tab2_figs_arch$Data[idx_arch_ch] <- "CHELSA"
tab2_figs_arch$Data[idx_arch_wc] <- "Worldclim"

## create a total areas fig and a percent areas fig 
## Global
Fig_Global <- tab2_figs_global %>%
  gather(Defs , value = Region , Murphy_and_Lugo:Aridity) %>%
  ggplot(aes(x = Defs , y = Region , fill = Data)) +
  geom_bar(stat="identity", width = .5 , position = position_dodge()) +
  labs(title = "Figure 13", 
       subtitle = "Global Climatic Definition Extents (km2)", 
       caption = "Data: CHELSA (Karger et al., 2017), Worldclim (Fick and Hijmans, 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#E69F00" , "#56B4E9"))
Fig_Global

## Subcontinents
Fig_Sub_CH_km <- tab2_figs_sub[idx_sub_ch,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 14", 
       subtitle = "Subcontinent Climatic Definition Extents (km2) CHELSA", 
       caption = "Data: CHELSA (Karger et al., 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8)) +
  scale_x_discrete(labels=c("Africa" = "Africa" ,
                            "North and Central America" = "N+C America" ,
                            "South America" = "S America" ,
                            "South Asia" = "S Asia" ,
                            "South East Asia and Asia Pacific" = "SE Asia + Pacific"))
Fig_Sub_CH_km

Fig_Sub_WC_km <- tab2_figs_sub[idx_sub_wc,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 15", 
       subtitle = "Subcontinent Climatic Definition Extents (km2) Worldclim", 
       caption = "Data: Worldclim (Fick and Hijmans, 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8)) +
  scale_x_discrete(labels=c("Africa" = "Africa" ,
                            "North and Central America" = "N+C America" ,
                            "South America" = "S America" ,
                            "South Asia" = "S Asia" ,
                            "South East Asia and Asia Pacific" = "SE Asia + Pacific"))
Fig_Sub_WC_km

## Biodiversity Hotspots
Fig_Bio_CH_km <- tab2_figs_bio[idx_bio_ch,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 16", 
       subtitle = "Biodiversity Hotspot Climatic Definition Extents (km2) CHELSA", 
       caption = "Data: CHELSA (Karger et al., 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8)) +
  scale_x_discrete(labels=c("Caribbean" = "Caribbean" ,
                            "East Melanesian" = "E Melanesia" ,
                            "Madagascar and Indian Ocean" = "Madagascar" ,
                            "New Caledonia" = "New Caledonia" ,
                            "Polynesia and Micornesia" = "Polynesia" ,
                            "Sundaland and Nicobar of India" = "Sundaland" , 
                            "Tumbes Choco Magdalena" = "Tumbes" ,
                            "Wallacea" = "Wallacea"))
Fig_Bio_CH_km

Fig_Bio_WC_km <- tab2_figs_bio[idx_bio_wc,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 17", 
       subtitle = "Biodiversity Hotspot Climatic Definition Extents (km2) Worldclim", 
       caption = "Data: Worldclim (Fick and Hijmans, 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8)) +
  scale_x_discrete(labels=c("Caribbean" = "Caribbean" ,
                            "East Melanesian" = "E Melanesia" ,
                            "Madagascar and Indian Ocean" = "Madagascar" ,
                            "New Caledonia" = "New Caledonia" ,
                            "Polynesia and Micornesia" = "Polynesia" ,
                            "Sundaland and Nicobar of India" = "Sundaland" , 
                            "Tumbes Choco Magdalena" = "Tumbes" ,
                            "Wallacea" = "Wallacea"))
Fig_Bio_WC_km

## Archipelagos
Fig_Arch_CH_km <- tab2_figs_arch[idx_arch_ch,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 18", 
       subtitle = "Archipelagos Climatic Definition Extents (km2) CHELSA", 
       caption = "Data: CHELSA (Karger et al., 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8))
Fig_Arch_CH_km

Fig_Arch_WC_km <- tab2_figs_arch[idx_arch_wc,] %>%
  rename(ML = Murphy_and_Lugo) %>%
  gather(Definitions , value , ML:Aridity) %>%
  ggplot(aes(x = Region , y = value , fill = Definitions)) +
  geom_bar(stat = "identity" , width = .5 , position = position_dodge()) + 
  labs(title = "Figure 19", 
       subtitle = "Archipelagos Climatic Definition Extents (km2) Worldclim", 
       caption = "Data: Worldclim (Fick and Hijmans, 2017)" , 
       x = '' , 
       y = 'Area (Km2)') + 
  scale_fill_manual(values = c("#999999" , "#56B4E9" , "#E69F00" , "#15D180")) +
  theme(axis.text.x = element_text(angle = 20 , hjust = 0.8))
Fig_Arch_WC_km