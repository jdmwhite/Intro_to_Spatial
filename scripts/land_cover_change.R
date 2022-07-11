# install.packages('ggalluvial') 

library(sf)
library(terra)
library(tidyverse)
library(ggalluvial)
library(patchwork)

#Read 2020 Land cover (raster)
lc2020 <- rast("data/LULC/SANLC_2020_COJ_extent.tif")
#Read the COJ boundary. we can omit "terra::" because this function name does not overlap with any other R package function.
coj <- vect("data/LULC/COJ_boundary.shp")

#
print(paste0("The CRS of both layers are the same: ", crs(coj)==crs(lc2020)))

coj <- project(coj, crs(lc2020))
print(paste0("The CRS of both layers are the same: ", crs(coj)==crs(lc2020)))

# Crop
lc2020_cropped = crop(lc2020, coj)
# Mask
lc2020_final = mask(lc2020_cropped, coj)
#Convert raster to Vector
lc2020_poly = as.polygons(lc2020_final, dissolve =T) %>% st_as_sf()
#Rename the column with the class name to "class"
colnames(lc2020_poly)<- c("class", "geometry")

#Load LC1990
lc1990 = rast("data/LULC/SANLC_1990_COJ_extent.tif")
#Reproject to match lc2020
lc1990<- terra::project(lc1990,lc2020)
# Crop
lc1990_cropped = crop(lc1990, coj)
# Mask
lc1990_final = mask(lc1990_cropped, coj)
#Convert raster to Vector (sf)
lc1990_poly = as.polygons(lc1990_final, dissolve =T) %>% st_as_sf()
#Rename the column with the class name to "class"
colnames(lc1990_poly)<- c("class", "geometry")
#sanity check
print(unique(lc2020_poly$class))


# reclassify 1990 and 2020 land cover

# Key
# 1= Water
# 2= Agriculture
# 3= Artificial surfaces
# 4 = Vegetation

#1990 reclasssify

lc1990_rcl<-dplyr::mutate(lc1990_poly, 
                          rcls_1990 = case_when(
                            (class %in% c(1, 2, 3, 37, 38)) ~ 1,
                            (class %in% c(10:31)) ~ 2,
                            (class %in% c(35, 36, 39:51, 53:56, 61:72)) ~ 3,
                            (class %in% c(4:9, 32:34, 52, 57:60)) ~ 4,
                            TRUE~0))
#2020 reclasssify
lc2020_rcl<-mutate(lc2020_poly, 
                   rcls_2020 = case_when(
                     (class %in% c(14:24)) ~ 1,
                     (class %in% c(32:46)) ~ 2,
                     (class %in% c(25:31, 47:60, 65:73)) ~ 3,
                     (class %in% c(1:13, 14:19, 20:25, 61:64)) ~ 4,
                     TRUE~0))

lc2020r  <- rasterize(vect(lc2020_rcl), lc2020_final, 'rcls_2020')

lc1990r  <- rasterize(vect(lc1990_rcl), lc1990_final, 'rcls_1990')

#stack land cover
landcover_stack <- c(lc2020r, lc1990r)
#Change analysis. long =T returns a data frame instead of a table.
Changes<-crosstab(landcover_stack, long=T)

Changes %>% mutate(
  area = Freq*900/1e6,
  lc1990 = case_when(
    rcls_1990 == 1 ~ 'Water',
    rcls_1990 == 2 ~ 'Agriculture',
    rcls_1990 == 3 ~ 'Artificial',
    rcls_1990 == 4 ~ 'Vegetation'),
  lc2020 = case_when(
    rcls_2020 == 1 ~ 'Water',
    rcls_2020 == 2 ~ 'Agriculture',
    rcls_2020 == 3 ~ 'Artificial',
    rcls_2020 == 4 ~ 'Vegetation')
) %>% select(lc1990, lc2020, area) -> Changes
####

alluv_plot <- ggplot(Changes, aes(axis1 = lc1990, axis2 = lc2020, y = area)) +
  geom_alluvium(aes(fill = lc1990)) +
  scale_fill_manual(values = c('#7E6148B2','#F39B7FB2','#00A087B2','#4DBBD5B2'), guide = 'none') +
  geom_stratum(fill = c('#4DBBD5B2','#00A087B2','#F39B7FB2','#7E6148B2','#4DBBD5B2','#00A087B2','#F39B7FB2','#7E6148B2'), col = NA, alpha = 0.8) +
  geom_text(stat = 'stratum', aes(label = paste(after_stat(stratum),'\n',round(after_stat(prop)*100,1))), size = 2.5) +
  scale_x_continuous(breaks = c(1, 2), labels = c('1990','2020'), position = 'top') +
  theme_void() +
  theme(axis.text.x = element_text())

alluv_plot

#### Plot maps
# upsample the rasters using aggregate
lc1990ra <- aggregate(lc1990r, 5, fun = 'median')
lc1990_df <- as.data.frame(lc1990ra, xy = TRUE)
names(lc1990_df)[3] <- 'land_cover'

lc2020ra <- aggregate(lc2020r, 5, fun = 'median')
lc2020_df <- as.data.frame(lc2020ra, xy = TRUE)
names(lc2020_df)[3] <- 'land_cover'

lc1990_plot <- ggplot(lc1990_df) +
  geom_tile(aes(x = x, y = y, fill = as.factor(land_cover))) +
  scale_fill_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'),
                    labels = c('Water', 'Agriculture', 'Artificial', 'Vegetation'), guide = 'none') +
  geom_sf(data = st_as_sf(coj), fill = NA, col = 'black', lwd = 0.2) +
  labs(title = '1990', fill = 'Land Cover') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) 
  
lc2020_plot <- ggplot(lc2020_df) +
  geom_tile(aes(x = x, y = y, fill = as.factor(land_cover))) +
  scale_fill_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'),
                    labels = c('Water', 'Agriculture', 'Artificial', 'Vegetation'), guide = 'none') +
  geom_sf(data = st_as_sf(coj), fill = NA, col = 'black', lwd = 0.2) +
  labs(title = '2020', fill = 'Land Cover') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) 


lc_plots <- lc1990_plot + lc2020_plot + alluv_plot & plot_annotation(tag_levels = 'a', tag_suffix = ')')

ggsave('output/figs/land_cover_plots.png', lc_plots,
       width = 180, height = 100, units = c('mm'), dpi = 'retina')
