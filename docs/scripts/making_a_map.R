#### Install and load packages ----
# install.packages('rgbif')
# install.packages('tidyverse')
# install.packages('CoordinateCleaner')
# install.packages('rnaturalearth')
# install.packages('sf')
# install.packages('mapview')
# install.packages('ggspatial')
# install.packages('patchwork')

library(tidyverse)
library(rgbif)
library(CoordinateCleaner)
library(rnaturalearth)
library(sf)
library(mapview)
library(ggspatial)
library(patchwork)

# # Alternative method to load many libraries
# if(!require("pacman")) install.packages("pacman")
# pacman::p_load(tidyverse, rgbif, CoordinateCleaner, rnaturalearth, sf, mapview, ggspatial, patchwork)

#### Download species data from GBIF ----

# Name your species
myspecies <- c("protea roupelliae")

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_download <- occ_data(scientificName = myspecies, hasCoordinate = TRUE, limit = 1000)

gbif_data <- gbif_download$data
names(gbif_data)

#### Clean GBIF data ----

gbif_data %>%
  filter(year >= 1900) %>%
  filter(!coordinateUncertaintyInMeters %in% c(301,3036,999,9999)) %>% # known inaccurate default values
  filter(!decimalLatitude == 0 | !decimalLongitude == 0) %>%
  cc_cen(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(lat = 'decimalLatitude', lon = 'decimalLongitude', buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea(lat = 'decimalLatitude', lon = 'decimalLongitude') %>% # remove from ocean 
  distinct(decimalLongitude,decimalLatitude,speciesKey,datasetKey, .keep_all = TRUE) -> p_roup

#### Download map layers ----
# using the rnaturalearth package, download and convert the download to an sf object
world <- ne_countries(scale = 'medium', returnclass = 'sf')
names(world)

# We can easily plot the world map layer using ggplot
ggplot() +
  geom_sf(data = world, aes(fill = gdp_md_est)) + # use the fill option within the aesthetics parameter
  scale_fill_viridis_c(trans = 'log', name = 'GDP (log)') # chose a palette for the colour gradient, transform the variable and add a name for the legend

# We can also change the projection of the coordinate reference system (CRS)
ggplot() +
  geom_sf(data = world, aes(fill = gdp_md_est)) +
  scale_fill_viridis_c(trans = 'log', name = 'GDP (log)') +
  coord_sf(crs = 3035)

# Let's move our focus back to our species map for Protea roupelliae
# Create a vector of names for countries we want to use (= southern African countries) to filter to our area of focus
countries <- c('South Africa', 'Lesotho', 'Swaziland', 'Namibia', 'Botswana', 'Zimbabwe', 'Mozambique')

# Now filter down the list of all countries to only include southern African countries
world %>% filter(admin %in% countries) -> sern_a

# Plot this using geom_sf from ggplot()
ggplot() +
  geom_sf(data = sern_a) 

# This looks a bit strange, it includes Marion Island in the sub-antarctic, which is part of South Africa!
# Remove this by setting limits on the x and y axes
base_plot1 <- ggplot() + # save the output of the map to an object
  geom_sf(data = sern_a) +
  scale_x_continuous(limits = c(16, 33)) +
  scale_y_continuous(limits = c(-35, -22))
base_plot1 # print this object, we can now reuse 'base_plot1' and add to it later

# Now we're focused in on South Africa, let's add the GBIF points to the map
# First clean the dataset by selecting only the columns we need, rename the variables
p_roup %>% select(key, year, basisOfRecord, decimalLongitude, decimalLatitude) %>% rename(lon = decimalLongitude, lat = decimalLatitude) -> pr_clean

# Now we want to convert & project our data to the correct Coordinate Reference System
pr_no_crs <- st_as_sf(pr_clean, coords = c('lon', 'lat'))
st_crs(pr_no_crs) # as we converted these directly from GPS points, there is no set CRS

# Let's try plot this
base_plot1 +
  geom_sf(data = pr_no_crs)

# Error! This doesn't work, so we need to give the points a coordinate reference system (CRS)
st_crs(sern_a) # first see what CRS the map layers have
pr_WGS_84 <- pr_no_crs # create a new object for the points
st_crs(pr_WGS_84) = 4326 # assign it the CRS 4326, which is short for the World Geodetic System 1984 (WGS 84), a spatial reference system. 
st_crs(pr_WGS_84) # check that this worked

# Let's try plot this now
base_plot1 +
  geom_sf(data = pr_WGS_84, size = 0.8) 

# Some of these points look a bit strange. It is always good to inspect these to see if they may be an issue later. We can make a simple interactive map to explore this...
mapview(pr_WGS_84)

# you should notice when you zoom in that many points are too regularly spaced to be exact locations of species sightings; rather, such points are likely to be centroids of (relatively large) grid cells on which particular surveys were based. This is not necessarily a problem, but remember to adjust the spatial resolution of your analysis accordingly!

# let's separate the countries with points from those without
sa_les_esw <- sern_a %>% filter(admin %in% c('South Africa', 'Lesotho', 'Swaziland'))
other_sern_a <- sern_a %>% filter(!admin %in% c('South Africa', 'Lesotho', 'Swaziland'))

# Filter out points pre 1980
pr_WGS_84 %>% filter(year > 1980) -> pr_1980

#### Let's prepare our final plot ----
base_plot2 <- ggplot() +
  geom_sf(data = other_sern_a, fill = 'gray90', lwd = 0.3) + # change the fill colour
  geom_sf(data = sa_les_esw, fill = 'white', lwd = 0.3) + # change the fill colour
  geom_sf(data = pr_1980, alpha = 0.6, size = 0.8) + # change the transparency and size of the points  
  scale_x_continuous(limits = c(16, 33), breaks = seq(16, 32, 4)) + # add in breaks for the labels using sequence function
  scale_y_continuous(limits = c(-35, -22), breaks = seq(-34, -22, 4)) + # add in breaks for the labels using sequence function
  xlab('Longitude') + ylab('Latitude') + # change the x and y axes labels
  annotation_scale(location = 'br', width_hint = 0.1, style = 'ticks', tick_height = 0) + # add in a scale bar
  annotation_north_arrow(location = 'br', pad_y = unit(1, 'cm'), style = north_arrow_minimal()) + # add in a north arrow, and change around the placement and arrow style
  annotate(geom = 'text', x = 23.5, y = -30.5, label = 'SOUTH AFRICA', colour = 'black', size = 4) + # add in text for SA
  annotate(geom = 'text', x = 28.25, y = -29.5, label = 'LESOTHO', colour = 'black', size = 2.6) + # add in text for Lesotho
  annotate(geom = 'text', x = 31.45, y = -26.6, label = 'ESWATINI', colour = 'black', size = 2.6) + # add in text for eSwatini
  theme_bw() + # edit the overall theme (there are several other options)
  theme(panel.grid = element_line(color = gray(.5), linetype = 'dashed', size = 0.3),
        panel.background = element_rect(fill = '#E8F6FC')) # lastly, change the grid line parameters and change the background (ocean) colour

base_plot2

#### Save the map ----
# first as a pdf, specify the output directory, the object and the width & height
ggsave('output/figs/making_a_map/map.pdf', base_plot2, width = 6.68, height = 5.76)
# next as a png, specify the same parameters, but now as it is a pixel-based file type, specify the resolution
ggsave('output/figs/making_a_map/map.png', base_plot2, width = 6.68, height = 5.76, dpi = 300)

# This is a reasonable place to end, however, if your map is of a small area (luckily not our case), it can often be useful to add an inset map to show where your site/study species is relative to other known features

#### Inset map ----
# Let's create a box that covers our map to show where our site is.
# We will use the limits of the previous map 

# A rectangle requires 5 points to make a complete polygon, so we will create a list with 2 columns (lon & lat) and 5 rows (5 points)
pol_lon = c(16, 33)
pol_lat = c(-35, -22)
coord_df = data.frame(pol_lon, pol_lat)

# Now make the list
coord_list <- list(
  cbind(
    coord_df$pol_lon[c(1,2,2,1,1)], 
    coord_df$pol_lat[c(1,1,2,2,1)])
)

# Convert the points to an sf polygon
base_plot2_bbox <- st_polygon(coord_list, dim = 'XY')
# Add a CRS that matches our previous systems
map_bbox <- st_sfc(base_plot2_bbox, crs = 4326)

# Plot to see if this is working
ggplot() +
  geom_sf(data = world) +
  geom_sf(data = map_bbox, fill = NA, col = 'red') # use fill = NA to make the middle of the polygon transparent

#### Let's now zoom in our inset map and try to match some of the other features of our basemap
# Download just Africa
africa <- ne_countries(continent = 'africa', returnclass = 'sf')
other_africa <- africa %>% filter(!admin %in% c('South Africa', 'Lesotho', 'Swaziland'))

# Create the inset map
inset_map <- ggplot() +
  geom_sf(data = other_africa, fill = 'gray90', lwd = 0.2) +
  geom_sf(data = sa_les_esw, fill = 'white', lwd = 0.2) +
  geom_sf(data = map_bbox, fill = NA, col = 'black', lwd = 0.8) +
  scale_x_continuous(limits = c(-20, 50)) +
  scale_y_continuous(limits = c(40, -35)) +
  theme_void() + # this theme removes ALL axes, grids and labels 
  theme(panel.border = element_rect(colour = 'black', fill = NA),
        panel.background = element_rect(fill = '#E8F6FC')) # we then want to add a border and the background colour back in

# We then use the patchwork package to add the inset map to our basemap
combined_maps <- base_plot2 + inset_element(inset_map, left = 0.025, bottom = 0.6, right = 0.4, top = 0.975) 
# we need to chose the positioning of the inset, assuming the plot values go from 0 to 1 for each axis.
combined_maps

# Save the output
ggsave('output/figs/making_a_map/combined_maps.png', combined_maps, width = 6.68, height = 5.76, dpi = 300)  

#### Export the data files for later ----
# write the P. roupelliae locations to a shapefile
write_sf(pr_1980, 'output/files/making_a_map/p_roup_gbif.shp')

#### END ####