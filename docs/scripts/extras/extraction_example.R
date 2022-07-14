# Goal: to extract environmental data using points or polygons

#### Libraries
library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)

##____________________________##
#### Load in environmental data ----
##____________________________##

# Load in worldclim data
worldclim <- rast('data/sdm/worldclim.tif')
# Subset the data to the first 6 variables just to make the example run faster
wc <- worldclim[[1:6]]
# check the variables
names(wc)

##____________________________##
#### Extract using polygon data ----
##____________________________##

# load in a country borders for Africa from online
# you could load in your own shapefile here
afr <- vect(ne_countries(continent = 'Africa'))
head(afr)
# there are lots of variables, let's keep only the country name
afr <- afr[,'name']
plot(afr)

# Check projections and reproject if needed
crs(afr) == crs(wc)
afr_prj <- project(afr, wc)
crs(afr_prj) == crs(wc)

# Crop & mask the worldclim layer to the Africa layer
wc_mask <- mask(crop(wc, afr_prj), afr_prj)
plot(wc_mask)

# Create a custom function to summarise the environmental raster data for each country
my_summary <- function(x) c(mean = mean(x, na.rm = T), min = min(x, na.rm=T), max = max(x, na.rm=T))

# Use terra::extract() to get multiple summary values
poly_ext <- terra::extract(wc_mask, afr_prj, fun = my_summary)
head(poly_ext)
# Bind the result back on to the polygons
afr_ext <- cbind(afr_prj, poly_ext)

# Convert the data to an sf object for plotting
afr_ext_sf <- afr_ext %>% st_as_sf()

# Plot
ggplot() +
  geom_sf(data = afr_ext_sf, aes(fill = bio1.mean/10)) +
  scale_fill_gradientn(colours = c('white', 'red'),
                       name = 'Mean ann. temp.') +
  theme_void()

##____________________________##
#### Extract using point data ----
##____________________________##

# Create 1000 random points over Africa
rand_pts <- terra::spatSample(x = afr, size = 1000, method = "random")
plot(afr)
points(rand_pts, cex = 0.5, col = 'red')

# Check projection
crs(rand_pts) == crs(wc)
rand_pts <- project(rand_pts, wc)
crs(rand_pts) == crs(wc)

# Write points to file and then read back in
writeVector(rand_pts, 'data/extraction_example/random_points.shp')

# Save as a data frame in csv format
xy <- geom(rand_pts, df = TRUE)[,c(1,3:4)]
write_csv(xy, 'data/extraction_example/random_points.csv')

# NOTE: if you have your own point data, this is how you can read it into R
# Read points back in
# 1: as a shapefile
rand_pts <- vect('data/extraction_example/random_points.shp')
head(rand_pts)
# 2: as a csv
rand_pts_df <- read_csv('data/extraction_example/random_points.csv')
rand_pts <- vect(rand_pts_df, geom = c('x','y'), crs = crs(wc))
head(rand_pts)

# So that's a few different ways to create or load in your point data in different formats, let's now run the extraction

# Extract data for points
# We don't need to specify a function, because these are just points and will extract data for one cell that they intersect with
pts_ext <- terra::extract(wc, rand_pts)
head(pts_ext)
rand_pts_ext <- cbind(rand_pts, pts_ext)

# Convert to sf for plotting
rand_pts_ext_sf <- rand_pts_ext %>% st_as_sf()
afr_sf <- afr %>% st_as_sf()

# Plot
ggplot() +
  geom_sf(data = afr_sf, fill = NA) +
  geom_sf(data = rand_pts_ext_sf, aes(fill = bio1/10), pch = 21, size = 1) +
  scale_fill_gradientn(colours = c('white', 'red'),
                       name = 'Mean ann. temp.') +
  theme_void()

##____________________________##
#### Buffer around points ----
##____________________________##

# Buffers are essentially polygons, so this will work in a very similar way to our first example.
buf_pts <- buffer(rand_pts, width = 10000) # 10km buffer
plot(afr)
plot(buf_pts, add = T)

# Extract data for buffers, we can use the same my_summary function as earlier because these are also polygons with multiple raster tiles intersecting with them. We therefore need to calculate summary values. 
buf_ext <- terra::extract(wc, buf_pts, fun = my_summary)
# warnings() produced due to some points falling in the ocean I presume...
buf_pts_ext <- cbind(buf_pts, buf_ext)

# Convert to sf for plotting
buf_pts_ext_sf <- buf_pts_ext %>% st_as_sf()

# Plot
ggplot() +
  geom_sf(data = afr_sf, fill = NA) +
  geom_sf(data = buf_pts_ext_sf, aes(col = bio1.mean/10), pch = 1, size = 1) +
  scale_colour_gradientn(colours = c('white', 'red'),
                       name = 'Mean ann. temp.') +
  theme_void()

#### END ####