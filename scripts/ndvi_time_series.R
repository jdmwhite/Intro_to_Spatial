#### Install packages ----
# install.packages('MODISTools')
# install.packages('lubridate')
# install.packages('animation')

# Find the MODISTools vignette here:
#https://cran.r-project.org/web/packages/MODISTools/vignettes/modistools-vignette.html

#### Load packages ----
library(MODISTools)
library(tidyverse)
library(terra)
library(lubridate)
library(animation)

#### Load in additional data----
coj <- vect("data/LULC/COJ_boundary.shp")

#### Explore modis products ----
# All products available
products <- mt_products()
head(products)

# bands of the vegetation indices product
bands <- mt_bands(product = "VNP13A1")
head(bands)

# dates available for a specific coordinate
dates <- mt_dates(product = "VNP13A1", lat = -26.2041, lon = 28.0473)
head(dates)

#### Download NDVI data for JHB ----
start_time <- Sys.time()
jhb_ndvi <- mt_subset(product = "VNP13A1",
                          lat = -26.183483,
                          lon =  27.975311,
                          band = c("500_m_16_days_NDVI",
                                   "500_m_16_days_pixel_reliability"),
                          start = "2021-01-01",
                          end = "2021-12-30",
                          km_lr = 25,
                          km_ab = 40,
                          site_name = "JHB",
                          internal = TRUE,
                          progress = TRUE)
end_time <- Sys.time()
# How long does it take to download 200 km2 of MODIS NDVI data at 500 m resolution for 2 bands?
end_time - start_time

#### Plot time series ----
# Summarise
jhb_ndvi %>% 
  filter(band == "500_m_16_days_NDVI") %>%
  group_by(calendar_date) %>%
  summarise(doy = yday(as_date(calendar_date)),
            ndvi_median = median(value * as.numeric(scale))) %>% 
  distinct(doy, .keep_all = TRUE) -> jhb_med_ndvi

# Plot
ggplot(jhb_med_ndvi, aes(x = doy, y = ndvi_median)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(x = 'Day of Year', y = 'NDVI (median)') +
  theme_classic() 
  
# Save plot
ggsave('output/figs/ndvi_time_series/jhb_med_ndvi.png',
       width = 120, height = 80, units = c('mm'), dpi = 'retina', bg = 'white')

#### Convert NDVI to raster ----
# filter to NDVI and then split the data frame into a list by date
jhb_ndvi_split <- jhb_ndvi %>% filter(band == "500_m_16_days_NDVI") %>% split(jhb_ndvi$calendar_date)

# convert each date into a raster (and reproject)
jhb_ndvi_rast_list <- lapply(jhb_ndvi_split, function(x) {rast(mt_to_raster(x, reproject = TRUE))})

# convert the list of rasters into a raster stack
jhb_ndvi_rasts <- rast(jhb_ndvi_rast_list)
jhb_ndvi_rasts <- mask(jhb_ndvi_rasts, coj)

# plot out two dates
plot(jhb_ndvi_rasts$`2021-01-01`)
plot(coj, add = T)
plot(jhb_ndvi_rasts$`2021-01-25`)
plot(coj, add = T)

#### Convert Pixel Reliability to raster ----
# filter to NDVI and then split the data frame into a list by date
jhb_pr_split <- jhb_ndvi %>% filter(band == "500_m_16_days_pixel_reliability") %>% split(jhb_ndvi$calendar_date)

# convert each date into a raster (and reproject)
jhb_pr_rast_list <- lapply(jhb_pr_split, function(x) {rast(mt_to_raster(x, reproject = TRUE))})

# convert the list of rasters into a raster stack
jhb_pr_rasts <- rast(jhb_pr_rast_list)
jhb_pr_rasts <- mask(jhb_pr_rasts, coj)

# plot out a few dates
plot(jhb_pr_rasts$`2021-01-01`)
plot(jhb_pr_rasts$`2021-01-25`)
# these are unreliable images
plot(jhb_pr_rasts$`2021-06-18`) # this is a reliable image
# scores below 3 are acceptable

#### Plot RGB ----
jhb_rgb <- mt_subset(product = "VNP13A1",
                     lat = -26.183483,
                     lon =  27.975311,
                     band = c("500_m_16_days_red_reflectance",
                               "500_m_16_days_green_reflectance",
                               '500_m_16_days_blue_reflectance'),
                      start = "2021-01-25",
                      end = "2021-01-25",
                      km_lr = 20,
                      km_ab = 20,
                      site_name = "JHB",
                      internal = TRUE,
                      progress = TRUE)

jhb_rgb_split <- split(jhb_rgb, jhb_rgb$band)

jhb_rgb_rasts <- lapply(jhb_rgb_split, function(x) {rast(mt_to_raster(x, reproject = TRUE))})

jhb_rgb_rast <- c(jhb_rgb_rasts[[1]], jhb_rgb_rasts[[2]], jhb_rgb_rasts[[3]])

plotRGB(jhb_rgb_rast, stretch = 'hist')
plot(coj, add = TRUE)

#### Download NDVI for multiple points ----
# Batch download
# Create a dataframe of sites and lat/lon or load in a csv
lat <- c(-26.2041, 5.6037, -1.2921)
lon <- c(28.0473, -0.1870, 36.8219)
site_name <- c('JHB', 'ACC', 'NAI')
coords <- data.frame(site_name, lat, lon)

# Use the mt_batch_subset function to load in data for multiple points
cities_ndvi <- mt_batch_subset(df = coords,
                      product = "VNP13A1",
                      band = c("500_m_16_days_NDVI"),
                      start = "2021-01-01",
                      end = "2021-12-30",
                      km_lr = 2,
                      km_ab = 2,
                      internal = TRUE)

# Summarise the values by site and date and pipe this directly into a ggplot call
# Use facets to display a smoothed NDVI curve for each site
cities_ndvi %>%
  group_by(site, calendar_date) %>%
  summarise(doy = yday(as_date(calendar_date)),
            ndvi_median = median(value * as.numeric(scale))) %>%
  distinct(doy, .keep_all = TRUE) %>%
  ggplot(aes(x = doy, y = ndvi_median)) +
    geom_point(size = 0.8, alpha = 0.8) +
    geom_smooth(method = 'loess', col = 'darkgreen') +
    labs(x = 'Day of Year', y = 'NDVI (median)') +
    theme_classic() +
    facet_grid(cols = vars(site))

# Save plot
ggsave('output/figs/ndvi_time_series/cities_med_ndvi.png',
       width = 180, height = 80, units = c('mm'), dpi = 'retina', bg = 'white')

# Animate ----
# Bonus - animate your NDVI rasters
saveGIF(animate(jhb_ndvi_rasts),
        interval = 0.1,
        ani.height = 200, ani.width = 300,
        movie.name = 'jhb_ndvi_animation.gif')

