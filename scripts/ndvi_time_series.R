#### Install packages ----
install.packages('MODISTools')
# Find the MODISTools vignette here:
#https://cran.r-project.org/web/packages/MODISTools/vignettes/modistools-vignette.html

#### Load packages ----
library(MODISTools)
library(tidyverse)
library(terra)
library(lubridate)

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
jhb_ndvi <- mt_subset(product = "VNP13A1",
                          lat = -26.2041,
                          lon =  28.0473,
                          band = c("500_m_16_days_NDVI",
                                   "500_m_16_days_pixel_reliability"),
                          start = "2021-01-01",
                          end = "2021-12-30",
                          km_lr = 20,
                          km_ab = 20,
                          site_name = "JHB",
                          internal = TRUE,
                          progress = TRUE)

#### Plot time series
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
  theme_minimal()

#### Convert to raster
jhb_ndvi %>% filter(band == "500_m_16_days_NDVI" & calendar_date == '2021-01-01') -> jhb_ndvi_01_01
jhb_ndvi_01_01_rast <- rast(mt_to_raster(jhb_ndvi_01_01, reproject = TRUE))
plot(jhb_ndvi_01_01_rast)

jhb_ndvi %>% filter(band == "500_m_16_days_NDVI" & calendar_date == '2021-01-25') -> jhb_ndvi_01_25
jhb_ndvi_01_25_rast <- rast(mt_to_raster(jhb_ndvi_01_25, reproject = TRUE))
plot(jhb_ndvi_01_25_rast)

# Check quality indicator
# range of indicator values
jhb_ndvi %>% filter(band == "500_m_16_days_pixel_reliability") %>% summarise(min(value), max(value))

jhb_ndvi %>% 
  filter(band == "500_m_16_days_pixel_reliability" & calendar_date == '2021-06-18') -> jhb_QI_01_01
table(jhb_QI_01_01$value)

jhb_QI_01_01_rast <- rast(mt_to_raster(jhb_QI_01_01, reproject = TRUE))
plot(jhb_QI_01_01_rast)

jhb_ndvi %>% 
  filter(band == "500_m_16_days_pixel_reliability" & calendar_date == '2021-01-25') -> jhb_QI_01_25
table(jhb_QI_01_25$value)

jhb_QI_01_25_rast <- rast(mt_to_raster(jhb_QI_01_25, reproject = TRUE))
plot(jhb_QI_01_25_rast)

#### Plot RGB ----
jhb_rgb <- mt_subset(product = "VNP13A1",
                      lat = -26.2041,
                      lon =  28.0473,
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

#### Download NDVI for multiple points (batch download)
lat <- c(-26.2041, 5.6037, -1.2921)
lon <- c(28.0473, -0.1870, 36.8219)
site_name <- c('JHB', 'ACC', 'NAI')
coords <- data.frame(site_name, lat, lon)

cities_ndvi <- mt_batch_subset(df = coords,
                      product = "VNP13A1",
                      band = c("500_m_16_days_NDVI"),
                      start = "2021-01-01",
                      end = "2021-12-30",
                      km_lr = 2,
                      km_ab = 2,
                      internal = TRUE)

cities_ndvi %>%
  group_by(site, calendar_date) %>%
  summarise(doy = yday(as_date(calendar_date)),
            ndvi_median = median(value * as.numeric(scale))) %>%
  distinct(doy, .keep_all = TRUE) %>%
  ggplot(aes(x = doy, y = ndvi_median)) +
    geom_point() +
    geom_smooth(method = 'loess') +
    labs(x = 'Day of Year', y = 'NDVI (median)') +
    theme_classic() +
    facet_grid(cols = vars(site))

