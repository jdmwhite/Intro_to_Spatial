# netCDF downloaded from:
# http://aphrodite.st.hirosaki-u.ac.jp/products.html

#### Install packages ----
# install.packages('animation')

#### Load packages ----
library(terra)
library(animation)

#### Load netCDF file ----
ma_temp_2015 <- rast('data/netCDF_example/APHRO_MA_TAVE_050deg_V1808.2015.nc')
ma_temp_2015
names(ma_temp_2015)

# warning provided for no extent, so we will need to assign this manually from the metadata
ext(ma_temp_2015)
ma_temp_2015 <- set.ext(ma_temp_2015, c(60, 150, -15, 55))
ext(ma_temp_2015)

# Plot a single raster
plot(ma_temp_2015[[1]])

# This looks upside down! 
# We can use flip to turn this the correct way around
ma_temp_2015_flip <- flip(ma_temp_2015)

# Plot the flipped raster
plot(ma_temp_2015_flip[[1]])

# Your netCDF file is now ready to analyse!
saveGIF(animate(ma_temp_2015_flip[[1:50]]),
        interval = 0.1,
        ani.height = 200, ani.width = 300,
        movie.name = 'netCDF_example.gif')

