# Take a look at the flexsdm webpage for more details:
# https://sjevelazco.github.io/flexsdm/index.html

#### Install.packages ----
# install.packages('raster')
# install.packages('remotes')
# devtools::install_github("sjevelazco/flexsdm")
# Select 3 and then No
# install.packages('corrplot')
# install.packages('blockCV')
# install.packages('SDMtune')

# #### Load libraries ----
library(terra) # raster manipulation
library(rnaturalearth) # country boundaries
library(tidyverse) # data manipulation and plotting
library(corrplot) # correlation plots
library(blockCV) # spatial block for cross-validation
library(sf) # vector manipulation
library(flexsdm) # SDM package
library(SDMtune) # SDM package
library(patchwork) # combine ggplots

#### Load in data ----
# Read species data 
protea <- vect("output/files/making_a_map/p_roup_gbif.shp")
# download and load SA boundary
sern_a <- ne_countries(scale = 'medium', country = c('South Africa', 'Lesotho', 'Swaziland'), returnclass = 'sf')
# dissolve to outer boundary
sern_a %>% group_by(level) %>% summarise() %>% vect() -> sern_a_dissolve
# Download bioclim data
# r <- rast(raster::getData("worldclim",var="bio",res=5))
# returns 19 variables
# or alternatively load from disk
r <- rast("data/sdm/worldclim.tif")

#### Visualise raw data ----
plot(r[[1]])
plot(sern_a_dissolve, add = T)
plot(protea, add=T)

#### Edit extent of Sern A vector 
plot(sern_a_dissolve)
sern_a_dissolve <- crop(sern_a_dissolve, ext(15, 33, -35, -20))
plot(sern_a_dissolve)

#### Check projections ----
crs(r) == crs(protea)
crs(r) == crs(sern_a_dissolve)
sern_a_dissolve <- project(sern_a_dissolve, r)
protea <- project(protea, r)
crs(r) == crs(sern_a_dissolve)
crs(r) == crs(protea)

#### Calibration area ----
# Create a 150km buffer around our species points to create a prediction region
ca_protea <- calib_area(data = as.data.frame(protea),
                        x = 'lon', y = 'lat',
                        method = c('buffer', width = 150000),
                        crs = crs(protea))

# Plot it out to see if this works
par(mfrow = c(1,1))
plot(ca_protea)
plot(protea, add = TRUE)

# Intersect this with our country boundary to remove areas in the ocean
aoi <- terra::intersect(ca_protea, sern_a_dissolve)
plot(ca_protea)
plot(aoi, col = 'red', add = TRUE)
plot(protea, add = TRUE)

# Mask the covariates to the area of interest
covariates <- mask(crop(r, aoi), aoi)
plot(covariates[[1]])
plot(aoi, add = TRUE)

# rename the worldclim bands to more reasonable names
names(covariates)
names(covariates) <- c("mean_ann_t","mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q','mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q')
names(covariates)

# Re-scale temperature values
covariates[[c(1:2,5:11)]] <- covariates[[c(1:2,5:11)]]/10
covariates[[3:4]] <- covariates[[3:4]]/100

#### Check for colinearity ----
# Using Pearson correlation
cov_colin <- correct_colinvar(covariates, method = c('pearson', th = "0.7"))
# Take a look at the correlations using corrplot
corrplot(cov_colin$cor_table, tl.cex = 0.6)

# Show which variables are correlated
cov_colin$cor_variables

# Select the subset we want
selected_vars <- c('min_t_cold_m', 'max_t_warm_m','isothermality','ann_p','p_seas')

# Subset the covariate data
cov_clean <- covariates[[selected_vars]]
cov_clean

#### Presence filtering ----
# select only the lon/lat and provide a unique ID to each row
protea_df <- as.data.frame(protea) %>% select(lon, lat)
protea_df$id <- 1:nrow(protea_df)

# Run a filter on presence values based on the similarity in values of the environmental covariates
occ_filt_10bin <- occfilt_env(
  data = protea_df,
  x = 'lon',
  y = 'lat',
  id = 'id',
  env_layer = cov_clean,
  nbins = 10 # run ?occfilt_env to find out more on the method
)
# This removes 125 records and only keeps 99 records due to similarity!

# Plot our old points and new filtered points
par(mfrow = c(1,1))
plot(cov_clean[[1]]); 
points(protea_df, pch = 19, cex = 0.3);
points(occ_filt_10bin[,2:3], pch = 19, cex = 0.3, col = 'red')
# most of the filtering is happening where the points are spatially clustered...

# Save the 10 bin environmental filtering
protea_filt_pres <- occ_filt_10bin[,2:3]

# Assign 1 to represent presence
protea_filt_pres$pr_ab <- 1

#### Spatial block cross-validation ----
# We now want to partition our data into different spatial blocks or "folds" to avoid spatial autocorrelation
# find the range value of block sizes by fitting variograms to each environmental raster to find the effective range of spatial autocorrelation
spat_range <- spatialAutoRange(
  rasterLayer = raster::raster(cov_clean),
  speciesData = st_as_sf(protea_filt_pres, coords = c('lon','lat'), crs = crs(cov_clean)),
  doParallel = TRUE,
  showPlots = TRUE)
# This suggest ~110km as an appropriate size for our blocks

# We can now create our spatial folds (k) for later cross-validation
# set the random seed, so this the output is the same across machines
spat_blocks1 <- spatialBlock(
  speciesData = st_as_sf(protea_filt_pres, coords = c('lon','lat'), crs = crs(cov_clean)),
  species = "pr_ab",
  rasterLayer = raster::raster(cov_clean),
  k = 4,
  theRange = 110469,
  border = st_as_sf(aoi),
  seed = 101
)

# Assign the folds (or partitions) back onto the presence dataset
protea_filt_pres$folds <- spat_blocks1$foldID
# Each point is now assigned to 1 of 4 partitions
head(protea_filt_pres)

# count the number of presences in each fold
protea_filt_pres %>% group_by(folds) %>% count()

# Rasterize the blocks
grid_env <- rasterize(vect(spat_blocks1$blocks), cov_clean, field = 'folds')
plot(grid_env)

#### Create pseudo-absences ----
# We can now create our pseudo-absences and make sure they are evenly spaced across each of our fold grids, based on the number of presence points we have in each fold
# We can also make sure that they are placed away from our presence points by a environmental or distance buffer - here we use a 20 km distance buffer
# pseudo-absences
pa <- lapply(1:4, function(x) {
  sample_pseudoabs(
  data = protea_filt_pres,
  x = 'lon',
  y = 'lat',
  n = sum(protea_filt_pres$folds == x),
  # method = c('env_const', env = cov_clean), # constrain to env less suitable places based on bioclim model of protea presences
  # method = 'random',
  method = c('geo_const', width = 20000),
  maskval = x,
  rlayer = grid_env,
  calibarea = aoi
)
}) %>% bind_rows() 

# Extract the partition number for each pseudo-absence point
pa <- sdm_extract(data = pa, x = "lon", y = "lat", env_layer = grid_env)
head(pa)

# Count the number of pseudo-absence & presence points in each fold and see if they are equal
pa %>% group_by(folds) %>% count() == protea_filt_pres %>% group_by(folds) %>% count()

#### Let's plot the presences and pseudo-absences and view which folds they fall into
ggplot() +
  geom_sf(data = st_as_sf(aoi)) +
  geom_sf(data = st_as_sf(spat_blocks1$blocks)) +
  geom_point(data = rbind(protea_filt_pres, pa), aes(x = lon, y = lat, col = as.factor(folds), pch = as.factor(pr_ab))) +
  theme_void()

#### Extract covariate values ----
# Prepare a SWD (Sample with Data), which is a class of data specifically used in the SDMtune package
SWDdata <- prepareSWD(
  species = 'Protea roupelliae',
  p = protea_filt_pres[,1:2],
  a = pa[,1:2],
  env = cov_clean
)

# For some reason ID is included, so we need to remove this...
SWDdata@data <- SWDdata@data[-1]

#### Fit Model ----
# We can now fit our models! Let's first create random folds so that we can compare the output with the spatial blocking folds

###___________________________###
# RandomForest with random folds
###___________________________###
rand_folds <- randomFolds(SWDdata, k = 4, only_presence = TRUE, seed = 25)

# Run a RandomForest model with default setting and random folds
rf_randcv <- train(method = 'RF', data = SWDdata, folds = rand_folds)
# Check the overall AUC and TSS values
paste0('Testing AUC: ', round(SDMtune::auc(rf_randcv, test = TRUE),2))
paste0('Testing TSS: ', round(SDMtune::tss(rf_randcv, test = TRUE),2))

###___________________________###
# RandomForest with spatial folds
###___________________________###

# Our initial spatial blocking did not include the pseudo-absence data, so let's re-run the spatialBlock function to include both presence (1) and pseudo-absences (0)
# we can now use our previously predefined blocks
spat_blocks2 <- spatialBlock(speciesData = st_as_sf(bind_rows(protea_filt_pres, pa), coords = c('lon','lat'), crs = crs(cov_clean)), 
                   species = "pr_ab",
                   rasterLayer = raster::raster(cov_clean)[[1]], 
                   selection = 'predefined',
                   k = 4,
                   blocks = spat_blocks1$blocks,
                   foldsCol = "folds",
                   seed = 101) 
# View the output
spat_blocks2$folds
# as we can see, we have 4 folds In each of these we have our training data [[1]] and our testing data [[2]]
# also, as you can see from the console output, the training/testing split is approximately 70% training, 30% testing

# Run a RandomForest model with default setting and spatial folds
rf_sbcv <- train(method = 'RF', data = SWDdata, folds = spat_blocks2)
# Check the overall AUC and TSS values
paste0('Testing AUC: ', round(SDMtune::auc(rf_sbcv, test = TRUE),2))
paste0('Testing TSS: ', round(SDMtune::tss(rf_sbcv, test = TRUE),2))

# Extract the ROC curve and AUC values for each model
source('scripts/functions/extract_roc_vals.R')
spec_sens_vals <- extract_spec_sens_vals(rf_sbcv, spat_blocks2, SWDdata)
auc_vals <- extract_auc_vals(rf_sbcv, spat_blocks2, SWDdata)
auc_vals$label <- paste0(auc_vals$model_no, ": ", round(auc_vals$auc,2))

# ROC curves with AUC values for each model 
ggplot(data = spec_sens_vals) +
  geom_abline(aes(slope = 1, intercept = 0), lty = 2) +
  geom_path(aes(x = 1- specificities, y = sensitivities, group = model_no, col = as.factor(model_no)), alpha = 0.8) +
  scale_colour_viridis_d(name = 'Model no. & AUC',
                         labels = auc_vals$label) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.25),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

#### Variable importance ----
# VI for the random fold RF model
vi_rf_randcv <- varImp(rf_randcv)
plotVarImp(vi_rf_randcv)

# VI for the spatial fold RF model
vi_rf_sbcv <- varImp(rf_sbcv)
plotVarImp(vi_rf_sbcv)

#### Response curves ----
plotResponse(rf_sbcv, var = "ann_p", marginal = TRUE, rug = TRUE) + labs(x = 'Ann. precip.') +
plotResponse(rf_sbcv, var = "max_t_warm_m", marginal = TRUE, rug = TRUE) + labs(x = 'Max. temp. warmest month')

#### Model prediction ----
# We can now predict which areas appear most suitable to our target species across our full area of interest using  predict() and our environmental layers
pred <- predict(rf_sbcv, data = raster::stack(cov_clean))

# using the SMDtune::plotPred can give us a nice quick map
plotPred(pred, lt = "Habitat\nsuitability", colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

# We can also customise this further in ggplot if we like
pred_df <- as.data.frame(pred, xy = TRUE)

ggplot() +
  geom_sf(data = st_as_sf(sern_a), fill = 'white', col = NA) + 
  geom_tile(data = pred_df, aes(x = x, y = y, fill = layer, col = layer)) +
  scale_colour_viridis_c(na.value = NA, option = 'C', breaks = seq(0,1,0.25),limits = c(0,1)) +
  scale_fill_viridis_c(na.value = NA, option = 'C', breaks = seq(0,1,0.25),limits = c(0,1)) +
  geom_sf(data = st_as_sf(sern_a), fill = NA, col = 'black', lwd = 0.25) + 
  geom_sf(data = st_as_sf(protea_filt_pres, coords = c('lon', 'lat'), crs = crs(sern_a)), size = 0.5, col = 'black', fill = 'white', pch = 21) +
  scale_x_continuous(limits = c(ext(cov_clean)[1], ext(cov_clean)[2]), breaks = seq(26,32,3)) + 
  scale_y_continuous(limits = c(ext(cov_clean)[3],ext(cov_clean)[4]), breaks = seq(-32,-22,5)) +
  labs(fill = 'Habitat\nsuitability', 
       col = 'Habitat\nsuitability',
       x = 'Longitude', y = 'Latitude') +
  theme_minimal() 

#### END ####
