# Take a look at the flexsdm webpage for excellent details:
# https://sjevelazco.github.io/flexsdm/index.html

#### Install.packages ----
# install.package('raster')
devtools::install_github("sjevelazco/flexsdm")
# Select 3 and then No

#### Load libraries ----
library(terra)
library(rnaturalearth)
library(flexsdm)
library(SDMtune)

#### Load in data ----
# Read species data 
protea<- vect("output/files/making_a_map/p_roup_gbif.shp")
# download and load SA boundary
sern_a <- ne_countries(scale = 'medium', country = c('South Africa', 'Lesotho', 'Swaziland'), returnclass = 'sf')
sern_a %>% group_by(level) %>% summarise() %>% vect() -> sern_a
# Download bioclim data
r <- rast(raster::getData("worldclim",var="bio",res=5))
# returns 19 variables
# or alternatively load from disk
# r <- rast("data/sdm/worldclim.tif")

#### Visualise raw data ----
plot(r[[1]])
plot(sern_a, add = T)
plot(protea, add=T)

#### Edit extent of Sern A vector 
plot(sern_a)
sern_a <- crop(sern_a, ext(15, 33, -35, -20))
plot(sern_a)

#### Check projections ----
crs(r) == crs(protea)
crs(r) == crs(sern_a)
sern_a <- project(sern_a, r)
crs(r) == crs(sern_a)

#### Calibration area ----
ca_protea <- calib_area(data = as.data.frame(protea),
                        x = 'lon', y = 'lat',
                        method = c('buffer', width = 150000),
                        crs = crs(protea))

plot(covariates[[1]])
plot(ca_protea, add = TRUE)
plot(protea, add = TRUE)

aoi <- terra::intersect(ca_protea, sern_a)
plot(covariates[[1]])
plot(aoi, add = TRUE)
plot(protea, add = TRUE)


# Mask the covariates to the area of interest
covariates <- mask(crop(r, aoi), aoi)
plot(covariates[[1]])
names(covariates)

# rename bands 
names(covariates) <- c("mean_ann_t","mean_diurnal_t_range","isothermality", "t_seas", 'max_t_warm_m','min_t_cold_m',"t_ann_range",'mean_t_wet_q','mean_t_dry_q','mean_t_warm_q','mean_t_cold_q','ann_p', 'p_wet_m','p_dry_m','p_seas','p_wet_q','p_dry_q','p_warm_q','p_cold_q')
names(covariates)

# Re-scale temperature values
covariates[[c(1:2,5:11)]] <- covariates[[c(1:2,5:11)]]/10
covariates[[3:4]] <- covariates[[3:4]]/100

#### Check for colinearity ----
# Using Pearson correlation
cov_colin <- correct_colinvar(covariates, method = c('pearson', th = "0.7"))
cov_colin$cor_table
cov_colin$cor_variables

selected_vars <- c('min_t_cold_m', 'max_t_warm_m','isothermality',
                   'ann_p','p_seas')

cov_clean <- covariates[[selected_vars]]

#### Presence filtering ----
# select only the lon/lat and provide a unique ID to each row
protea_df <- as.data.frame(protea) %>% select(lon, lat)
protea_df$id <- 1:nrow(protea_df)

occ_filt_10bin <- occfilt_env(
  data = protea_df,
  x = 'lon',
  y = 'lat',
  id = 'id',
  env_layer = cov_clean,
  nbins = 10
)

par(mfrow = c(2,1))
plot(cov_clean[[1]]); points(protea_df)
plot(cov_clean[[1]]); points(occ_filt_10bin[,2:3])

# Save the 10 bin environmental filtering
protea_filt_pres <- occ_filt_10bin[,2:3]
protea_filt_pres$pr_ab <- 1

#### Spatial block cross-validation ----
spat_range <- blockCV::spatialAutoRange(
  rasterLayer = raster::raster(cov_clean),
  speciesData = st_as_sf(protea_filt_pres, coords = c('lon','lat'), crs = crs(cov_clean)),
  doParallel = TRUE,
  showPlots = TRUE)

SB <- blockCV::spatialBlock(
  speciesData = st_as_sf(protea_filt_pres, coords = c('lon','lat'), crs = crs(cov_clean)),
  species = "pr_ab",
  rasterLayer = raster::raster(cov_clean),
  k = 4,
  iteration = 100,
  theRange = 110469
)

SB$plots
protea_filt_pres$part <- SB$foldID

# select the new dataset with partitions
protea_part <- protea_filt_pres
protea_part %>% group_by(part) %>% count()

# grid_env <- get_block(env_layer = cov_clean, best_grid = sb_cv$grid)
grid_env <- rasterize(vect(SB$blocks), cov_clean, field = 'folds')
plot(grid_env)

#### Create pseudo-absences ----
# pseudo-absences
set.seed(101)
pa <- lapply(1:4, function(x) {
  sample_pseudoabs(
  data = protea_part,
  x = 'lon',
  y = 'lat',
  n = sum(protea_part$part == x),
  method = c('env_const', env = cov_clean), # constrain to env less suitable places based on bioclim model of protea presences
  maskval = x,
  # method = 'random',
  rlayer = grid_env,
  calibarea = aoi
  )
}) %>% bind_rows() 

protea_part %>% group_by(part) %>% count()
names(pa)[4] <- 'part'
pa %>% group_by(part) %>% count()


# pa <- sdm_extract(data = pa, x = "lon", y = "lat", env_layer = grid_env)


# bg %>% group_by(.part) %>% count()

par(mfrow = c(1,1))
plot(grid_env) 
points(protea_filt_pres, cex = 0.2)
points(pa, cex = 0.2, col = 'red')
# points(bg, cex = 0.2, col = 'blue')

#### Extract covariate values ----
# Prepare a Sample with Data object for SDMtune package
SWDdata <- prepareSWD(
  species = 'Protea roupelliae',
  p = protea_filt_pres[,1:2],
  a = pa[,1:2],
  env = cov_clean
)
# For some reason ID is included, so we need to remove this...
SWDdata@data <- SWDdata@data[-1]

#### Fit Model ----
# Random folds
rand_folds <- randomFolds(SWDdata, k = 4, only_presence = TRUE, seed = 25)
rf_randcv <- train(method = 'RF', data = SWDdata, folds = rand_folds)
c('Testing AUC: ', auc(rf_randcv, test = TRUE))

# spatial blocking by specified range and random assignment
sb <- blockCV::spatialBlock(speciesData = sp_df, # sf or SpatialPoints
                   species = "SWDdata@pa", # the response column (binomial or multi-class)
                   rasterLayer = raster::raster(cov_clean)[[1]], # a raster for background (optional)
                   theRange = 110469, # size of the blocks in meters
                   k = 4, # number of folds
                   selection = "random",
                   iteration = 100) # find evenly dispersed folds

# Spatial blocking folds
rf_sbcv <- train(method = 'RF', data = SWDdata, folds = sb)
c('Testing AUC: ', auc(rf_mod_sbcv, test = TRUE))

SDMtune::confMatrix(rf_1)

train_test <- trainValTest(SWDdata, test = 0.2, seed = 25)

rf <- train(method = 'RF', data = SWDdata)

rf_1 <- rf_sbcv@models[[1]]
rf_2 <- rf_sbcv@models[[2]]

train_1 
SWDdata_1 <- SWDdata
SWDdata_1@coords$X <- SWDdata_1@coords$X[rf_sbcv@folds$train[,1]]

SWDdata_1@coords[rf_sbcv@folds$train[,1],]

SWDdata@pa

SWDdata

plotROC(rf_1, test = TRUE)



#### Variable importance ----
vi_rf_randcv <- varImp(rf_randcv)
plotVarImp(vi_rf_randcv)

vi_rf_sbcv <- varImp(rf_sbcv)
plotVarImp(vi_rf_sbcv)

#### Response curves ----
plotResponse(rf_sbcv, var = "max_t_warm_m", marginal = TRUE, rug = TRUE) + labs(x = 'Maximum temperature of warmest month')
plotResponse(rf_sbcv, var = "min_t_cold_m", marginal = TRUE, rug = TRUE) + labs(x = 'Minimum temperature of coldest month')

#### Model prediction ----
pred <- predict(rf_sbcv, data = raster::stack(cov_clean),
                type = 'prob')

plotPred(pred, lt = "Habitat\nsuitability", colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))



# # bind together the presence and pseudo-absence data
# protea_pts <- bind_rows(protea_part, pa)


# pts_extract <- sdm_extract(
#   data = protea_pts,
#   x = 'lon',
#   y = 'lat',
#   env_layer = cov_clean
# )
# head(pts_extract)


raf_mod <- fit_raf(
  data = pts_extract,
  response = 'pr_ab',
  predictors = selected_vars,
  # background = bg,
  partition = 'part',
  thr = c('max_sens_spec')
)

class(raf_mod$model)

perf <- rf_mod$performance
perf$AUC_mean
perf$TSS_mean

#### Validate 
pred_sdm <- sdm_predict(
  models = rf_mod,
  pred = cov_clean,
  thr = 'max_sens_spec'
)

plot(rast(pred_sdm))

plot(pred_sdm$raf[[1]])
plot(pred_sdm$raf[[2]])





#Assign a value of 1 to all presence points. Later, we will assign 0 to background points
df$Presence<-1
protea$ID <- df$ID


#### Background values
set.seed(42)
bg <- spatSample(covariates, length(df$Presence), "random", na.rm=TRUE, as.points=TRUE, ext=boundary)
plot(covariates$Temp)
plot(bg, add=T)
plot(protea, add=T, col = 'red')
bg$Presence<-0

pairs(covariates, cex=0.1)

#Extract covariate values at each randomly sampled background point
bgdf<- terra::extract(covariates, bg)
bg$ID <- bgdf$ID
# Assign 0 to all background points
bgdf$Presence<- 0
# combine both presence and background points and their extracted covariates
sdmdata<-data.frame(rbind(df, bgdf))
# convert the binary presence/absence points to factor variables.
sdmdata$Presence<- factor(sdmdata$Presence)
# Drop the ID column
# sdmdata$ID<-NULL
# Take a peek at the output data
glimpse(sdmdata)


ggplot(sdmdata) +
  geom_point(aes(x = Temp, y = Prec, col = Presence)) +
  ggplot(sdmdata) +
  geom_point(aes(x = Prec, y = ndvi, col = Presence)) +
  ggplot(sdmdata) +
  geom_point(aes(x = Temp, y = ndvi, col = Presence)) + plot_layout(guides = 'collect') & 
  theme_bw() & theme(legend.position = 'bottom')

# Create a train/test data partition. 
sdm_split <- initial_split(sdmdata, prop = 0.7, strata = 'Presence')

#Take a look at the training data
sdm_split %>%training() %>% glimpse()

# select the train data
traindf<- sdm_split %>% training()
# Fit a random forest model
RF <- rand_forest(trees = 100, mode = "classification") %>%
  set_engine("randomForest", importance=T) %>%
  fit(Presence ~ ., data = traindf)

predict.model_fit(RF, covariates)
# Select test data
testdf<- sdm_split %>% testing()
# Estimate accuracy
RF %>%
  predict(testdf) %>%
  bind_cols(testdf) -> RF_test

# Confusion matrix
tb <- table(RF_test$.pred_class, RF_test$Presence)
knitr::kable(tb) 
# Get a probability output
sdm_probs <- RF %>%
  predict(testdf, type = "prob") %>%
  bind_cols(testdf)

glimpse(sdm_probs)

# Plot a ROC curve
sdm_probs%>%
  roc_curve(Presence, .pred_0) %>%
  autoplot(col = .threshold) -> roc_plot
roc_plot

ggplot(roc_plot$data) +
  geom_point(aes(x = 1 - specificity, y = sensitivity, col = .threshold)) +
  # geom_path(aes(x = 1 - specificity, y = sensitivity, col = .threshold), lwd = 1) +
  scale_color_viridis_c() +
  geom_abline(aes(intercept = 0, slope = 1), lty = 3) +
  theme_bw()

ggsave('output/roc_points_w_thresh.png')
ggsave('output/roc_line_w_thresh.png')

# AUC value
sdm_probs %>% roc_auc(Presence, .pred_0)

# Obtain a dataframe with both probabilities and class output.
predict(RF, testdf, type = "prob") %>%
  bind_cols(predict(RF, testdf)) %>%
  bind_cols(select(testdf, Presence)) %>%
  glimpse()

# Run inference
r2 <- predict(covariates, RF, na.rm = T, type = 'prob')

r3 <- r2 > 0.5
plot(r3[[2]])

pred1_df <- setNames(as.data.frame(r3$.pred_1, xy = TRUE), c('x', 'y','prob_occ'))

bg_id <- bg[,c('ID','Presence')]
protea$Presence <- 1
protea_id <- protea[,c('ID','Presence')]
spat_points <- rbind(bg_id, protea_id)

points_sf <- sf::st_as_sf(spat_points)
points_sf$Presence <- factor(points_sf$Presence)
str(points_sf)

points_sf %>% left_join(testdf, by = c('ID', 'Presence')) %>% drop_na() -> test_sf


# Create a base ggplot
ggplot() +
  geom_sf(data = SA, fill = 'white', col = 'black', lwd = 0.2) +
  geom_tile(data = pred1_df, aes(x = x, y = y, fill = prob_occ), lwd = 0) +
  scale_fill_discrete(name = 'Classes', labels = c('0','1')) +
  # new_scale_fill() +
  # geom_sf(data = test_sf, aes(fill = Presence), col = 'white', pch =21) +
  # scale_fill_discrete() +
  scale_x_continuous(limits = c(22, 33)) +
  scale_y_continuous(limits = c(-34.5, -22)) +
  
  annotation_scale(location = 'br', width_hint = 0.1, style = 'ticks', tick_height = 0) + # add in a scale bar
  annotation_north_arrow(location = 'br', pad_y = unit(1, 'cm'), style = north_arrow_minimal()) + # add in a north arrow, and change around the placement and arrow style
  
  xlab('Longitude') + ylab('Latitude') +
  guides(colour = 'none') +
  coord_sf() +
  theme_bw()

ggsave('output/test_classes_protea_no_test.png')


rf_recipe <-recipe(formula=Presence ~ ., data = traindf)
rf_spec<- rand_forest(trees = 100, mode = "classification") %>%
  set_engine("randomForest", importance=T)

workflow() %>% 
  add_recipe(rf_recipe) %>% 
  add_model(rf_spec) %>% 
  fit(traindf) %>% 
  extract_fit_parsnip() %>% 
  vip(train=traindf, target='Presence',metric='accuracy',pred_wrapper = function(object, newdata) predict(object, newdata), method='permute',aesthetics = list(alpha = 0.8, fill = "midnightblue"))


##### AUC explanation
sdm_probs %>% select(Presence, .pred_1) %>% rename(Actual = Presence, `Prediction Prob.` = .pred_1) %>% 
  mutate(
`> 0.0` = case_when(
  `Prediction Prob.` > 0.0 ~ 1,
  `Prediction Prob.` <= 0.0 ~ 0), 
`> 0.1` = case_when(
  `Prediction Prob.` > 0.1 ~ 1,
  `Prediction Prob.` <= 0.1 ~ 0), 
`> 0.2` = case_when(
  `Prediction Prob.` > 0.2 ~ 1,
  `Prediction Prob.` <= 0.2 ~ 0), 
`> 0.3` = case_when(
  `Prediction Prob.` > 0.3 ~ 1,
  `Prediction Prob.` <= 0.3 ~ 0), 
`> 0.5` = case_when(
  `Prediction Prob.` > 0.5 ~ 1,
  `Prediction Prob.` <= 0.5 ~ 0), 
`> 0.7` = case_when(
  `Prediction Prob.` > 0.7 ~ 1,
  `Prediction Prob.` <= 0.7 ~ 0), 
`> 0.8` = case_when(
  `Prediction Prob.` > 0.8 ~ 1,
  `Prediction Prob.` <= 0.8 ~ 0),
`> 0.9` = case_when(
  `Prediction Prob.` > 0.9 ~ 1,
  `Prediction Prob.` <= 0.9 ~ 0),
`= 1` = case_when(
  `Prediction Prob.` == 1 ~ 1,
  `Prediction Prob.` < 1 ~ 0)
) -> auc_explained

tb0 <- table(auc_explained$Actual, auc_explained$`> 0.0`)
tb0_5 <- table(auc_explained$Actual, auc_explained$`> 0.5`)
tb1 <- table(auc_explained$Actual, auc_explained$`= 1`)

tb0
TN0 <- tb0[4]/(tb0[4] + tb0[3])
FP0 <- 1 - tb0[3]/(tb0[3] + tb0[1])

tb0_5
TN0_5 <- tb0_5[4]/(tb0_5[4] + tb0_5[3])
FP0_5 <- 1 - tb0_5[3]/(tb0_5[3] + tb0_5[1])

tb1
TN1 <- tb1[4]/(tb1[4] + tb1[3])
FP1 <- 1 - tb1[3]/(tb1[3] + tb1[1])

roc_manual <- as.data.frame(rbind(
  cbind(TN0, FP0, 0),
  cbind(TN0_5, FP0_5, 0.5),
  cbind(TN1, FP1, 1)
))

names(roc_manual) <- c('TN', 'FP', 'Threshold')

ggplot(roc_manual, aes(x = FP, y = TP)) +
  geom_point() +
  geom_path() +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))



remove.packages('cli')
remove.packages('rlang')
remove.packages('tidymodels')

install.packages('tidymodels')

.rs.restartR()
library(tidymodels)




data("spp")
spp
spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)

