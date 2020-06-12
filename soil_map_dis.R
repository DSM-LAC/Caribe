#load a library
library(raster)
library(RStoolbox)
library(clhs)
library(caret)
library(elevatr)
library(spatialEco)
#load a shapefile
sh <- shapefile('input/soil_map/Soils.shp')
salina <- sh[sh$DOMSOL=='Salina',]
plot(sh[!sh$DOMSOL=='Salina',], col='blue')
plot(sh[sh$DOMSOL=='Salina',], add=TRUE, col='red')
#get elevation data (9=~175m grids)
elev <- get_elev_raster(salina, prj = projection(salina), z = 10, clip = "locations")
#mask for the study area
pts <- rasterToPoints(elev, spatial=TRUE)
#load a multiple rasters with the same CRS and spatial support
predictors <- stack(list.files('input/raster_predictors/', full.name=TRUE))
#check projections
projection(predictors)==projection(sh) 
#project raster 
predictors <- projectRaster(predictors, crs=projection(sh))

train_sal <- data.frame(data.frame(extract(predictors, pts), pts@coords))

predictors_b <- mask(predictors, sh[!sh$DOMSOL=='Salina',])
#salina_predictors <- crop(predictors, salina)
#salina_predictors <- mask(salina_predictors, salina)
#apply a feature extraction strategy (PCA) to the covariates
covs_pca <- rasterPCA(predictors_b)
#converto to points the covariate file
covs_pts <- rasterToPoints(covs_pca$map, spatial=TRUE)
#define a seed for reprudicing results
set.seed(260)
#run a sampling strategy based on a latin hypercube
s.clhs <- clhs(covs_pts, size = length(pts), progress = FALSE, iter = 1000, simple = FALSE)
#get the sampled points
pts_covs <- s.clhs$sampled_data

train_no_sal <- data.frame(data.frame(extract(predictors, pts_covs), pts_covs@coords))

train_no_sal$target <- 0 
train_sal$target <- 1

reg_mat <- rbind(train_sal, train_no_sal)

reg_mat  <- na.omit(reg_mat )

#y <- as.factor(as.character(reg_mat$target))
y <- reg_mat$target
x <- reg_mat[9:12]

#train a model applied to feature extraction
set.seed(7)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=5, repeats=5)
#parallel
library(doMC)
registerDoMC(3)
names(train)
# run the RFE algorithm
results <- rfe(x=x, y=y, sizes=c(1:5), rfeControl=control)
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))
#perform a model prediction
prediction <- predict(predictors, results)
#standardize the map from 0 (no SAL) to 1 (SAL)
prediction_trans <- raster.transformation(prediction, trans = "stretch", smin = 0, smax = 1)
#plot the prediction map
plot(prediction_trans, col=colorRampPalette(c("white", "red"))(255), colNA='gray')
#add the limit of the country
lim <- getData('GADM', country='LC', level=1)
#project limit
lim <- spTransform(lim, CRS(projection(prediction)))
#plot lim
plot(lim, add=TRUE)
