###
# Project: Molly Thesis
# Data:    HMSC Outputs
# Task:    Predict species distributions and other outputs
# author:  Kingsley Griffin, Molly Gerrans
# date:    July-Aug 2022
##

library(Hmsc)
library(ggplot2)
library(raster)

######## MAKE THE GRID FOR PREDICTION ---------
# spatial setup
wgscrs  <- CRS("+proj=longlat +datum=WGS84")
sppcrs  <- CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs")     # crs for sp objects

# load data layers we need
p_hab            <- readRDS("data/spatial/broad_habitat_predictions.rds")
p_rel            <- readRDS("data/spatial/predicted_relief_raster.rds")
p_rel[p_rel < 0] <- 0

# extract relief data and combine with other environmental covariate and predicted habitat info
plot(p_rel)
head(p_hab)

habsp <- SpatialPointsDataFrame(coords = cbind(p_hab[, 1:2]), data = p_hab)
plot(habsp, add = T)
habsp$meanrelief <- extract(p_rel, habsp)
head(habsp)

covdat  <- as.data.frame(habsp@data)
covdat  <- na.omit(covdat)
covrast <- rasterFromXYZ(covdat)
plot(covrast)

covrast$pmacroalg <- covrast$pmacroalg + covrast$pkelps                         # in your project these distributions are summed
covrast           <- covrast[[c(5, 7, 10, 12:15, 17)]]                          # get rid of all the irrelevant variables
plot(covrast)
names(covrast) <- c("slope", "detrended", "depth", "macroalgae", 
                    "unconsolidated", "consolidated", 
                    "spongegarden", "mean.relief")                              # fix the names

# start with just the northern site - this is to reduce prediction time
nsiteext <- extent(100000, 140000, 6985000, 6995000)
covrast  <- crop(covrast, nsiteext)
plot(covrast)

# fix the coordinate system so it matches our sample crs
proj4string(covrast) <- sppcrs
covrast <- projectRaster(covrast, crs = wgscrs)
covdf   <- as.data.frame(covrast, xy = TRUE, na.rm = TRUE)
head(covdf)
saveRDS(covdf, "output/covariate_grid.rds")

# save this as the grid file so you dont have to run the code above

### prepare the grid and make predictions
covdf      <- readRDS("output/covariate_grid.rds")

# load model output
mod.dir  <- "output/hmsc_model_data"
nChains  <- 4
samples  <- 1000
thin     <- 100
filename <- file.path(paste(mod.dir), 
                      paste0("model_chains_", as.character(nChains),
                             "_samples_", as.character(samples),
                             "_thin_", as.character(thin), sep = ""))
load(filename)

# setup grid and gradient
xy_grid    <- as.matrix(cbind(covdf$y, covdf$x))
XData_grid <- covdf[ , 3:10]
summary(XData_grid)
XGrad <- prepareGradient(m, 
                         XDataNew = XData_grid, 
                         sDataNew = list(sample = xy_grid))

# make predictions

# some notes on this:
# this takes a very long time and I think that's because of how I've coded it.
# potentially we can get it to perform resampling of the random effects which could speed things up.
# I think it's best to leave the re-running of this until I get back!
# also not sure how we can force it to know which location we are predicting for.

predictsp   <- predict(m, 
                       # post = poolMcmcChains(m$postList),
                       Gradient = XGrad, 
                       expected = FALSE, 
                       nParallel = 8,
                       # studydesign = studyDesign, 
                       # mcmcStep = 1000,
                       type = "response")

length(predictsp)
length(predictsp[[1]])
nrow(XData_grid)
length(predictsp[[1]])/nrow(XData_grid)
# ok - the dimensions of these estimates are really confusing me, but lets see if we can just get some predictions onto a map.

# calculate the mean prediction for each species from up to 4000 estimates - assume that's what has happened here.. why 4000 I'm not sure
# there will be a much quicker way to do this but I'm just more familiar with loops and too short on time to go hunting

for(i in 1:4000){
  predi  <- cbind(covdf[1:2], predictsp[[i]])                                   # join the coordinates and a single set of predictions
  prasti <- rasterFromXYZ(predi)                                                # make this a raster
  if(i == 1){prast <- prasti}else{                                              # if it's the first run, save it as a new object
    for(j in 1:21){                                                             # otherwise go through each layer (species)
      if(is.na(sum(prasti[[j]][]))){next}                                       # if the layer is full of na values (no prediction), skip it
      prast[[j]] <- prast[[j]] + prasti[[j]]                                    # add the new predictions to the saved predictions
    }
  }
}
prast <- prast/4000                                                             # divide by the number of estimates to get the mean
plot(prast)

# obviously the values are incorrect but looking at the relative distributions, the values look quite cool.
# the huge values could be because we used scaled (transformed) covariates in the model but didn't scale the covariates we predicted with
# I will review this fully when I get back.

# some of the species had no non-zero predictions - I assume these were quite rare?
# lets remove those for now and just look at species with full predicted distributions
rastsumm <- summary(prast)
rastsumm

predicted_sp <- dropLayer(prast, c(1,2,5,7,8,11,12,14,15, 18:21))
plot(predicted_sp)

saveRDS(predicted_sp, "output/hmsc_predictions/nthsite_predicted_spdist.rds")

## load predictions and plot

predicted_sp <- readRDS("output/hmsc_predictions/nthsite_predicted_spdist.rds")

plot(predicted_sp)

## ask claude for some tips for pretty plotting if you like or I can help when I return


