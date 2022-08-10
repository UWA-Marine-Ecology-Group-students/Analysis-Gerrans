###
# Project: Molly Thesis
# Data:    HMSC Outputs
# Task:    Predict species distributions and other outputs
# author:  Kingsley Griffin, Molly Gerrans
# date:    July-Aug 2022
##

library(Hmsc)
library(ggplot2)

# load model output
mod.dir  <- "output/hmsc_model_data"
nChains  <- 4
samples  <- 1000
thin     <- 100
filename <- file.path(paste(mod.dir), 
                      paste0("model_chains_", as.character(nChains),
                             "_samples_", as.character(samples),
                             "_thin_", as.character(thin), sep = ""), sep = "")
load(filename)


######## GRID ---------

grid <- read.csv("data/spatial/grid.csv")

head(grid)
colnames(grid)

xy.grid = as.matrix(cbind(grid$lat,grid$long))

XData.grid = data.frame(dominant.habitat=grid$dominant.habitat, 
                        location=grid$location, 
                        relief = grid$relief, 
                        depth = grid$depth,
                        stringsAsFactors = TRUE)

# We next use the prepareGradient function to convert the environmental and 
# spatial predictors into a format that can be used as input for the predict 
# function

Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(sample=xy.grid))

####UP TO HERE -----

# We are now ready to compute the posterior predictive distribution (takes a minute to compute it)
nParallel=4
predY = predict(m, Gradient=Gradient, expected = FALSE, nParallel=nParallel)


# Note that we used expected = TRUE to predict occurrence probabilities 
# (e.g. 0.2) instead of occurrences (0 or 1) 
# Note also that if you have very large prediction grid, you can use the 
# predictEtaMean = TRUE option to speed up the computations
# predY = predict(m, Gradient=Gradient, predictEtaMean = TRUE, expected = TRUE)


# Let's explore the prediction object.

class(predY)

# It is a list... 

length(predY)

# ...of length 200, if you fitted two chains with 100 samples from each

dim(predY[[1]])

# Each prediction is a matrix with dimensions 951 x 50, 
# as there are 951 prediction locations and 50 species.

head(predY[[1]])

# Each matrix is filled in with occurrence probabilities
# We may simply by ignoring parameter uncertainty and just looking at 
# the posterior mean prediction. 

EpredY=Reduce("+",predY)/length(predY)
dim(EpredY)

# EpredY is a 951 x 50 matrix of posterior mean occurrence probabilities
# The next step is to post-process the predictions to those community features
# that we wish to illustrate over the prediction space. With the script below,
# we derive from the predictions the occurrence probability of C. monedula 
# (species number 50),
# the species richness, and community-weighted mean traits.
# We also include data on habitat type and climatic conditions to the 
# dataframe mapData that
# includes all the information we need to visualize the predictions as maps

Cm = EpredY[,50]
S=rowSums(EpredY)
CWM = (EpredY%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
xy = grid[,1:2]
H = XData.grid$hab
C = grid$AprMay
mapData=data.frame(xy,C,S,Cm,CWM,H, stringsAsFactors=TRUE)

# We will use the ggplot function from the ggplot2 package, so let's load the data
library(ggplot2)

# We first plot variation in the habitat and climatic conditions on which 
# the predictions are based on.

ggplot(data = mapData, aes(x=x, y=y, color=H))+geom_point(size=2) + ggtitle("Habitat") + scale_color_discrete() + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=2) + ggtitle("Climate") + scale_color_gradient(low="blue", high="red") + coord_equal()

# We then exemplify prediction for one focal species, here C. monedula 
# over Finland

ggplot(data = mapData, aes(x=x, y=y, color=Cm))+geom_point(size=2) + ggtitle(expression(italic("Corvus monedula")))+ scale_color_gradient(low="blue", high="red") + coord_equal()

# This prediction is reassuringly very similar to that based on the 
# single-species model of Chapter 5.7
# We next plot predicted species richness, which is highest in Southern Finland

ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=2) + ggtitle("Species richness")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the proportion of resident species, also highest in Southern Finland 

ggplot(data = mapData, aes(x=x, y=y, color=MigrationR))+geom_point(size=2) + ggtitle("Proportion of resident species")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the community-weighted mean log-transformed body size, 
# which is highest in Northern Finland

ggplot(data = mapData, aes(x=x, y=y, color=LogMass))+geom_point(size=2) + ggtitle("Mean log-transformed body mass") + scale_color_gradient(low="blue", high="red") + coord_equal()
