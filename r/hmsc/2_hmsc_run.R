###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    HMSC data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##


######## fit models -------
library(Hmsc)
set.seed(1)

# read in data
bruv_covs   <- readRDS("data/bruv_covariates_wide.rds")
bruv_maxn   <- readRDS("data/bruv_maxn_wide.rds")
bruv_traits <- readRDS("data/bruv_traits_my_species.rds")
bruv_grid <- readRDS("data/bruv_grid.rds")

# set data structure for HMSC
bruv_covs$sample   <- as.factor(bruv_covs$sample)
bruv_covs$location <- as.factor(bruv_covs$location)

# prepare our data for HMSC model structure
Y           <- bruv_maxn
XData       <- bruv_covs
TrData      <- bruv_traits
XFormula    <- ~ depth + location
TrFormula   <- ~ feeding.guild
studyDesign <- data.frame(sample = as.factor(XData$sample),
                          location = as.factor(XData$location)) # moved location to study design area (can consider moving back)
rL1 <- HmscRandomLevel(units = levels(studyDesign$sample))
rL2 <-  HmscRandomLevel(units = levels(studyDesign$location))   # added study design as nested random factor

# form the data structure required for HMSC modelling
m <- Hmsc(Y = bruv_maxn, 
          XData       = bruv_covs, 
          XFormula    = XFormula, 
          TrData      = bruv_traits,
          TrFormula   = TrFormula,
          distr       = "poisson", 
          studyDesign = studyDesign)

# cross-check what we have set up before running the model
head(m$X)
head(m$XScaled)
head(m$Tr)
head(m$TrScaled)

# we should use scaled versions of the bruv covariates - will chat about why sometime

# setup and run the actual analysis
model.directory <- "output/hmsc_model_data"

# thin 10
nChains   <- 4
nParallel <- 4
samples   <- 1000
for (thin in c(1,10)){
  transient <- 50*thin
  m <- sampleMcmc(m, thin = thin, samples = samples, 
                  transient = transient, nChains = nChains, 
                  initPar = "fixed effects", nParallel = nParallel)
  filename <- file.path(model.directory, 
                       paste0("model_chains_", as.character(nChains), 
                              "_samples_", as.character(samples), 
                              "_thin_", as.character(thin), sep = ""))
  save(m, file = filename)
}

######## Evaluating convergence ----
set.seed(1)

#new written code 
list.files("output/hmsc_model_data")


#thin 10
nChains = 4
samples = 1000
thin = 10
filename=file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

#thin 100
nChains = 4
samples = 1000
thin = 100
filename=file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)


# We restrict here the study of MCMC convergence to the examination of 
# the potential scale reduction factor of the beta parameters.

mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), 
                            covNamesNumbers = c(T,F))
str(mpost)

psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
psrf.beta


# Visual chain tests for different coefficients of interest 

# Beta estimates are the influence of environment on species
plot(mpost$Beta)

# if you review these you'll see that some species resolve well
# that is to say that the line of best fit is mostly flat and 
# the y-value range of the fuzzy line is small - the hairy caterpillar!
# eg lethrinus minatus, chrysophrys auratus
# but some others wander about (not flat) or even trend off course
# this probably comes down to how many observations we have


# use the gelman diagnostic to objectively review convergence
gelman.diag(mpost$Beta[,1:50])
# ideally values are close to 1 - ours look fairly good really

# check out the estimate values themselves
postBeta = getPostEstimate(m, parName = "Beta")
par(mar=c(5,11,2.5,0))
plotBeta(m,
         post = postBeta, 
         plotTree = F,
         spNamesNumbers = c(T,F))



plot(mpost$Gamma) #traits covariates

gelman.diag(mpost$Beta[,1:50]) #establish convergence

postBeta = getPostEstimate(m,parName = "Beta")
par(mar = c(5,11,2.5,0))
plotBeta(m, 
         post = postBeta,
         plotTree = F,
         spNamesNumbers = c(T,F))

plotBeta(m, 
         post = postBeta,
         param = "Mean",
         plotTree = F,
         spNamesNumbers = c(T,F))

postGamma = getPostEstimate(m,parName = "Gamma")
plotGamma(m, post = postGamma, supportLevel = 0.2)

VP = computeVariancePartitioning(m)
plotVariancePartitioning(m, VP = VP, las =2, horiz = F)

#gradient depth
Gradientd = constructGradient(m,focalVariable = "depth")

head(Gradientd$XDataNew)

predY = predict(m,
                XData = Gradientd$XDataNew, 
                studyDesign = Gradientd$studyDesignNew,
                ranLevels = Gradientd$rLNew,
                expected = TRUE)
predY = predict(m, Gradient=Gradientd, expected = TRUE)


plotGradient(m,
             Gradientd,
             pred=predY,
             measure="S",
             showData = TRUE)

plotGradient(m, Gradientd, pred=predY, measure="Y", index = 50, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)


#gradient location
Gradientl = constructGradient(m,focalVariable = "location")

head(Gradientl$XDataNew)

predY = predict(m,
                XData = Gradientl$XDataNew, 
                studyDesign = Gradientl$studyDesignNew,
                ranLevels = Gradientl$rLNew,
                expected = TRUE)
predY = predict(m, Gradient=Gradientl, expected = TRUE)


plotGradient(m,
             Gradientl,
             pred=predY,
             measure="S",
             showData = TRUE)

plotGradient(m, Gradientl, pred=predY, measure="Y", index = 17, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)


head(bruv_grid)

#grid = droplevels(subset(grid,!(Habitat=="Ma")))

# We next construct the objects xy.grid and XData.grid that have the 
# coordinates and environmental predictors, named similarly (hab and clim) as
# for the original data matrix (see m$XData) 

xy.grid = as.matrix(cbind(bruv_grid$latitude,bruv_grid$longitude))

#XData.grid = data.frame(hab=grid$Habitat, clim=grid$AprMay, stringsAsFactors = TRUE)

# choose XData grid points to include into grid data

# We next use the prepareGradient function to convert the environmental and 
# spatial predictors into a format that can be used as input for the predict 
# function

Gradient = prepareGradient(mpost, XDataNew = XData.grid, sDataNew = list(Route=xy.grid))

# We are now ready to compute the posterior predictive distribution (takes a minute to compute it)
nParallel=2
predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel=nParallel)

# Note that we used expected = TRUE to predict occurrence probabilities 
# (e.g. 0.2) instead of occurrences (0 or 1) 
# Note also that if you have very large prediction grid, you can use the 
# predictEtaMean = TRUE option to speed up the computations
# predY = predict(m, Gradient=Gradient, predictEtaMean = TRUE, expected = TRUE)


##### end new draft -------

# Here we take the subset of species pairs. 
# We loop over the 2 MCMC chains.

for(i in 1:length(tmp)){ 
  tmp[[i]] = tmp[[i]][,sel]
}
psrf.omega = gelman.diag(tmp,multivariate=FALSE)$psrf

par(mfrow=c(1,2))
hist(psrf.beta, xlab = "psrf (beta)")
hist(psrf.omega, xlab = "psrf (Omega)")

# The MCMC convergence diagnostics can be considered satisfactory if for
# most parameters the potential scale reduction factor is close to the 
# ideal value of one.

# For example, if most of the values are smaller than the value of 1.1. 
# If you are evaluating this script with thin=1, the largest values you 
# see are likely to be greater than ten.
# This means that MCMC convergence is not achieved, which means that the 
# posterior sample
# can be a very bad approximation of the true posterior distribution. 
# However, this should not
# stop one from making a preliminary exploration of the results, 
# just keeping in mind
# that the results can be very misleading. Thus, even if the MCMC did not 
# converge yet,
# you may proceed to the next scripts of exploring parameter estimates 
# (S3) and making predictions (S4).
# Then re-run those scripts once you have a better approximation of 
# the posterior distribution.


######## Explore model fit ------------

set.seed(1)

nChains = 2
samples = 100
thin = 1000
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)

nParallel = 2
run.cross.validation = TRUE

filename=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin)))
if(run.cross.validation){
  partition = createPartition(m, nfolds = 2, column = "Route")
  preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
  MFCV = evaluateModelFit(hM=m, predY=preds)
  save(partition,MFCV,file=filename)
} else {
  load(filename)
}

nParallel = 2
run.cross.validation = FALSE

filename=file.path(model.directory, paste0("CV_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(m$thin)))
if(run.cross.validation){
  partition = createPartition(m, nfolds = 2, column = "Route")
  preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
  MFCV = evaluateModelFit(hM=m, predY=preds)
  save(partition,MFCV,file=filename)
} else {
  load(filename)
}

tmp = c(mean(MF$AUC), mean(MFCV$AUC), mean(MF$TjurR2), mean(MFCV$TjurR2))
names(tmp)=c("AUC","AUC (CV)","TjurR2","TjurR2 (CV)")
tmp

plot(MF$AUC,MFCV$AUC)
abline(0,1)

WAIC = computeWAIC(m)
WAIC

######## explore parameter estimates -----------

set.seed(1)

nChains = 2
samples = 100
thin = 1000
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

head(m$X)

groupnames = c("habitat","climate")
group = c(1,1,1,1,1,2,2)
VP = computeVariancePartitioning(m, group = group, groupnames = groupnames)
plotVariancePartitioning(m,VP)
op = par(mar=c(9,4,1,1)+.1)
plotVariancePartitioning(m,VP, las=2)
par(op)

postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, plotTree = TRUE, spNamesNumbers=c(FALSE, FALSE))


range(m$XData$clim)

-postBeta$mean[6,]/2/postBeta$mean[7,]

postGamma = getPostEstimate(m, parName="Gamma")
plotGamma(m, post=postGamma, supportLevel = 0.95)
plotGamma(m, post=postGamma, supportLevel = 0.85)

VP$R2T$Beta
VP$R2T$Y

mpost = convertToCodaObject(m)
round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))[[2]],2)

library(corrplot)
OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)
plotOrder = corrMatOrder(OmegaCor[[1]]$mean,order="AOE")
corrplot(toPlot[plotOrder,plotOrder], method = "color", tl.cex=0.5,
         col=colorRampPalette(c("blue", "white", "red"))(255))

summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))


######## predictions ---------

set.seed(1)

nChains = 2
samples = 100
thin = 1000
filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

Gradient = constructGradient(m,focalVariable = "clim")
predY = predict(m, Gradient=Gradient, expected = TRUE)

# Occurrence probability of Corvus monedula (to see other species, change index)
plotGradient(m, Gradient, pred=predY, measure="Y", index = 50, showData = TRUE)

# Occurrence probability of Numenius arquata with a peak on the
plotGradient(m, Gradient, pred=predY, measure="Y", index = 47, showData = TRUE)

# location of the peak from the posterior mean of Beta coefficients
postBeta = getPostEstimate(m, parName="Beta")
abline(v = -postBeta$mean[6,47]/2/postBeta$mean[7,47], col="red")

# Species richness
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE)

# Proportion of species with resident migratory strategy
plotGradient(m, Gradient, pred=predY, measure="T", index = 2, showData = TRUE)

# Community mean weighted log-transformed body mass
plotGradient(m, Gradient, pred=predY, measure="T", index = 4, showData = TRUE)

Gradient = constructGradient(m,focalVariable = "hab")
predY = predict(m, Gradient=Gradient, expected = TRUE)
plotGradient(m, Gradient, pred=predY, measure="Y", index = 50, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradient, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)
 

grid = read.csv(file.path(data.directory, "grid_1000.csv"), stringsAsFactors=TRUE)
head(grid)
grid = droplevels(subset(grid,!(Habitat=="Ma")))

# We next construct the objects xy.grid and XData.grid that have the 
# coordinates and
# environmental predictors, named similarly (hab and clim) as for the original 
# data matrix (see m$XData) 

xy.grid = as.matrix(cbind(grid$x,grid$y))
XData.grid = data.frame(hab=grid$Habitat, clim=grid$AprMay, stringsAsFactors = TRUE)

# We next use the prepareGradient function to convert the environmental and 
# spatial
# predictors into a format that can be used as input for the predict function

Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(Route=xy.grid))

# We are now ready to compute the posterior predictive distribution (takes a minute to compute it)
nParallel=2
predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel=nParallel)

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
