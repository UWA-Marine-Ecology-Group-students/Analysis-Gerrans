###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    HMSC data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##


######## Fit models -------
library(Hmsc)
set.seed(1)

# read in data
bruv_covs   <- readRDS("data/bruv_covariates_wide.rds")
bruv_maxn   <- readRDS("data/bruv_maxn_wide.rds")
bruv_traits <- readRDS("data/bruv_traits_my_species.rds")
bruv_xy     <- readRDS("data/bruv_xy.rds")

# set data structure for HMSC
bruv_covs$sample   <- as.factor(bruv_covs$sample)
bruv_covs$location <- as.factor(bruv_covs$location)
str(bruv_covs)

# prepare our data for HMSC model structure
Y           <- bruv_maxn
XData       <- bruv_covs
TrData      <- bruv_traits
XFormula    <- ~ depth + macroalgae + mean.relief + spongegarden + consolidated + unconsolidated + slope + detrended
TrFormula   <- ~ splength + trophic_group + complexity + substrate_group + water_column + night_day
studyDesign <- data.frame(sample = as.factor(XData$sample), location = as.factor(XData$location))
rL          <- HmscRandomLevel(sData = bruv_xy)
#rL1 <- HmscRandomLevel(units = levels(studyDesign$sample))
#rL2 <-  HmscRandomLevel(units = levels(studyDesign$location))   # added study design as nested random factor


# form the data structure required for HMSC modelling
m <- Hmsc(Y = bruv_maxn, 
          XData       = bruv_covs, 
          XFormula    = XFormula, 
          TrData      = bruv_traits,
          TrFormula   = TrFormula,
          distr       = "poisson", 
          studyDesign = studyDesign,
          ranLevels   = list(sample = rL))

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
for (thin in c(100)){
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
thin = 100
filename=file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# #thin 100
# nChains = 4
# samples = 1000
# thin = 100
# filename=file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
# load(filename)


# We restrict here the study of MCMC convergence to the examination of 
# the potential scale reduction factor of the beta parameters.

mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), 
                            covNamesNumbers = c(T,F))
str(mpost)

psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
psrf.beta


######## Explore model fit --------

# Visual chain tests for different coefficients of interest 

# Beta estimates are the influence of environment on species
# plot(mpost$Beta)

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

plotBeta(m, 
         post = postBeta,
         param = "Mean",
         plotTree = F,
         spNamesNumbers = c(T,F))


plot(mpost$Gamma) #traits covariates

gelman.diag(mpost$Beta[,1:50]) #establish convergence


postGamma = getPostEstimate(m,parName = "Gamma")
plotGamma(m, post = postGamma, supportLevel = 0.2)

# variance partitioning
VP = computeVariancePartitioning(m)
dev.off()
plotVariancePartitioning(m, VP = VP, las =2, horiz = F)
