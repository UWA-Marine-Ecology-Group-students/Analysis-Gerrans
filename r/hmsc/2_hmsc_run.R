###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    HMSC data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##

library(Hmsc)
set.seed(1)

# read in data
bruv_covs   <- readRDS("data/bruv_covariates_wide.rds")
bruv_maxn   <- readRDS("data/bruv_maxn_wide.rds")
bruv_traits <- readRDS("data/bruv_traits_my_species.rds")

# set data structure for HMSC
bruv_covs$sample   <- as.factor(bruv_covs$sample)
bruv_covs$location <- as.factor(bruv_covs$location)

# prepare our data for HMSC model structure
XFormula    <- ~ depth + location
TrFormula   <- ~ feeding.guild
studyDesign <- data.frame(sample = as.factor(XData$sample))

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
model.directory <- "output/hmsc/model_data"

nChains = 2
nParallel = 2 # optional setting of nParallel
samples = 100
for (thin in c(1,10)) #,100,1000))
{
  transient = 50*thin
  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, initPar = "fixed effects",
                 nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(m,file=filename)
}



