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
Y <- bruv_maxn
XData <- bruv_covs
TrData <- bruv_traits
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


# thin 10
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


#thin 100
nChains = 2
nParallel = 2 # optional setting of nParallel
samples = 100
for (thin in c(1,100)) #,100,1000))
{
  transient = 50*thin
  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, initPar = "fixed effects",
                 nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(m,file=filename)
}


#thin 1000
nChains = 2
nParallel = 2 # optional setting of nParallel
samples = 100
for (thin in c(1,1000)) #,100,1000))
{
  transient = 50*thin
  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, initPar = "fixed effects",
                 nParallel = nParallel)
  filename=file.path(model.directory, paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(m,file=filename)
}

##### Evaluating convergence ----

#new writtten code 
list.files("output/hmsc/model_data")

#thin 1
nChains = 2
samples = 100
thin = 1 
filename=file.path("output/hmsc/model_data", paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

#thin 10
nChains = 2
samples = 100
thin = 10
filename=file.path("output/hmsc/model_data", paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

#thin 100
nChains = 2
samples = 100
thin = 100
filename=file.path("output/hmsc/model_data", paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

#thin 1000
nChains = 2
samples = 100
thin = 1000
filename=file.path("output/hmsc/model_data", paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We restrict here the study of MCMC convergence to the examination of 
# the potential scale reduction factor
# of the beta parameters and the Omega parameters. For the latter, we 
# take a subsample of 200 randomly selected
# species pairs to avoid excessive computations.

mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), 
                            covNamesNumbers = c(T,F))
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
tmp = mpost$Omega[[1]]
z = ncol(tmp[[1]])
sel = sample(z, size=200) #ERROR

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
