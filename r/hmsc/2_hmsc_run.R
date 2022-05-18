###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    HMSC data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##

library(Hmsc)
set.seed(1)

# read in all files
bruv_covariates <- read.csv("data/staging/2021-05_Abrolhos_bruv_covariates_wide.csv", stringsAsFactors=TRUE)
bruv_maxn <- read.csv("data/staging/2021-05_Abrolhos_bruv_maxn_wide.csv",  stringsAsFactors=TRUE)
bruv_traits <- read.csv("data/staging/2021-05_Abrolhos_bruv_traits_my_species.csv", stringsAsFactors=TRUE)

str(bruv_covariates)

nrow(bruv_covariates)

str(bruv_maxn)
str(bruv_traits)

head(bruv_traits)
head(bruv_maxn)

bruv_maxn <- bruv_maxn[ , colnames(bruv_maxn) %in% c(bruv_traits$scientific)]

#prepare files for HMSC

Y = bruv_maxn[,-1]

XData = data.frame(Sample = bruv_covariates$sample, depth=bruv_covariates$depth, location = bruv_covariates$location)

XFormula = ~ depth + poly(location,degree = 2,raw = TRUE)

TrData = data.frame(Scientific = bruv_traits$scientific, Marine_region=bruv_traits$marine.region, Feeding = bruv_traits$feeding.guild)
rownames(traits) <- traits$Scientific

TrFormula = ~Marine_region + Feeding

studyDesign = data.frame(Sample = XData$Sample)



# run HMSC

m = Hmsc(Y=Y, XData = XData, XFormula=XFormula, TrData = TrData, TrFormula = TrFormula,
         distr="probit", studyDesign = studyDesign)






