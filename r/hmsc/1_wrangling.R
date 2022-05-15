###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    Wrangling data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##

library(dplyr)
library(reshape2)

#### wrangle maxn data into wide format ----

# read in maxn data
bruv_maxn <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.maxn.csv")
head(bruv_maxn)
bruv_maxn_w <- reshape2::dcast(bruv_maxn[,1:3], 
                          sample ~ scientific, fun = sum)
head(bruv_maxn_w)

write.csv(bruv_maxn_w, "data/staging/2021-05_Abrolhos_bruv_maxn_wide.csv", row.names = FALSE)


# wrangle habitat and environmental covariate info into wide format ----
colnames(bruv_maxn) # the columns of the original data that we can choose covariates from
bruv_covs <- select(bruv_maxn, c("sample", "depth", "location")) # collate all covariates we're interested in
head(bruv_covs)

write.csv(bruv_covs, "data/staging/2021-05_Abrolhos_bruv_covariates_wide.csv", row.names = FALSE)

# generate traits table including each species ----
# read in traits table

alltrait <- read.csv("data/traits/life_history.csv")
colnames(alltrait)

# reduce traits table to just our species

# get unique species from maxn table
bruv_species <- data.frame("scientific" = c(paste(bruv_maxn$genus, bruv_maxn$species, sep = " ")))
bruv_species <- unique(bruv_species)

# subset traits data to just our species
bruv_traits <- alltrait[alltrait$scientific %in% bruv_species$scientific, ]


write.csv(bruv_traits, "data/staging/2021-05_Abrolhos_bruv_traits_my_species.csv", row.names = FALSE)

# find species without traits :(
nrow(bruv_traits)

bruv_notraits <- alltraits[!alltrait$scientific %in% bruv_species$scientific, ]

nrow(bruv_notraits)



## MOLLY FIXES
# 1) species without traits = spp and Pempheris tominagi 
# (best to estimate averages for these fish based on genus information?)
#
# 2) find species without traits script


# fix!

