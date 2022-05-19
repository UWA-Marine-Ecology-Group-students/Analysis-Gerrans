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
colnames(bruv_maxn)

# remake scientific name with just genus and species for matching up with other dfs
bruv_maxn$genus_species <- c(paste(bruv_maxn$genus, 
                                   bruv_maxn$species, sep = "_"))
head(bruv_maxn$genus_species)
bruv_maxn_w  <- reshape2::dcast(bruv_maxn[, c(1, 3, 20)], 
                                sample ~ genus_species, 
                                fun = sum, value.var = "maxn")
rownames(bruv_maxn_w) <- bruv_maxn_w$sample                                     # make sample ids row names
rownames(bruv_maxn_w) <- gsub("\\.", "_", bruv_maxn_w$sample)
bruv_maxn_w <- bruv_maxn_w[, -1]                                                # drop sample column
head(bruv_maxn_w)

# wrangle habitat and environmental covariate info into wide format ----
bruv_meta <- read.csv("data/raw/em export/2021-05_Abrolhos_stereo-BRUVs_Metadata.csv.csv")
colnames(bruv_meta)                                                             # the columns of the original data that we can choose covariates from
bruv_covs <- select(bruv_meta, c("Sample", "Latitude", "Longitude", 
                                 "Depth", "Location"))                          # collate all covariates we're interested in
head(bruv_covs)

# collapse rows to make it one row per sample (match with bruv_maxn)
bruv_covs <- unique(bruv_covs)
nrow(bruv_covs)
nrow(bruv_maxn_w)
# there are now more maxn rows than there are metadata - maybe we didn't complete a few drops?

# for now remove those extra drop records from the covariates
bruv_covs$Sample <- gsub("\\.", "_", bruv_covs$Sample)                          # the periods were causing issues later i think
bruv_covs <- bruv_covs[bruv_covs$Sample %in% rownames(bruv_maxn_w), ]

# finally fix column names
colnames(bruv_covs) <- tolower(colnames(bruv_covs))

head(bruv_covs)
nrow(bruv_covs) == nrow(bruv_maxn_w)                                            # do row numbers match in maxn and covariates?

# generate traits table including each species ----
# read in traits table
alltrait <- read.csv("data/traits/life_history.csv")
colnames(alltrait)

# reduce traits table to just our species
bruv_species        <- colnames(bruv_maxn_w)                                    # species from maxn column names
alltrait$scientific <- gsub(" ", "_", alltrait$scientific)                      # adding . to species names for consistency
bruv_traits         <- alltrait[alltrait$scientific %in% bruv_species, ]

# how many species are wemissing traits for?
length(bruv_species) - nrow(bruv_traits)

# reduce traits data to just our covariates of interest - use ::summary to choose traits with most data
summary(bruv_traits)
interesting_traits <- c("scientific", "feeding.guild", "fb.vulnerability")
bruv_traits <- bruv_traits[ , colnames(bruv_traits) %in% interesting_traits]

# clean up/remove species from traits data if they have NA in any trait info
bruv_traits <- na.omit(bruv_traits)
head(bruv_traits)

# clean traits data itself
bruv_traits$feeding.guild <- tolower(bruv_traits$feeding.guild)                 # make all lowercase
bruv_traits$feeding.guild <- sapply(bruv_traits$feeding.guild, 
               FUN = function(x)paste(unique(unlist(strsplit(x, "/"))), 
                                      collapse = " "))                          # remove punctuation and duplicates
unique(bruv_traits$feeding.guild)

# make species names row names
rownames(bruv_traits) <- bruv_traits$scientific
bruv_traits <- bruv_traits[ , -1]

# drop bruv maxn info for species that are lacking any traits
bruv_maxn_w <- bruv_maxn_w[ , colnames(bruv_maxn_w) %in% c(rownames(bruv_traits))]
dim(bruv_maxn_w)

# write to RDS to preserve row names
saveRDS(bruv_maxn_w, "data/bruv_maxn_wide.rds")
saveRDS(bruv_covs,   "data/bruv_covariates_wide.rds")
saveRDS(bruv_traits, "data/bruv_traits_my_species.rds")


# list species without traits :(
bruv_notrait <- bruv_species[(bruv_species %in% bruv_traits$scientific) == FALSE]
length(bruv_notrait)
bruv_notrait

## MOLLY FIXES
# 1) species without traits = spp and Pempheris tominagi 
# (best to estimate averages for these fish based on genus information?)
# 2) find species without traits script


# fix!

# clear environment
rm(list=ls())
