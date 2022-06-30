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
rownames(bruv_maxn_w) <- bruv_maxn_w$sample                         # make sample ids row names
rownames(bruv_maxn_w) <- gsub("\\.", "_", bruv_maxn_w$sample)
bruv_maxn_w <- bruv_maxn_w[, -1]                                    # drop sample column
head(bruv_maxn_w)


#### wrangle habitat and environmental covariate info into wide format ----

bruv_meta    <- read.csv("data/raw/em export/2021-05_Abrolhos_stereo-BRUVs_Metadata.csv.csv")
bruv_habitat <- read.csv("data/tidy/2021-05_Abrolhos_BRUVs_random-points_percent-cover_broad.habitat.csv")
colnames(bruv_meta)                                                             # the columns of the original data that we can choose covariates from
bruv_covs <- select(bruv_meta, c("Sample", "Latitude", "Longitude", 
                                 "Depth", "Location"))                          # collate all covariates we're interested in
bruv_covs <- unique(bruv_covs)                                                  # collapse rows to make it one row per sample (match with bruv_maxn)
head(bruv_covs)

# clean habitat data
head(bruv_habitat)
colnames(bruv_habitat)    <- gsub("broad.", "", colnames(bruv_habitat))         # shorten column names
bruv_habitat$spongegarden <- rowSums(bruv_habitat[, colnames(bruv_habitat) %in% 
                                                     c("ascidians", "bryozoa", 
                                                       "hydroids", "sponges",
                                                       "invertebrate.complex", 
                                                       "octocoral.black")])     # collapse invertebrate reef tags

# add habitat to BRUV covars
head(bruv_habitat)
bruv_covs <- merge(bruv_covs, bruv_habitat[, colnames(bruv_habitat) %in%
                                             c("sample", "consolidated",
                                               "macroalgae", "unconsolidated",
                                               "spongegarden", "mean.relief")], 
                   by.x = "Sample", by.y = "sample")
bruv_covs[, 6:10] <- round(bruv_covs[, 6:10], 2)
head(bruv_covs)

# check dimensions against drops
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


#### generate traits table including each species ----

# read body length and mass
# bodylength <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.length.csv") # both measures appear in mass table
bodymass   <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.mass.csv")

# clean up/remove species from length/mass data if they have NA in any measurement info
# bodylength <- na.omit(bodylength)
bodymass   <- na.omit(bodymass)

# calculate mean body measures per fish
bodydim <- summarise(group_by(bodymass, family, genus, species),
                     meanlength = sum(length)/sum(number),
                     meanmass = sum(mass.g)/sum(number))
bodydim$scientific  <- interaction(bodydim$genus, bodydim$species, sep = "_")
head(bodydim)

# check correlation among length and mass
plot(log(bodydim$meanlength), log(bodydim$meanmass))
cor(bodydim$meanlength, bodydim$meanmass) 

# read in full traits table, clean column names and shorten from entire sheet to just our species
alltrait   <- read.csv("data/traits/Australia.life.history - australia.life.history.csv")
colnames(alltrait) <- gsub("\\.", "_", colnames(alltrait))
colnames(alltrait) <- tolower(colnames(alltrait))
colnames(alltrait)

bruv_species        <- colnames(bruv_maxn_w)                                    # species from maxn column names
head(alltrait)
alltrait$scientific <- gsub(" ", "_", alltrait$scientific)                      # adding _ to species names for consistency
bruv_traits         <- alltrait[alltrait$scientific %in% bruv_species, ]

# how many species are we missing traits for?
length(bruv_species) - nrow(bruv_traits)
#17 = 16 spp and pempheris tomanagi


# look at number of maxn for individual species and see if there are a lot of them or not
# see if similar traits across genus that can be used
# use any of species names thought it was but couldnt be 100% sure 
# eg chromis westaustralis


# reduce traits data to just our covariates of interest - use ::summary to choose 
# possible traits to include,
###### feeding guild, fb.vulnerability, iucn ranking, body mass

summary(bruv_traits)
interesting_traits <- c("scientific", 
                        "rls_trophic_group")
bruv_traits <- bruv_traits[ , colnames(bruv_traits) %in% interesting_traits]

# overall traits table for all species but just our traits

alltrait_sub1 <- subset(alltrait, select = "scientific")
alltrait_sub3 <- subset(alltrait, select = "rls_trophic_group")
alltrait_sub0 <- cbind(alltrait_sub1, alltrait_sub3)

# clean up/remove species from traits data if they have NA in any trait info
bruv_traits <- na.omit(bruv_traits)
head(bruv_traits)

# clean traits data itself
#bruv_traits$feeding.guild <- tolower(bruv_traits$feeding.guild)                 # make all lowercase
#bruv_traits$feeding.guild <- sapply(bruv_traits$feeding.guild, 
              # FUN = function(x)paste(unique(unlist(strsplit(x, "/"))), 
                                      #collapse = " "))                          # remove punctuation and duplicates
#unique(bruv_traits$feeding.guild)

# add body data columns and tidy
bruv_traits <- merge(bruv_traits, bodydim[, c(4:6)], by = 'scientific')
bruv_traits$rls_trophic_group <- tolower(bruv_traits$rls_trophic_group)
bruv_traits$rls_trophic_group <- as.factor(bruv_traits$rls_trophic_group)
head(bruv_traits)
summary(bruv_traits)

# make species names row names
rownames(bruv_traits) <- bruv_traits$scientific
bruv_traits <- bruv_traits[ , -1]

# drop bruv maxn info for species that are lacking any traits
bruv_maxn_w <- bruv_maxn_w[ , colnames(bruv_maxn_w) %in% c(rownames(bruv_traits))]
dim(bruv_maxn_w)

# list species without traits :(
bruv_species_traits <- rownames(bruv_traits)

bruv_notrait <- bruv_species[(bruv_species %in% bruv_species_traits) == FALSE]
length(bruv_notrait)
bruv_notrait


#### wrangle spatial context data -----

str(bruv_covs)

bruv_xy <- subset(bruv_covs, select = sample:longitude)
rownames(bruv_xy) <- bruv_xy$sample
bruv_xy <- bruv_xy[,-1]

#### write RDS to preserve row names ------

# write to RDS to preserve row names 
saveRDS(bruv_maxn_w, "data/bruv_maxn_wide.rds")
saveRDS(bruv_covs,   "data/bruv_covariates_wide.rds")
saveRDS(bruv_traits, "data/bruv_traits_my_species.rds")
saveRDS(bruv_xy,     "data/bruv_xy.rds")


# fix!

# clear environment
rm(list=ls())
