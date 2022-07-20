###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    Wrangling data for HMSC formats
# author:  Kingsley Griffin, Molly Gerrans
# date:    May 2022
##

library(dplyr)
library(reshape2)
library(raster)
library(stars)
library(starsExtra)

#### wrangle maxn data into wide format ----

# read in maxn data
bruv_maxn <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.maxn.csv")
colnames(bruv_maxn)

# filter out rare species
bruv_maxn$pres_abs <- ifelse(bruv_maxn$maxn > 0, 1, 0)

# calculate number of sightings (sum of 1's from above)
for(speciesi in unique(bruv_maxn$scientific)){
  bruv_maxn$sp_n_sightings[bruv_maxn$scientific == speciesi] <- 
    sum(bruv_maxn$pres_abs[bruv_maxn$scientific == speciesi])
}
head(bruv_maxn)

# how many unique species do we start with?
length(unique(bruv_maxn$scientific))

# how many if we cut off at n sightings?
length(unique(bruv_maxn$scientific[bruv_maxn$sp_n_sightings > 5]))

sptokeep <- unique(bruv_maxn$scientific[bruv_maxn$sp_n_sightings > 5])
# exclude low sighting species

bruv_maxn <- bruv_maxn[bruv_maxn$scientific %in% sptokeep, ]

# remake scientific name with just genus and species for matching up with other dfs
bruv_maxn$genus_species <- c(paste(bruv_maxn$genus, 
                                   bruv_maxn$species, sep = "_"))
bruv_maxn$genus_species <- gsub(" ", "_", bruv_maxn$genus_species) 

head(bruv_maxn$genus_species)

bruv_maxn_w  <- reshape2::dcast(bruv_maxn[, c(1, 3, 22)],     
                                sample ~ genus_species, 
                                fun = sum, value.var = "maxn")      # columns included are sample, maxn and genus_species
rownames(bruv_maxn_w) <- bruv_maxn_w$sample                         # make sample ids row names
rownames(bruv_maxn_w) <- gsub("\\.", "_", bruv_maxn_w$sample)
bruv_maxn_w <- bruv_maxn_w[, -1]                                    # drop sample column
head(bruv_maxn_w)


#### wrangle habitat and environmental covariate info into wide format ----

bruv_meta    <- read.csv("data/raw/em export/2021-05_Abrolhos_stereo-BRUVs_Metadata.csv.csv")
bruv_habitat <- read.csv("data/tidy/2021-05_Abrolhos_BRUVs_random-points_percent-cover_broad.habitat.csv")
colnames(bruv_meta)                                                             # the columns of the original data that we can choose covariates from
bruv_covs <- dplyr::select(bruv_meta, c("Sample", "Latitude", "Longitude", 
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

# get the bathymetric derivatives
# make drop data into a spatial dataset
bruvcov_sp <- SpatialPointsDataFrame(coords = cbind(bruv_covs$longitude, 
                                                    bruv_covs$latitude), 
                                     data = bruv_covs)
plot(bruvcov_sp)

# bring in bathymetry
bathy              <- raster("data/spatial/rasters/WA_500m_bathy.tif")
wgscrs             <- CRS("+proj=longlat +datum=WGS84")
proj4string(bathy) <- wgscrs
bathy              <- crop(bathy, buffer(bruvcov_sp, 5000))
names(bathy)       <- "depth"
plot(bathy)
plot(bruvcov_sp, add = T)

# generate terrain variables for project area
slope <- terrain(bathy, opt = "slope", unit = "degrees", neighbours = 8)
zstar <- st_as_stars(bathy)                                                     # convert to stars obj
detre <- detrend(zstar, parallel = 8)                                           # detrend bathymetry
detre <- as(object = detre, Class = "Raster")                                   # back to raster
names(detre)     <- c("detrended", "lineartrend")
envcov           <- stack(bathy, slope, detre[[1]])
plot(envcov)

bruv_covs <- cbind(bruv_covs, extract(envcov[[2:3]], bruvcov_sp))
head(bruv_covs)

saveRDS(envcov, "data/spatial/rasters/bathy_derivatives.rds")

#### generate traits table including each species ----

# read body length and mass
# bodylength <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.length.csv") # both measures appear in mass table
bodymass   <- read.csv("data/tidy/2021-05_Abrolhos_stereo-BRUVs_complete.mass.csv")

# clean up/remove species from length/mass data if they have NA in any measurement info
#bodylength <- na.omit(bodylength)
bodymass   <- na.omit(bodymass)

# calculate mean body measures per fish
bodydim <- summarise(group_by(bodymass, family, genus, species),
                     meanlength = sum(length)/sum(number),
                     meanmass = sum(mass.g)/sum(number))
bodydim$scientific  <- interaction(bodydim$genus, bodydim$species, sep = "_")
bodydim$scientific <- gsub(" ", "_", bodydim$scientific)     
head(bodydim)

# check correlation among length and mass
plot(log(bodydim$meanlength), log(bodydim$meanmass))
cor(bodydim$meanlength, bodydim$meanmass) 

# read in full traits table, clean column names and shorten from entire sheet to just our species
alltrait   <- read.csv("data/traits/Australia.life.history - australia.life.history.csv")
sptraits <- read.csv("data/traits/traits_spp.csv")
colnames(sptraits) <- gsub("\\.", "_", colnames(sptraits))
colnames(sptraits) <- tolower(colnames(sptraits))
colnames(sptraits)

bruv_species        <- colnames(bruv_maxn_w)     # species from maxn column names
head(sptraits)
# View(bruv_species)
sptraits$scientific <- gsub(" ", "_", sptraits$scientific)                      # adding _ to species names for consistency
bruv_traits         <- sptraits[sptraits$scientific %in% bruv_species, ]

# how many species are we missing traits for?
length(bruv_species) - nrow(bruv_traits)
#17 = 16 spp and pempheris tomanagi


# look at number of maxn for individual species and see if there are a lot of them or not
# see if similar traits across genus that can be used
# use any of species names thought it was but couldnt be 100% sure 
# eg chromis westaustralis

# merge body length and mass data - only one species taken from RLS
bruv_traits <- merge(bruv_traits, bodydim[, 4:6], by = "scientific", all.x = TRUE)
head(bruv_traits)

bruv_traits$splength <- ifelse(is.na(bruv_traits$meanlength), 
                               as.numeric(bruv_traits$body_length), 
                               bruv_traits$meanlength)
bruv_traits$spmass <- ifelse(is.na(bruv_traits$meanmass), 
                               as.numeric(bruv_traits$body_mass), 
                               bruv_traits$meanmass)

# reduce traits data to just our covariates of interest - use ::summary to choose 
# possible traits to include,
###### feeding guild, fb.vulnerability, iucn ranking, body mass

summary(bruv_traits)
head(bruv_traits)

interesting_traits <- c("scientific", 
                        "trophic_group",
                        "complexity",
                        "substrate_group",
                        "water_column",
                        "night_day",
                        "spmass",
                        "splength")

bruv_traits <- bruv_traits[ , colnames(bruv_traits) %in% interesting_traits]

head(bruv_traits)
summary(bruv_traits)

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

bruv_traits$trophic_group   <- tolower(bruv_traits$trophic_group)
bruv_traits$substrate_group <- tolower(bruv_traits$substrate_group)
bruv_traits$complexity      <- tolower(bruv_traits$complexity)
bruv_traits$water_column    <- tolower(bruv_traits$water_column)
bruv_traits$night_day       <- tolower(bruv_traits$night_day)

bruv_traits$trophic_group   <- as.factor(bruv_traits$trophic_group)
bruv_traits$substrate_group <- as.factor(bruv_traits$substrate_group)
bruv_traits$complexity      <- as.factor(bruv_traits$complexity)
bruv_traits$water_column    <- as.factor(bruv_traits$water_column)
bruv_traits$night_day       <- as.factor(bruv_traits$night_day)
head(bruv_traits)
summary(bruv_traits)

# make species names row names
rownames(bruv_traits) <- bruv_traits$scientific
bruv_traits           <- bruv_traits[ , -1]

# drop bruv maxn info for species that are lacking any traits
bruv_maxn_w <- bruv_maxn_w[ , colnames(bruv_maxn_w) %in% c(rownames(bruv_traits))]
dim(bruv_maxn_w)

# list species without traits :(
bruv_species_traits <- rownames(bruv_traits)
bruv_notrait        <- bruv_species[(bruv_species %in% bruv_species_traits) == FALSE]
length(bruv_notrait)
bruv_notrait

#### wrangle spatial context data -----
str(bruv_covs)

bruv_xy           <- subset(bruv_covs, select = sample:longitude)
rownames(bruv_xy) <- bruv_xy$sample
bruv_xy           <- bruv_xy[, -1]

#### write RDS to preserve row names ------
saveRDS(bruv_maxn_w, "data/bruv_maxn_wide.rds")
saveRDS(bruv_covs,   "data/bruv_covariates_wide.rds")
saveRDS(bruv_traits, "data/bruv_traits_my_species.rds")
saveRDS(bruv_xy,     "data/bruv_xy.rds")


# clear environment
rm(list=ls())
