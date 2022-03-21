### Merge EventMeasure database output tables into maxn and length files
### OBJECTIVES ###
# combine database tables into single Metadata, MaxN and Length files for subsequent validation and data analysis.

### Please forward any updates and improvements to tim.langlois@uwa.edu.au & brooke.gibbons@uwa.edu.au or raise an issue in the "globalarchive-query" GitHub repository

rm(list=ls()) # Clear memory

## Load Libraries ----
# To connect to GlobalArchive
library(devtools)
install_github("UWAMEGFisheries/GlobalArchive") #to check for updates
library(GlobalArchive)
# To connect to GitHub
library(RCurl)
library(R.utils)
# To tidy data
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

## Set Study Name ----
# Change this to suit your study name. This will also be the prefix on your final saved files.
study <- "database.tables.example" 

## Save directory name to use later----
download.dir <- paste(getwd(),"/data/raw/em export",sep="/")

# Combine all data----
# The below code will find all files that have the same ending (e.g. "_Metadata.csv") and bind them together.
# The end product is three data frames; metadata, maxn and length.

# Metadata ----
# You will need a metadata file.
# If using a .csv the file name MUST end in "_Metadata.csv"
# Need to match the global archive format
# See the user manual: https://globalarchivemanual.github.io/ for the correct format

# For csv file ----
metadata <- ga.list.files("_Metadata.csv") %>% # list all files ending in "_Metadata.csv"
  purrr::map_df(~ga.read.files_em.csv(.)) %>% # combine into dataframe
  dplyr::mutate(campaignid = study) %>%
  dplyr::select(campaignid,sample,latitude,longitude,date,time,location,status,site,depth,observer,successful.count,successful.length) %>% # This line ONLY keep the 15 columns listed. Remove or turn this line off to keep all columns (Turn off with a # at the front).
  glimpse()

unique(metadata$campaignid) # check the number of campaigns in metadata, and the campaign name

write.csv(metadata, paste("data/staging/", study,"_metadata.csv", sep = ""), row.names = FALSE)

## Combine Points and Count files into maxn ----
maxn<-ga.create.em.maxn() %>%
  dplyr::inner_join(metadata) %>%
  dplyr::filter(successful.count == "Yes") %>%
  dplyr::filter(maxn > 0)

unique(maxn$sample)

# Save MaxN file ----
write.csv(maxn, paste("data/staging/", study,"_maxn.csv", sep=""), row.names = FALSE)

## Combine Length, Lengths and 3D point files into length3dpoints----
length3dpoints<-ga.create.em.length3dpoints() %>%
  dplyr::select(-c(time, comment)) %>% # take time out as there is also a time column in the metadata
  dplyr::inner_join(metadata) %>%
  dplyr::filter(successful.length == "Yes") %>%
  glimpse()

## Save length files ----
setwd(staging.dir)
write.csv(length3dpoints, paste("data/staging/", study,"_length3dpoints.csv", sep=""), row.names = FALSE)