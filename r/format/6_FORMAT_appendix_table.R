###
# Project: Molly Thesis
# Data:    BRUVS Data from Abrolhos
# Task:    Appendix table 
# author:  Kingsley Griffin, Molly Gerrans
# date:    September 2022
##

library(devtools)
devtools::install_github('UWAMEGFisheries/GlobalArchive')
library(GlobalArchive)

# To connect to GlobalArchive


library(httpuv)

install.packages('googlesheets4')

library(googlesheets4)

# To tidy data

library(tidyr)

library(plyr)

library(dplyr)

library(stringr)

library(readr)

library(ggplot2)



# Read in life history

url <- "https://docs.google.com/spreadsheets/u/1/d/1SMLvR9t8_F-gXapR2EemQMEPSw_bUbPLcXd3lJ5g5Bo/edit?usp=drive_web&ouid=100340010373917954123"


master<-read_sheet(url)%>%
  
  ga.clean.names()%>%
  
  filter(grepl('Australia', global.region))%>%
  
  filter(grepl('NW', marine.region))%>%
  
  dplyr::select(family,genus,species,iucn.ranking,fishing.mortality,fishing.type,australian.common.name)%>%
  
  distinct()%>%
  
  glimpse()



names(master)



# Create Species list ----



species.table <- bruv_maxn%>%
  
  group_by(family,genus,species,scientific)%>%
  
  summarise_at(vars(matches("maxn")),funs(sum,mean,sd,se=sd(.)/sqrt(n())))%>%
  
  ungroup()%>%
  
  mutate(mean=round(mean,digits=2))%>%
  
  mutate(sd=round(sd,digits=2))%>%
  
  mutate(se=round(se,digits=2))%>%
  
  mutate(genus.species=paste(genus,species,sep=" "))%>%
  
  arrange(family)%>%
  
  left_join(master)%>%
  
  dplyr::select(-c(scientific))%>%
  
  dplyr::mutate(mean.relative.abundance.per.deployment.plus.minus.SE=paste(mean,"+/-",se,sep=" "))%>%
  
  dplyr::rename(total.relative.abundance = sum)%>%
  
  ungroup()



unique(species.table$fishing.type)



cleaned<-species.table%>%
  
  dplyr::select(family,genus.species,australian.common.name,fishing.type,iucn.ranking,mean.relative.abundance.per.deployment.plus.minus.SE,total.relative.abundance)%>%
  
  ## fix up variables
  
  mutate(fishing.type=ifelse(fishing.type%in%c("C/R","C","B/C"),"Commercial",""))

## Make names nicer for table

unique(cleaned$fishing.type)

write.csv(cleaned, "plots/appendix_table.csv")
