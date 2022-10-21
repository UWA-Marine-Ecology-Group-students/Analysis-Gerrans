###
# Project: Molly Thesis
# Data:    HMSC model data
# Task:    Plot beta effects
# author:  Kingsley Griffin, Molly Gerrans
# date:    Oct 2022
##

library(Hmsc)
library(dplyr)
library(reshape2)
library(ggplot2)

# load HMSC model data
model.directory <- "output/hmsc_model_data"

nChains = 4
samples = 1000
thin = 1000
filename=file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

# extract beta mean and CI in a single dataframe
postBeta         <- getPostEstimate(m, parName = "Beta", q = c(0.05, 0.95))
bmeans           <- melt(as.data.frame(postBeta$mean))
bmeans$covariate <- c("Intercept", "Depth", "Macroalgae", "Mean Relief", 
                      "Spongegarden", "Consolidated", "Unconsolidated", 
                      "Slope", "Detrended")
colnames(bmeans) <- c("species", "mean", "covariate")
bmeans$spcov <- interaction(bmeans$species, bmeans$covariate)
head(bmeans)
dim(bmeans)

quants   <- as.data.frame(unlist(postBeta$q))
quants$q <- row.names(quants)
quants   <- dcast(melt(quants), variable ~ q)
head(quants)
dim(quants)
quants$covariate <- rep(c("Intercept", "Depth", "Macroalgae", "Mean Relief", 
                        "Spongegarden", "Consolidated", "Unconsolidated", 
                        "Slope", "Detrended"), 21)
quants$species   <- rep(unique(bmeans$species), each = 9)
quants$spcov <- interaction(quants$species, quants$covariate)
head(quants)

beta_ci <- merge(bmeans, quants, by = "spcov")
beta_ci <- beta_ci[, -c(1, 5, 8, 9)]
colnames(beta_ci) <- c("species", "mean", "covariate", "lower", "upper")
beta_ci$sig       <- ifelse(beta_ci$lower < 0 & beta_ci$upper > 0, 0.8, 0.9)
head(beta_ci)

ggplot(beta_ci[beta_ci$covariate != "Intercept", ], 
       aes(species, mean, alpha = as.numeric(sig))) + 
  geom_hline(yintercept = 0, lty = 3, alpha = 0.8) +
  geom_point() +
  coord_flip() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.125) +
  guides(alpha = "none") +  
  labs(x = NULL, y = "Estimate") +
  facet_grid( ~ covariate, scales = "free") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        plot.title = element_text(size=20)) +
  theme(axis.text.y = element_text(face = "italic"))+
  theme(text = element_text(size = 24))
