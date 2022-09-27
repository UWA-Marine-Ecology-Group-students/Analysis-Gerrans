
mixing <- as.mcmc(m, parameters = "paramX")                                 # just mixing object
head(mixing)


mix_df <- as.data.frame(mixing$Y)                                               # Convert the mixing object to a matrix
dat <- melt(mixingDF)

# boxplots of estimates other than intercepts
p <- ggplot(dat[! dat$variable %in% unique(dat$variable)[1:5], ], aes(variable, value))
p + geom_boxplot() + geom_hline(yintercept = 0, lty = 3) + coord_flip()

# table 95% credible intervals for species~X parameters ----
average   <- apply(model$results$estimation$paramX, 1:2, mean)
CI.025    <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.025)
CI.975    <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.975)
CI        <- cbind(as.vector(CI.025),as.vector(CI.975))
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)), 
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average), 
                                     each = nrow(average)), "_", 
                                 rep(rownames(average), 
                                     ncol(average)), sep="")
paramXCITable

ptab <- as.data.frame(paramXCITable)
ptab$response <- sapply(strsplit(rownames(ptab), split = "_"), tail, 1)
ptab$predictor <- "NA"
for (i in 1:nrow(ptab)){
  ptab$predictor[i] <- removeWords(rownames(ptab)[i], ptab$response[i], "_")
}
ptab$sig <- "NA"
ptab$sig[(ptab$lowerCI < 0 & ptab$upperCI < 0)|
           (ptab$lowerCI > 0 & ptab$upperCI > 0)] <- 1
ptab$sig[ptab$sig != 1] <- 0

ptab$predictor <- recode(ptab$predictor, depth5m = "Depth", 
                         distmouth5m = "Estuary", urchins5m = "Urchins",
                         anchor30m = "Anchoring", recactivity30m = "Recreation",
                         fishers30m = "Fishing", auto1 = "Autocorrelation")
ptab$response <- recode(ptab$response, algcanop = "Other Canopy", 
                        algeckl = "Ecklonia", algunders = "Understorey",
                        algenc = "Encrusting", algturf = "Turfs")

# saveRDS(ptab, 'output/hmsc_summary_df.rds')

ggplot(ptab[ptab$predictor != "(Intercept)", ], aes(predictor, paramX, alpha = as.numeric(sig)+4/5)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.125) +
  geom_hline(yintercept = 0, lty = 3, alpha = 3/4) +
  guides(size = FALSE, shape = FALSE, alpha = FALSE) +  
  labs(colour = "estimate                    ") + 
  facet_grid( ~ response, scales = "free_x") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_blank(), 
        plot.title = element_text(hjust = -0.15, size=11),
        strip.background = element_blank())
