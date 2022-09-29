 ###
# Project: Molly Thesis
# Data:    HMSC Outputs
# Task:    Investigate influence of covariates
# author:  Kingsley Griffin, Molly Gerrans
# date:    July 2022
##

library(Hmsc)

# load model outputs
model.directory <- "output/hmsc_model_data"

#thin 10
nChains = 4
samples = 1000
thin = 1000
filename = file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

## plot distribution of species-specific explanatory power (r2)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
x_expression <- expression(R^2)
hist(MF$SR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$SR2),2)), xlab = x_expression)

# note that for many species r2 was low but there are a few for which this 
# model has good explanatory power (>0.7)
# up to you if you want to use this but may be something we add to 
# supplementary materials

######## Explore parameter predictions ---------
#####gradient depth-----
Gradientd = constructGradient(m,focalVariable = "depth")
head(Gradientd$XDataNew)
predY = predict(m,
                XData = Gradientd$XDataNew, 
                studyDesign = Gradientd$studyDesignNew,
                ranLevels = Gradientd$rLNew,
                expected = TRUE)
par(mar=c(5,6,3,2))

#species
plotGradient(m, Gradientd, pred=predY, measure="Y", index = 9, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientd, pred=predY, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientd, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)


#####gradient slope-----
Gradients = constructGradient(m,focalVariable = "slope")
head(Gradients$XDataNew)
predYs = predict(m,
                XData = Gradients$XDataNew, 
                studyDesign = Gradients$studyDesignNew,
                ranLevels = Gradients$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradients, pred=predYs, measure="Y", index = 6, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradients, pred=predYs, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradients, pred=predYs, measure="T", index = 2, showData = TRUE, jigger = 0.1)

#####gradient detrended-----
Gradientde = constructGradient(m,focalVariable = "detrended")
head(Gradientde$XDataNew)
predYde = predict(m,
                XData = Gradientde$XDataNew, 
                studyDesign = Gradientde$studyDesignNew,
                ranLevels = Gradientde$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientde, pred=predYde, measure="Y", index = 6, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientde, pred=predYde, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientde, pred=predYde, measure="T", index = 5, showData = TRUE, jigger = 0.1)

#####gradient macroalgae-----
Gradientm = constructGradient(m,focalVariable = "macroalgae")
head(Gradientm$XDataNew)
predYm = predict(m,
                XData = Gradientm$XDataNew, 
                studyDesign = Gradientm$studyDesignNew,
                ranLevels = Gradientm$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientm, pred=predY, measure="Y", index = 8, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientm, pred=predY, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientm, pred=predY, measure="T", index = 5, showData = TRUE, jigger = 0.1)

#####gradient mean.relief------
Gradientmr = constructGradient(m,focalVariable = "mean.relief")
head(Gradientmr$XDataNew)
predYmr = predict(m,
                XData = Gradientmr$XDataNew, 
                studyDesign = Gradientmr$studyDesignNew,
                ranLevels = Gradientmr$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientmr, pred=predYmr, measure="Y", index = 9, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientmr, pred=predYmr, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientmr, pred=predYmr, measure="T", index = 5, showData = TRUE, jigger = 0.1)

#####gradient spongegarden-----
Gradientsg = constructGradient(m,focalVariable = "spongegarden")
head(Gradientsg$XDataNew)
predYsg = predict(m,
                XData = Gradientsg$XDataNew, 
                studyDesign = Gradientsg$studyDesignNew,
                ranLevels = Gradientsg$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientsg, pred=predYsg, measure="Y", index = 8, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientsg, pred=predYsg, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientsg, pred=predYsg, measure="T", index = 2, showData = TRUE, jigger = 0.1)

#####gradient consolidated-----
Gradientc = constructGradient(m,focalVariable = "consolidated")
head(Gradientc$XDataNew)
predYc = predict(m,
                XData = Gradientc$XDataNew, 
                studyDesign = Gradientc$studyDesignNew,
                ranLevels = Gradientc$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientc, pred=predYc, measure="Y", index = 8, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientc, pred=predYc, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientc, pred=predYc, measure="T", index = 2, showData = TRUE, jigger = 0.1)

#####gradient unconsolidated-----
Gradientuc = constructGradient(m,focalVariable = "unconsolidated")
head(Gradientuc$XDataNew)
predYuc = predict(m,
                XData = Gradientuc$XDataNew, 
                studyDesign = Gradientuc$studyDesignNew,
                ranLevels = Gradientuc$rLNew,
                expected = TRUE)

#species
plotGradient(m, Gradientuc, pred=predYuc, measure="Y", index = 8, showData = TRUE,  jigger = 0.1)
#totalcount
plotGradient(m, Gradientuc, pred=predYuc, measure="S",showData = TRUE, jigger = 0.1)
#traits
plotGradient(m, Gradientuc, pred=predYuc, measure="T", index = 2, showData = TRUE, jigger = 0.1)

