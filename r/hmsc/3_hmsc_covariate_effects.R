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
thin = 100
filename = file.path(paste(model.directory), paste0("model_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

######## Explore parameter predictions ---------

#gradient depth
Gradientd = constructGradient(m,focalVariable = "depth")

head(Gradientd$XDataNew)

predY = predict(m,
                XData = Gradientd$XDataNew, 
                studyDesign = Gradientd$studyDesignNew,
                ranLevels = Gradientd$rLNew,
                expected = TRUE)
# predY = predict(m, Gradient=Gradientd, expected = TRUE)


plotGradient(m,
             Gradientd,
             pred=predY,
             measure="S",
             showData = TRUE)

plotGradient(m, Gradientd, pred=predY, measure="Y", index = 50, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientd, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)


#gradient location
Gradientl = constructGradient(m,focalVariable = "macroalgae")

head(Gradientl$XDataNew)

predY = predict(m,
                XData = Gradientl$XDataNew, 
                studyDesign = Gradientl$studyDesignNew,
                ranLevels = Gradientl$rLNew,
                expected = TRUE)
# predY = predict(m, Gradient=Gradientl, expected = TRUE)


plotGradient(m,
             Gradientl,
             pred=predY,
             measure="S",
             showData = TRUE)

plotGradient(m, Gradientl, pred=predY, measure="Y", index = 17, showData = TRUE,  jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="S", showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="T", index = 2, showData = TRUE, jigger = 0.1)
plotGradient(m, Gradientl, pred=predY, measure="T", index = 4, showData = TRUE, jigger = 0.1)

