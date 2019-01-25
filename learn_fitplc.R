 
laselva_vc <- read.csv("laselva_vcurves.csv")

##REMKO, can you run through some examples of how to get the most out of fitplc
## We de-gas the stipes so Kmax is 0.0 MPa and we start with K 
## Getting PLC first doesnt seem necessary with fitplcs but I havent gotten it to work

##I was having issues with using raw K values:
#Error in fitplc(dfr, varnames = c(PLC = "plc", WP = varnames[["WP"]]),  : 
#                  Missing values found in PLC or WP - remove first!

#the much larger dataset will be split into terrestrial and epiphtic factors.
#is there a way in fitplc to plot the two factor groups with CI's?
#could randomly assign some factor level below to give an example of needed

library(fitplc)

#test
test <- laselva_vc[1:7,]
plc_test <- fitcond(test, varnames=c(K="K",
                                     WP="MPa"),  WP_Kmax = 0.0)
##why the missing values error? I thought K was used in fitcond...

laselva_fits <- fitconds(laselva_vc, varnames=c(K="conductivity",
                         WP="MPa"), group="sample_id", WP_Kmax = 0.0)

plot(laselva_fits, onepanel=TRUE, plotci=FALSE, px_ci="none", pxlinecol="dimgrey")