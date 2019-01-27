
laselva <- read.csv("laselva_vcurves.csv") %>%
  mutate(curve_id = paste(genusspecies,individual,sep="-"))
var_names <- c(K = "K", WP = "MPa")

# latest from bitbucket: important changes will be pushed to CRAN soon
if(FALSE){
  remotes::install_bitbucket("remkoduursma/fitplc")
}
library(fitplc)
if(packageVersion("fitplc") < "1.3"){
  stop("did I not tell you to update?")
}


# from plot
library(ggplot2)
ggplot(laselva, aes(MPa,K)) +
  geom_line() +
  facet_wrap(~curve_id) +
  theme_minimal()
# it is obvious that K does not plateau at high water potentials,
# this means a typical S-curve cannot be fit

# global fit
# for fitcond you have to specify what you think is the asymptoting value,
# i.e. K at a water potential of zero. There can be lots of debate about this,
# and the water physiologists often disagree!
# HEre I use rescale_Px=TRUE, see ?fitcond, and fit a non-parametric curve
# because your data are not S-shaped
f1 <- fitcond(laselva, model = "loess", rescale_Px=TRUE, varnames=var_names)

# see ?plot.fitplc for options
plot(f1)


# now by species
# this takes a while because every curve is bootstrapped
# (I think this overestimates the variance of these fits by the way,
# i can look into it)
f2 <- fitconds(laselva, "curve_id",
               model = "loess", rescale_Px=TRUE, varnames=var_names)

par(mfrow=c(3,4))
plot(f2)

