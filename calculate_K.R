## laselva vulnerability curve format
library(plyr)

## read in datasheet with flow parameters

laselva_times <- read.csv("vcurve_raw_data/laselva/vcurve_times_laselva.csv")
  laselva_times$water_potential <- round(laselva_times$water_potential,1)
  laselva_times$sample_id <- with(laselva_times,
                                  paste(species, individual, 
                                  format(water_potential,digits=1), sep="_"))


##Read in flow meter data files into a list and create unique ID to match times dfr

#create a list of clean file names that to use latter
vcurves <- list.files(path="vcurve_raw_data/laselva/",
                      pattern="mpa",full.names=TRUE)
#extract the genusspecies and treatment info
vcurves_names <- gsub("vcurve_raw_data/laselva/", "", vcurves)
vcurves_names <- gsub("_mpa.csv", "", vcurves_names)

##read in all data files from laselva
#skip first 14 rows from sensirion flow meter csv files
vcurve_files <- llply(list.files(path="vcurve_raw_data/laselva/",
                                 pattern="mpa",full.names=TRUE),
                                 function(filename){
                                 dat=read.csv(filename, header=TRUE,skip=14)
                                 })

#add new variable with unique ID
for(i in seq_along(vcurve_files)){
  vcurve_files[[i]]$sample_id <- vcurves_names[i]
}

## function to calculate conductivity from stipe flow rates--------
vcurve_function <- function(dfr, timesdfr,
                            massflow_constant1 = -0.2323,
                            massflow_constant2 = 1002.9,
                            acceleratebygrav = 0.09806) {
  
  x <- data.frame(dfr)
  
  #subset times dataframe to match sample id of an individual curve
  y <- timesdfr[timesdfr$sample_id == unique(x$sample_id),]
  
  #density of H20 calculation
  y$h20_dens <- (massflow_constant1 * y$air_temp_C )+ massflow_constant2
  
  #isolate backgrounds (x2) and flow from timesdfr to trim flow data sets
  #ackward variable names are raw names from flow meter
  background1 <- dfr[dfr$Relative.Time.s. > (y$back_first_initial-1)
                     & dfr$Relative.Time.s. < y$back_first_final ,]
  
  background2 <- x[x$Relative.Time.s. > (y$back_second_initial-1)
                   & x$Relative.Time.s. < y$back_second_final ,]
  
  flow <- x[x$Relative.Time.s. > (y$flow_initial-1)
            & x$Relative.Time.s. < y$flow_final ,]
  
  #mass flow unit conversion for each set of flow data
  background1$background1_mgsec <- (((background1$Flow..ul.min. * 10^-9)
                                     * y$h20_dens) * 1000000)/60
  
  
  background2$background2_mgsec <- (((background2$Flow..ul.min. * 10^-9)
                                     * y$h20_dens) * 1000000)/60
  
  flow$flow_mgsec <- (((flow$Flow..ul.min. * 10^-9)
                       * y$h20_dens) * 1000000)/60
  
  ##average mass flow rates or each data set
  background1_mean <- mean(background1$background1_mgsec)
  background2_mean <- mean(background2$background2_mgsec)
  flow_mean <- mean(flow$flow_mgsec)
  
  #average backgrounds
  background_agg <- (background1_mean + background2_mean)/2
  
  #calculate corrected mass flow rate (mg sec-1)
  corr_flow_rate <- flow_mean - background_agg
  
  #convert to conductivity
  pressurehead <- y$pressure_head_cm * acceleratebygrav
  
  conductivity <- corr_flow_rate/(pressurehead/y$stipe_length_mm) #stipe must be in mm
  #mg mm KPa-1s-1
  
  id_cond <- data.frame(sample_id = unique(x$sample_id), 
                        genusspecies = y$species, individual = y$individual,
                        K = conductivity, MPa = y$water_potential)
  return(id_cond)
}

#test function with simple data frames
testdata <- vcurve_files[[6]]
test1 <- vcurve_function(dfr = testdata, timesdfr = laselva_times)

#test function with larger list
laselva_cond <- lapply(vcurve_files, function(x) 
  vcurve_function(x,  timesdfr=laselva_times))

laselva_cond_genspe <- rbind.fill(laselva_cond)

write.csv(laselva_cond_genspe, "laselva_vcurves.csv", row.names = FALSE)

