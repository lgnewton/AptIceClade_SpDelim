##### Aptostichus icenoglei SDMs for the 3 clades #####

##### Load all needed libraries #####
library(rJava)
library(ENMeval)
library(raster)
library(rgdal)
library(microbenchmark)
library(ENMTools)
library(ggplot2)
library(maps)
library(dismo)
#library(maxnet)

##### Trimming bioclimatic variable rasters to geographic area of interest #####
# Read in occurrence records; only 2 columns = Longitude and Latitude 
Aptice_occ_recs <- read.csv("~/Desktop/CA_Mygals_ClimateChange/Aptice_ClimateChange/Aptice_occ_recs.csv")
IceNorth_occ_recs <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceNorth_MaxEnt.csv")
IceCentral_occ_recs <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceCentral_MaxEnt.csv")
IceSouth_occ_recs <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceSouth_MaxEnt.csv")
IceSC_occ_recs <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceSC_MaxEnt.csv")

# Read in current bioclim variables from worldclim.org; 0.5 min resolution
Aptice_biovars <- getData("worldclim", var="bio", res=0.5, lon=-120, lat=32)
#Aptice_biovars_2.5 <- getData("worldclim", var="bio", res=2.5)
# Lists information - class, dimensions, extent
unlist(Aptice_biovars)

# Plot to visuaize what each climate layer of the tile looks like 
plot(Aptice_biovars)

# Defines the extent of the range for Aptostichus icenoglei
# You don't want to include too little or too much area (could potentially be misleading in subsequent niche-based distribution modeling)
Aptice_extent <- extent(-119.75, -115, 31.5, 35.25)
#Aptice_extent <- extent(-120, -113, 30, 37)
#Aptice_extent <- extent(x= c(min(Aptice_occ_recs$Longitude)-1, max(Aptice_occ_recs$Longitude)+1, min(Aptice_occ_recs$Latitude)-1, max(Aptice_occ_recs$Latitude))+1)

# Clip the tile to include only the extent of Apt. ice range 
Aptice_trim <- crop(Aptice_biovars, Aptice_extent)

# Plot to see if the study area is now limited to specified extent
plot(Aptice_trim)

# Change names; will be important if predicting future scenarios downstream
names(Aptice_trim) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')

# Identify highly correlated variables (Pearson correlation coefficient >= 0.8)
Aptice_trim_bc_corr <- raster.cor.matrix(Aptice_trim)
write.csv(Aptice_trim_bc_corr, file="~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_trim_bc_corr.csv")

# Only include variables not highly correlated
Aptice_trim <- Aptice_trim[[c("bio2", "bio4", "bio8", "bio9", "bio11", "bio14", "bio15", "bio17", "bio18", "bio19")]]

# Makes an asc file for each variable
setwd('~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/ClimLayers/')
raster_AptIce <- writeRaster(Aptice_trim, filename="Aptice_0.5.grd", format="ascii",bylayer=T, overwrite=T)

##### MaxEnt Models #####
#### Generate background (i.e., pseudo-absence) points ####
Aptice_background <- randomPoints(mask = Aptice_trim[[1]],     # Provides resolution of sampling points
                                  n = 10000,      # Number of random points
                                  ext = Aptice_extent, # Spatially restricts sampling
                                  extf = 1)
# had issue where background column names didn't match occurrence data; struggled to change column names in R so edited with Excel instead
setwd('~/Desktop/Ch3_IceClade_Dec2021/ENM_files/')
write.csv(as.data.frame(Aptice_background),file = "Aptice_background.csv")
#changed column names to Longitude and Latitude in Excel
#read back in
Aptice_background <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_background.csv")

### Run models with 5 diff regularization multipliers and 6 diff feature classes; method can change depending on distribution pattern (e.g., checkerboard = continuous; jackknife = disjunct) ###
setwd("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/")

#### Ice North SDM ####
IceNorth_eval <- ENMevaluate(IceNorth_occ_recs, Aptice_trim, Aptice_background, method=c('jackknife'), RMvalues=c(1, 1.5, 2, 2.5, 3), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm ='maxent.jar')

# Write model evaluations to csv
IceNorth_models <- as.data.frame(IceNorth_eval@results)
write.csv(IceNorth_models, file = "IceNorth_NoCorr_models.csv", row.names = F)

# Get best model
IceNorth_bestmod <- which(IceNorth_eval@results$delta.AICc==0)

# Write best model results to csv (I'm most interested in variable contribution since model eval stats are in previous csv)
write.csv(as.data.frame(IceNorth_eval@models[[IceNorth_bestmod]]@results), file="IceNorth_NoCorr_bestmodel.csv")

# Project the current SDM
IceNorth_SDM <- predict(Aptice_trim, IceNorth_eval@models[[IceNorth_bestmod]])
IceNorth_SDM_df <- as.data.frame(IceNorth_SDM, xy=T)

# Write current raster
writeRaster(IceNorth_SDM, file="IceNorth_NoCorr_SDM.asc", overwrite=T)

# Export current plot
IceNorth_plot <- ggplot() +
  geom_raster(data = IceNorth_SDM_df, aes(x = x, y = y, fill = layer)) + coord_quickmap() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), legend.title = element_blank()) + scale_fill_viridis_c(option = "viridis", na.value = "grey90") + xlab("Longitude") + ylab("Latitude")
pdf(file="IceNorth_NoCorr_plot.pdf")
IceNorth_plot
dev.off()

#### Ice Central SDM ####
IceCentral_eval <- ENMevaluate(IceCentral_occ_recs, Aptice_trim, Aptice_background, method=c('jackknife'), RMvalues=c(1, 1.5, 2, 2.5, 3), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm='maxent.jar')

# Write model evaluations to csv
IceCentral_models <- as.data.frame(IceCentral_eval@results)
write.csv(IceCentral_models, file = "IceCentral_NoCorr_models.csv", row.names = F)

# Get best model
IceCentral_bestmod <- which(IceCentral_eval@results$delta.AICc==0)

# Write best model results to csv (I'm most interested in variable contribution since model eval stats are in previous csv)
write.csv(as.data.frame(IceCentral_eval@models[[IceCentral_bestmod]]@results), file="IceCentral_NoCorr_bestmodel.csv")

# Project the current SDM
IceCentral_SDM <- predict(Aptice_trim, IceCentral_eval@models[[IceCentral_bestmod]])
IceCentral_SDM_df <- as.data.frame(IceCentral_SDM, xy=T)

# Write current raster
writeRaster(IceCentral_SDM, file="IceCentral_NoCorr_SDM.asc", overwrite=T)

# Export current plot
IceCentral_plot <- ggplot() +
  geom_raster(data = IceCentral_SDM_df, aes(x = x, y = y, fill = layer)) + coord_quickmap() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), legend.title = element_blank()) + scale_fill_viridis_c(option = "viridis", na.value = "grey90") + xlab("Longitude") + ylab("Latitude")
pdf(file="IceCentral_NoCorr_plot.pdf")
IceCentral_plot
dev.off()

#### Ice South SDM ####
IceSouth_eval <- ENMevaluate(IceSouth_occ_recs, Aptice_trim, Aptice_background, method=c('jackknife'), RMvalues=c(1, 1.5, 2, 2.5, 3), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm='maxent.jar')

# Write model evaluations to csv
IceSouth_models <- as.data.frame(IceSouth_eval@results)
write.csv(IceSouth_models, file = "IceSouth_NoCorr_models.csv", row.names = F)

# Get best model
IceSouth_bestmod <- which(IceSouth_eval@results$delta.AICc==0)

# Write best model results to csv (I'm most interested in variable contribution since model eval stats are in previous csv)
write.csv(as.data.frame(IceSouth_eval@models[[IceSouth_bestmod]]@results), file="IceSouth_NoCorr_bestmodel.csv")

# Project the current SDM
IceSouth_SDM <- predict(Aptice_trim, IceSouth_eval@models[[IceSouth_bestmod]])
IceSouth_SDM_df <- as.data.frame(IceSouth_SDM, xy=T)

# Write current raster
writeRaster(IceSouth_SDM, file="IceSouth_NoCorr_SDM.asc", overwrite=T)

# Export current plot
IceSouth_plot <- ggplot() +
  geom_raster(data = IceSouth_SDM_df, aes(x = x, y = y, fill = layer)) + coord_quickmap() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), legend.title = element_blank()) + scale_fill_viridis_c(option = "viridis", na.value = "grey90") + xlab("Longitude") + ylab("Latitude")
pdf(file="IceSouth_NoCorr_plot.pdf")
IceSouth_plot
dev.off()

#### Ice Central+South SDM ####
IceSC_eval <- ENMevaluate(IceSC_occ_recs, Aptice_trim, Aptice_background, method=c('jackknife'), RMvalues=c(1, 1.5, 2, 2.5, 3), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm='maxent.jar')

# Write model evaluations to csv
IceSC_models <- as.data.frame(IceSC_eval@results)
write.csv(IceSC_models, file = "IceSC_NoCorr_models.csv", row.names = F)

# Get best model
IceSC_bestmod <- which(IceSC_eval@results$delta.AICc==0)

# Write best model results to csv (I'm most interested in variable contribution since model eval stats are in previous csv)
write.csv(as.data.frame(IceSC_eval@models[[IceSC_bestmod]]@results), file="IceSC_NoCorr_bestmodel.csv")

# Project the current SDM
IceSC_SDM <- predict(Aptice_trim, IceSC_eval@models[[IceSC_bestmod]])
IceSC_SDM_df <- as.data.frame(IceSC_SDM, xy=T)

# Write current raster
writeRaster(IceSC_SDM, file="IceSC_NoCorr_SDM.asc", overwrite=T)

# Export current plot
IceSC_plot <- ggplot() +
  geom_raster(data = IceSC_SDM_df, aes(x = x, y = y, fill = layer)) + coord_quickmap() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), legend.title = element_blank()) + scale_fill_viridis_c(option = "viridis", na.value = "grey90") + xlab("Longitude") + ylab("Latitude")
pdf(file="IceSC_NoCorr_plot.pdf")
IceSC_plot
dev.off()

##### ENMTools Background Test #####
### upload SDMs ###
IceNorth_allbiovars_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_AllBioVars/IceNorth_allbiovars_SDM.asc")
IceCentral_allbiovars_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_AllBioVars/IceCentral_allbiovars_SDM.asc")
IceSouth_allbiovars_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_AllBioVars/IceSouth_allbiovars_SDM.asc")
IceSC_allbiovars_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_AllBioVars/IceSC_allbiovars_SDM.asc")

IceNorth_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/IceNorth_NoCorr_SDM.asc")
IceCentral_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/IceCentral_NoCorr_SDM.asc")
IceSouth_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/IceSouth_NoCorr_SDM.asc")
IceSC_SDM <- raster("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/IceSC_NoCorr_SDM.asc")

### Load environmental data generated from previous script run ###
env.files <- list.files(path = "~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/ClimLayers/", pattern = ".asc", full.names = TRUE)
env <- stack(env.files)
env <- setMinMax(env)

### Set coordinate system; otherwise will get error running background test ###
# First check your coordinate system in your asc files
crs(env)
# Add crs from asc files to your environmental data object
crs(env) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Look at loaded environmental data to make sure it looks correct
plot(env)

### Create species object for each lineage; can add each bit manually like here or combine into 1 line of code ###
# Central #
IceCentral_spo <- enmtools.species()
IceCentral_spo$species.name <- "IceCentral"
IceCentral_spo$presence.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceCentral_MaxEnt.csv")
IceCentral_spo$range <- background.raster.buffer(IceCentral_spo$presence.points, 50000, mask = env)
IceCentral_spo$background.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/NicheSimilarity/ArcGIS_BackgroundPoints/Central_BGpoints_RasterPolygonsHSS75.csv")

# South #
Ice_South_spo <- enmtools.species()
Ice_South_spo$species.name <- "IceSouth"
Ice_South_spo$presence.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceSouth_MaxEnt.csv")
Ice_South_spo$range <- background.raster.buffer(Ice_South_spo$presence.points, 50000, mask = env)
Ice_South_spo$background.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/NicheSimilarity/ArcGIS_BackgroundPoints/South_BGpoints_RasterPolygonsHSS75.csv")

# North #
IceNorth_spo <- enmtools.species()
IceNorth_spo$species.name <- "IceNorth"
IceNorth_spo$presence.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceNorth_MaxEnt.csv")
IceNorth_spo$range <- background.raster.buffer(IceNorth_spo$presence.points, 50000, mask = env)
IceNorth_spo$background.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/NicheSimilarity/ArcGIS_BackgroundPoints/North_BGpoints_RasterPolygonsHSS75.csv")

# South + Central #
Ice_SC_spo <- enmtools.species()
Ice_SC_spo$species.name <- "IceSC"
Ice_SC_spo$presence.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/AptIceSC_MaxEnt.csv")
Ice_SC_spo$range <- background.raster.buffer(Ice_SC_spo$presence.points, 50000, mask = env)
Ice_SC_spo$background.points <- read.csv("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/NicheSimilarity/ArcGIS_BackgroundPoints/SC_BGpoints_RasterPolygonsHSS75.csv")

### Background tests with varying km polygons ###
#asymmetric vs symmetric = the background test here can be conducted against a null hypothesis that is generated from "asymmetric" (species.1 vs species.2 background) or "symmetric" (species.1 background vs. species.2 background) comparisons

#### South vs Central ####
setwd("~/Desktop/Ch3_IceClade_Dec2021/ENM_files/Aptice_NoCorr/NicheSimilarity/ENMTools_R/CvS_SvC/")

bg_SvC <- background.test(species.1 = Ice_South_spo, species.2 = IceCentral_spo, env = env, type = "mx", nreps = 100, test.type = "asymmetric" )
bg_CvS <- background.test(species.1 = IceCentral_spo, species.2 = Ice_South_spo, env = env, type = "mx", nreps = 100, test.type = "asymmetric" )


write.csv(bg_CvS$reps.overlap,file="CvS_RasPolysHSS75.csv")
write.csv(bg_SvC$reps.overlap,file='SvC_RasPolysHSS75.csv')

png(filename="bg_SvC_RasPolysHSS75.png", height=800, width=700)
bg_SvC
dev.off()

png(filename="bg_CvS_RasPolysHSS75.png", height=800, width=700)
bg_CvS
dev.off()

bg_SvC$p.values
bg_CvS$p.values

#### North vs South+Central ####
bg_NvSC <- background.test(species.1 = IceNorth_spo, species.2 = Ice_SC_spo, env = env, type = "mx", nreps = 100, test.type = "asymmetric" )
bg_SCvN <- background.test(species.1 = Ice_SC_spo, species.2 = IceNorth_spo, env = env, type = "mx", nreps = 100, test.type = "asymmetric" )

write.csv(bg_NvSC$reps.overlap,file="NvSC_RasPolysHSS75.csv")
write.csv(bg_SCvN$reps.overlap,file="SCvN_RasPolysHSS75.csv")

png(filename="NvSC_RasPolysHSS75.png", height=800, width=700)
bg_NvSC
dev.off()

png(filename="SCvN_RasPolysHSS75.png", height=800, width=700)
bg_SCvN
dev.off()

bg_NvSC$p.values
bg_SCvN$p.values


