###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## RESISTANCE DISTANCES ----

## Clean Global Environment 
rm(list=ls())

## Required R libraries
library(gdistance)
library(sp)
library(raster)
library(rgdal)
library(geosphere)
library(rgeos)
library(tidyverse)

## Set projection settings before running the script
A=make_EPSG()
LL_prj <- as.character(subset(A, A$code=="4326")[3]) 
UTM_prj <- as.character(subset(A, A$code=="5383")[3]) 

## Set project name
project_name = "charinus_SS"

## COORDINATES ----
## Load
coordinates <- read.csv("./coords/cha/coords_charinus_93ind.csv", h=T)

## Set coords projection  
coord_UTM <- coordinates[, 3:4]
coordinates(coord_UTM) <- coord_UTM
projection(coord_UTM) <- crs(UTM_prj)
plot(coord_UTM, axes=T)
  

## 0. GEOGRAPHIC DISTANCE ----
var_name = "geo_distance"

r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

# Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16
# directions = 16 combine king's and knight's moves - The section on distance transforms in de Smith, Goodchild, and Longley (2009) also discusses 16-cell neighborhoods.
# Connecting in 16 directions may increase the accuracy of the calculations.]

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)


## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance 
# (the number of cells transversed during a random walk from a 
# cell i on the grid to a cell j and back to i (Chandra et al. 1996))

rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)


## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
          melt %>%
          na.omit %>%
          setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End


## 1. ELEVATION ----
var_name = "elevation"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End


### 2. ROUGHNESS ----
var_name = "roughness"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)


## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End


## 3. GEOENVIRONMENT ----
var_name = "geoenvironment"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End


### 4. FOREST COVER ----
var_name = "forestcover"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End



### 5. BIO10 ----
var_name = "bio10"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)


## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End



## 6. BIO15 ----
var_name = "bio15"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End



### 7. BIO3 ----
var_name = "bio3"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End



### 8. PETdq ----
var_name = "PETdq"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End




### 9. PETs ----
var_name = "PETs"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End




### 10. TopoWet ----
var_name = "topoWet"

## Load raster
r <- raster(paste0("./rasters/resist_surfaces/cha_SS/", var_name, "_", project_name, ".tif"))
crs(r)

## Check plots
plot(r)
points(coord_UTM)

## Create transition layer

## Chose  the number of directions to build a neighborhood graph      
chosen_directions = 16

## Build transition layer
tr_lyr <- transition(r, transitionFunction = mean, directions = chosen_directions) 

## Check transition layer
plot(raster(tr_lyr), main= paste0(var_name), xlab="Longitude", ylab="Latitude") 
points(coord_UTM)

## Correction for random walk - commute distance
tr_lyrR <- geoCorrection(tr_lyr, "r", scl=TRUE)

## Calculate commute distance
rw_Dist <- as.matrix(commuteDistance(tr_lyrR, coord_UTM)/1000)

## Convert distance matrix into data frame for analyses
rw_Dist = rw_Dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("cd_", var_name)))

rw_Dist = rw_Dist  %>% 
  filter(ID1!=ID2)

## Check
head(rw_Dist)
hist(rw_Dist[,3])
nrow(rw_Dist)

## Save
var_name
write.csv(rw_Dist, paste0("./5-LandGenAnalysis/results/resist_dist/cha_SS/resdist_", var_name,".csv"), row.names = F)

## End



## FIX IDs BEFORE ANALYSIS ----

## Load all resitance distances
dir <- "./5-LandGenAnalysis/results/resist_dist/cha_SS/"
list <- list.files(path=dir, pattern =".csv", full.names=TRUE)
list ## 11 variables

v1 <- read.csv(list.files(path=dir, pattern ="bio10", full.names=TRUE), h=T)
v2 <- read.csv(list.files(path=dir, pattern ="bio15", full.names=TRUE), h=T)
v3 <- read.csv(list.files(path=dir, pattern ="bio3", full.names=TRUE), h=T)
v4 <- read.csv(list.files(path=dir, pattern ="geo_dist", full.names=TRUE), h=T)
v5 <- read.csv(list.files(path=dir, pattern ="elev", full.names=TRUE), h=T)
v6 <- read.csv(list.files(path=dir, pattern ="forest", full.names=TRUE), h=T)
v7 <- read.csv(list.files(path=dir, pattern ="geoenv", full.names=TRUE), h=T)
v8 <- read.csv(list.files(path=dir, pattern ="rough", full.names=TRUE), h=T)
v9 <- read.csv(list.files(path=dir, pattern ="PETs", full.names=TRUE), h=T)
v10 <- read.csv(list.files(path=dir, pattern ="PETdq", full.names=TRUE), h=T)
v11 <- read.csv(list.files(path=dir, pattern ="topoWet", full.names=TRUE), h=T)


## Bind dataframes
final_data <- bind_cols(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11)

head(final_data)

## Select variables and rename
final <- final_data %>% dplyr::select (ID1...1, 
                                       ID2...2,
                                       cd_bio10,
                                       cd_bio15,
                                       cd_bio3,
                                       cd_elevation,
                                       cd_forestcover,
                                       cd_geo_distance,
                                       cd_geoenvironment,
                                       cd_roughness,
                                       cd_PETdq,
                                       cd_PETs,
                                       cd_topoWet) %>%
  dplyr::rename(ID1 = ID1...1,
                ID2 = ID2...2)



## Prepare meta based on original coordinates file used to compute distances
head(coordinates)
meta <- coordinates %>% dplyr::select(sample_name, local_ID) # make new df with sample and cave ids

## Sample IDs
id <- as.character(coordinates$sample_name)
unique(id)

## Replace number IDs by correct sample names
length(unique(final$ID1))

final$ID1 <- plyr::mapvalues(final$ID1, from = sort(unique(final$ID1)), to = id)
final$ID2 <- plyr::mapvalues(final$ID2, from = sort(unique(final$ID2)), to = id)


## Final dataframe: join, filter and rename
df <- final %>% 
  left_join(meta, by=c("ID1" = "sample_name")) %>% # merge dfs by sample name (ID1)
  left_join(meta, by=c("ID2" = "sample_name")) %>% # merge dfs by sample name (ID2)
  dplyr::rename(Cave1=local_ID.x, Cave2=local_ID.y) %>%
  filter(ID1!= ID2) %>%
  unite("ID_1_2", ID1:ID2, sep = "_", remove = FALSE) %>%
  unite("ID_2_1", ID2:ID1, sep = "_", remove = FALSE)


nrow(df) # 8556


## Load relatedness dataframe to use pairs as reference
rel <- read.csv("./5-LandGenAnalysis/results/gen_dist/cha/charinus_relatedness_93samples.csv", sep = ",", h=T)
head(rel)
nrow(rel) ## 4371

## Relatedness: join, filter and rename
gendist.ind <- rel[,1:3] %>% 
  left_join(meta, by=c("ID1" = "sample_name")) %>% # merge dfs by sample name (ID1)
  left_join(meta, by=c("ID2" = "sample_name")) %>% # merge dfs by sample name (ID2)
  filter(ID1!= ID2) %>%
  unite("ID_1_2", ID1:ID2, sep = "_", remove = FALSE) %>%
  unite("ID_2_1", ID2:ID1, sep = "_", remove = FALSE)


head(as_tibble(gendist.ind))
nrow(gendist.ind) ## 4278

## Filtering
## Filter IND pairs in relatedness datraframe
fil1 <- df[as.factor(df$ID_2_1) %in% 
             as.factor(gendist.ind$ID_1_2), ]

nrow(fil1) ## 4278
head(fil1)

## Select and rename
final <- as_tibble(fil1) %>% 
  dplyr:: select(-"ID_1_2") %>%
  dplyr:: rename(ID.ind = ID_2_1)

head(final)
str(final)
nrow(final) ##  4278

## Set diff
setdiff(as.factor(gendist.ind$ID_1_2), as.factor(final$ID.ind)) ## no difference

## Check ids
identical(as.factor(gendist.ind$ID_1_2), as.factor(final$ID.ind)) ## TRUE

## Select variables to save
f <- final %>% dplyr::select(ID1, ID2, Cave1, Cave2, 
                             cd_bio10,
                             cd_bio15,
                             cd_bio3,
                             cd_elevation,
                             cd_forestcover,
                             cd_geo_distance,
                             cd_geoenvironment,
                             cd_roughness,
                             cd_PETdq,
                             cd_PETs,
                             cd_topoWet)
head(f)
str(f)
nrow(f) ##  4278


## Save
write.csv(f, 
          paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/final_ResDist_", project_name, ".csv"), 
          row.names = F)



## END OF SCRIPT ----



