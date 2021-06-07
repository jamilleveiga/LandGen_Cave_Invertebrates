###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

### Topographic Distance ----

## Clean Global Environment 
rm(list=ls())

## Load libraries
#devtools::install_github("ianjwang/topoDistance", force = T)
library(topoDistance)
library(sp)
library(raster)
library(rgdal)
library(geosphere)
library(tidyverse)
library(reshape2)

## Set projection settings before running the script ####
A=make_EPSG()
LL_prj <- as.character(subset(A, A$code=="4326")[3]) 
UTM_prj <- as.character(subset(A, A$code=="5383")[3]) 


## Set names
project_name = "charinus_SS"

## 1. COORDINATES ----
## Load
coordinates <- read.csv("./coords/cha/coords_charinus_93ind.csv", h=T)
head(coordinates)

## Set coords projection
coord_UTM <- coordinates[, 3:4]
coordinates(coord_UTM) <- coord_UTM
projection(coord_UTM) <- crs(UTM_prj)
plot(coord_UTM, axes=T)

## Reproject coordinates if necessary
coord_LongLat <- spTransform(coord_UTM, CRSobj = CRS(LL_prj))
crs(coord_LongLat)
plot(coord_LongLat, axes=T)


## 2. CALCULATE GEOGRAPHIC DISTANCE ----
eucl_dist = distm(coord_LongLat, fun=distGeo)
eucl_dist

## Converting in KM
eucl_dist = eucl_dist/1000

## Convert distance matrix into data frame for analyses
eucl_dist = eucl_dist %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", "geoDist")) %>%
  filter(ID1!= ID2)

## Check
head(eucl_dist)
summary(eucl_dist)
hist(eucl_dist[,3])
nrow(eucl_dist) ## 8556
min(eucl_dist$geoDist)
max(eucl_dist$geoDist) ## 15.86804

## 3. LOAD RASTER ----
var_name = "elevation"

## Load raster
r <- raster("./rasters/crop_rasters/cha_SS/30m/elevation_UTM_30m_charinus_SS.tif")
crs(r)
plot(r)

## 4. CALCULATE TOPOGRAPHIC DISTANCE ----

start_time <- Sys.time()
#topo8 <- topoDist(r, coord_UTM, directions = 8, paths = T, zweight = 1)

## Save RData  
save(topo8, file = "./5-LandGenAnalysis/RData/topo_dist_cha_SS.RData")
end_time <- Sys.time()
end_time - start_time  # 2.907572 hours

load("./5-LandGenAnalysis/RData/topo_dist_cha_SS.RData")


# 3. CONVERT DISTANCE MATRIX TO DATAFRAME

## Convert distance matrix into data frame for analyses
topdist8 = topo8[[1]] %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", "topodist"))

topdist8 = topdist8  %>% 
           filter(ID1!=ID2)

## Check
head(topdist8)
hist(topdist8[,3])
nrow(topdist8) ## 8372
topdist8$topodist <- topdist8$topodist/1000
min(topdist8$topodist) ## 0
max(topdist8$topodist) ## 17.19169


## 5. FIX IDs BEFORE ANALYSIS ----

## Rename
final = bind_cols(topdist8, eucl_dist) %>%
        dplyr::select(ID1...1, ID2...2, topodist, geoDist) %>%
        dplyr::rename(ID1 = ID1...1,
                      ID2 = ID2...2)
head(final)

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
  dplyr::rename(Cave1 = local_ID.x, Cave2 = local_ID.y) %>%
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
final <- as_tibble(fil1) %>% dplyr:: select(-"ID_1_2") %>%
  dplyr:: rename (ID.ind = ID_2_1)

head(final)

## Set diff
setdiff(as.factor(gendist.ind$ID_1_2), as.factor(final$ID.ind)) ## no difference

## Check ids
identical(as.factor(gendist.ind$ID_1_2), as.factor(final$ID.ind)) ## TRUE

## Select variables to save
head(final)
f <- final %>% dplyr::select(ID1, ID2, Cave1, Cave2, topodist, geoDist)
str(f)
names(f)


## Save
write.csv(f, 
          paste0("./5-LandGenAnalysis/results/topo_dist/cha_SS/topdist_8_", project_name, ".csv"), row.names = F)

## Save
write.csv(f, 
          paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/final_TopDist_", project_name, ".csv"), 
          row.names = F)


## END OF SCRIPT ----
