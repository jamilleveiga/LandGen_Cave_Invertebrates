###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

### RECLASSIFY RASTERS ----

## Clean Global Environment 
rm(list=ls())

## Load packages
library(sp)
library(raster)
library(rgdal)
library(geosphere)

## Set projection settings before running the script 
A=make_EPSG()
LL_prj <- as.character(subset(A, A$code=="4326")[3])  ## decimal degrees
UTM_prj <- as.character(subset(A, A$code=="5383")[3]) ## utm

## Set names
project_name = "charinus_SS"


## COORDINATES ----
## Load
coordinates <- read.csv("./coords/cha/coords_charinus_93ind.csv", h=T)
head(coordinates)

## Set coords projection
coord_UTM <- coordinates[, 3:4]
coordinates(coord_UTM) <- coord_UTM
projection(coord_UTM) <- crs(UTM_prj)
plot(coord_UTM, axes=T)


### BUILD RESISTANCE SURFACES ----

### 0. GEOGRAPHIC DISTANCE ----
### Use any raster to build a geographic distance surface
var_name = "geo_distance"

r <- raster("./rasters/crop_rasters/cha_SS/1km/elevation_UTM_1km_charinus_SS.tif")
projection(r) <- UTM_prj

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Substituir todos os valores por 0.5
r[!is.na(r)] <- 0.5

## Check values
unique(r[])
plot(r)
points(coord_UTM)
dim(r) # 34 46  1
res(r) # 922 922

## Rename raster
names(r) <- "geographic_distance"

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
            var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)



### 1. ELEVATION ----
var_name = "elevation"

r <- raster("./rasters/crop_rasters/cha_SS/1km/elevation_UTM_1km_charinus_SS.tif")
projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
dim(r) # 34 46  1     

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)

## End



### 2. ROUGHNESS ----
var_name = "roughness"

r <- raster("./rasters/crop_rasters/cha_SS/1km/roughness_UTM_1km_charinus_SS.tif")
projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
dim(r) # 34 46  1 
res(r) # 922 922

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)


## 3. GEOENVIRONMENT ----
## We used resistanceGA to optimize geoenvironment surface
var_name = "geoenvironment"

## Load data
r <- raster("./5-LandGenAnalysis/results/resistGA/cha_SS/mean_optimized_resist_values_geoenvrionment.tif")
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
sort(unique(r[]))
res(r) # 922 922
dim(r) # 114 141   1    

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)


### 4. FOREST COVER ----

## Set names
var_name = "forestcover"

## Load data
r <- raster("./rasters/crop_rasters/cha_SS/1km/forestcover_UTM_1km_charinus_SS.tif")
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Higher forest cover values will decrease resistance to gene flow
r[] <- abs(r[] - max(r[]))
unique(r[])
summary(r[])

## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 34 42  1     

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)



### 5. BIOCLIMATICS ----

## Load data
dir <- "./rasters/crop_rasters/cha_SS/1km/"
list <- list.files(path=dir, pattern ="bioclim_", full.names=TRUE)
list ## 3 variables


### BIO10 ----
list[1] ## bio10
#BIO10 = Mean Temperature of Warmest Quarter

## Set name
var_name = "bio10"

## Load data
r <- raster(list[1])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 31 38  1     

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)

## End

### BIO15 ----
list[2] ## bio15
#BIO15 = Precipitation Seasonality (Coefficient of Variation)

## Set name
var_name = "bio15"

## Load data
r <- raster(list[2])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 31 38  1   

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)


## End

### BIO3 ----
list[3] ## bio3
#BIO3 = Isothermality (BIO2/BIO7) (* 100)

## Set name
var_name = "bio3"

## Load data
r <- raster(list[3])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 31 38  1     

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)


### 6. ENVIREM ----

## Load data
dir <- "./rasters/crop_rasters/cha_SS/1km/"
list <- list.files(path=dir, pattern ="envirem_", full.names=TRUE)
list ## 3 variables


### PETDriestQuarter ----
list[1] ## PETdq

## Set name
var_name = "PETdq"

## Load data
r <- raster(list[1])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 31 38  1

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)

## End

### PETseasonality ----
list[2] ## PETseasonality

## Set name
var_name = "PETs"

## Load data
r <- raster(list[2])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 31 38  1     

## Check name and save
var_name
map_to_save = r
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)

## End


### TopoWet ----
list[3] ## topoWet

## Set name
var_name = "topoWet"

## Load data
r <- raster(list[3])
crs(r)

projection(r) <- UTM_prj
hist(r[])
summary(r[])

###Plot
plot(r)
points(coord_UTM)
unique(r[])

### Save values according to predicted effects (hypothesis table)
## Low topographic wetness will decrease resistance to gene flow
r[] <- abs(r[] - max(r[]))
unique(r[])
summary(r[])

## Check values
unique(r[])
plot(r)
points(coord_UTM)
res(r) # 922 922
dim(r) # 34 46  1    

## Check name and save
map_to_save = r

var_name
writeRaster(map_to_save, paste0("./rasters/resist_surfaces/cha_SS/", 
                                var_name, "_", project_name, ".tif"), format = "GTiff", overwrite=T)



## END OF SCRIPT ----

