###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## RASTER PREPARATION ----

## Clean Global Environment 
rm(list=ls())

## Load packages
library(raster)
library(rgeos)
library(rgdal)
library(sp)
library(sf)
library(PerformanceAnalytics)

## SETTINGS ----
## Set projection settings before running the script
A=make_EPSG()
LL_prj <- as.character(subset(A, A$code=="4326")[3])  ## decimal degrees
UTM_prj <- as.character(subset(A, A$code=="5383")[3]) ## utm


## Raster resolution settings

## Settings to 1km x 1km rasters
res_1km <- raster(nrow = 34, ncol = 42)
extent(res_1km) <- extent(550701, 585883.2, 9281264, 9309754)
projection(res_1km) <- UTM_prj
res(res_1km) <- c(922, 922)

## Settings for 250m x 250m resolution
res_250m <- raster(nrow = 125, ncol = 170)
extent(res_250m) <- extent(550701, 585883.2, 9281264, 9309754)
projection(res_250m) <- UTM_prj
res(res_250m) <- c(250, 250)
res_250m

## Settings for 30m x 30m resolution
res_30m <- raster(nrow = 7228, ncol = 7234)
extent(res_30m) <- extent(550701, 585883.2, 9281264, 9309754)
projection(res_30m) <- UTM_prj
res(res_30m) <- c(30.7, 30.7)
res_30m

## Set names
project_name = "charinus_SS"


## COORDINATES ----
## Load
coordinates <- read.csv("./coords/cha/coords_charinus_92ind.csv", h=T)
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

## Coordinates and extent in UTM
coord = coord_UTM

## Based on your default projection, choose a buffer size value
v1 = 10000 ## correspond to ~ 10000 meters; use this value when setting a UTM projection
v2 = 0.1 ## correspond to ~ 10000 meters; use this value when setting a Decimal Degrees projection

buff_value = v1

## Set raster extent + buffer
a <- extent(coord)[1] - buff_value
b <- extent(coord)[2] + buff_value
c <- extent(coord)[3] - buff_value
d <- extent(coord)[4] + buff_value
E_UTM <- extent(a,b,c,d)
E_UTM


## Coordinates and extent in deciam degrees LongLat
coord = coord_LongLat
buff_value = v2

## Set raster extent + buffer
a <- extent(coord)[1] - buff_value
b <- extent(coord)[2] + buff_value
c <- extent(coord)[3] - buff_value
d <- extent(coord)[4] + buff_value
E_LongLat <- extent(a,b,c,d)
E_LongLat


## SPATIAL LAYERS ----

## 1. ELEVATION ----
var_name = "elevation"

## Load data
r <- raster("./rasters/raw_rasters/elev/DEM_Carajas_UTM.tif")
crs(r)

## Crop raster
rC <- crop(r, E_UTM)

## Plot cropped raster and check points
plot(rC, axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord_UTM, pch= 19, col="black", cex = 1.2)

## Replace NA if necessary
rC[is.na(rC[])] <- 0.1 # meaning 10x less conductance in those areas
summary(rC) ## No NA's anymore

## Check rasters resolution
res(r) ## 30.7 30.7
res(rC) ## 30.7 30.7
plot(rC)

## Resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_30m ## pre-set raster
Map2 <- resample(rC, res_30m, method='bilinear')
summary(Map2)

## Replace NA if necessary
Map2[is.na(Map2[])] <- 0.1 # meaning 10x less conductance in those areas
summary(Map2) ## No NA's anymore

## Choose map to save:
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "30m"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/30m/", 
            var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)


## If you want to decrease raster resolution, please continue:

## Resample if necessary
## Choose the raster you want to resample
Map = rC

### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_1km ## pre-set raster
Map2 <- resample(Map, res_1km, method='bilinear')
summary(Map2) ## No NA's anymore

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)

## Save map
## Choose map to save:
## reproj_r ## cropped and reprojected
## Map2 ## cropped, reprojected and/or resampled
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "1km"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/1km/", 
           var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
           format = "GTiff", overwrite=T)

## End




## 2. ROUGHNESS ----
var_name = "roughness"

## Load data
r <- raster("./rasters/raw_rasters/elev/Roughness_TRI_Carajas_UTM.tif")
crs(r)

## Crop raster
rC <- crop(r, E_UTM)

## Plot cropped raster and check points
plot(rC, axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord_UTM, pch= 19, col="black", cex = 1.2)

## Check rasters resolution
res(r) ## 30.7 30.7
res(rC) ## 30.7 30.7

## Replace NA if necessary
rC[is.na(rC[])] <- 0.1 # meaning 10x less conductance in those areas
summary(rC) ## No NA's anymore

## Resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_30m ## pre-set raster
Map2 <- resample(rC, res_30m, method='bilinear')
summary(Map2)

## Replace NA if necessary
Map2[is.na(Map2[])] <- 0.1 # meaning 10x less conductance in those areas
summary(Map2) ## No NA's anymore

## Choose map to save:
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "30m"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/30m/", 
            var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)



## If you want to decrease raster resolution, please continue:
## Resample if necessary
## Choose the raster you want to resample
Map = rC

### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_1km ## pre-set raster
Map2 <- resample(Map, res_1km, method='bilinear')

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)

## Save map

## Choose map to save:
## reproj_r ## cropped and reprojected
## Map2 ## cropped, reprojected and/or resampled
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "1km"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/1km/", 
            var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)
## End




## 3. GEOENVIRONMENT ----
## Set names
var_name = "geoenvironment"

## Load data
r <- raster("./rasters/raw_rasters/geoenv/reprojected_geoamb+mina_2019_30m_SS.tif")
crs(r)
plot(r)

## Crop raster
rC <- crop(r, E_UTM)
sort(unique(rC[]))

## Plot cropped raster and check points
plot(rC, axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord_UTM, pch= 19, col="black", cex = 1.2)


## Resample to 30m x 30m resolution
## Choose a map to resample
Map = rC

### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_250m  ## pre-set raster
Map2 <- resample(Map, res_250m, method='ngb')
summary(Map2)

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)

## Replace NA if necessary
unique(Map2[])
summary(Map2)
Map2[is.na(Map2[])] <- 0.1 # meaning 10x less conductance in those areas
summary(Map2) ## No NA's anymore
plot(Map2)

##Choose a map
# reproj_r ## cropped and reprojected
# Map2 ## cropped, reprojected and resampled
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "250m"
projection = "UTM"


## Save
writeRaster(map_to_save, paste0("./5-LandGenAnalysis/results/resistGA/cha_SS/asc/ASC_", 
            var_name, "_", projection, "_", resolution,"_", project_name, ".asc"), 
            format = "ascii", overwrite=T)



## Resample to 1km x 1km

## Choose a map to resample
Map = rC

### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_1km ## pre-set raster
Map2 <- resample(Map, res_1km, method='ngb')

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)

## Replace NA if necessary
unique(Map2[])
summary(Map2)
Map2[is.na(Map2[])] <- 0.1 # meaning 10x less conductance in those areas
summary(Map2) ## No NA's anymore

##Choose a map
# reproj_r ## cropped and reprojected
# Map2 ## cropped, reprojected and resampled
map_to_save = Map2

## Check variable name and save final raster
var_name
resolution = "1km"
projection = "UTM"

writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/1km/", 
            var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)

# End




### 4. FOREST COVER ----

## Set names
var_name = "forestcover"

## Load data
r <- raster("./rasters/raw_rasters/fcover/vegetationcover_2019_30m_SS.tif")
crs(r)

## Crop raster
rC <- crop(r, E_LongLat)
unique(rC[])

## Plot cropped raster and check points
plot(rC, axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord_LongLat, pch= 19, col="black", cex = 1.2)

## Reproject raster if necessary
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
reproj_r <- projectRaster(rC, crs = UTM_prj, method = 'ngb') ## set new projection and reproject
## usei ngb para manter os valores inteiros de cobertura de floresta, pois com metodo 'bilinear' ele faz uma media

crs(reproj_r)
plot(reproj_r)
points(coord_UTM)
summary(reproj_r)

## Replace NA if necessary
reproj_r[is.na(reproj_r[])] <- 0.1 # meaning 10x less conductance in those areas
summary(reproj_r) ## No NA's anymore
unique(reproj_r[])

## Resample if necessary
## Choose the raster you want to resample
Map = reproj_r

## Resample to 1km x 1km resolution
### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_1km
Map2 <- resample(Map, res_1km, method='ngb')
summary(Map2)

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)


## Save map

## Choose map to save:
# reproj_r ## cropped and reprojected
# Map2 ## cropped, reprojected and/or resampled
map_to_save = Map2

## Check var_name and save final raster
var_name
resolution = "1km"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/1km/", var_name, "_", 
            projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)


## Resample to 30m x 30m resolution

## Choose the raster you want to resample
Map = reproj_r

### Decrease raster resolution using resample
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
res_30m
Map2 <- resample(Map, res_30m, method='ngb')
summary(Map2)

## Plot resampled raster
plot(Map2, axes=T)
points(coord_UTM)
summary(Map2)

## Saved map

## Choose map to save:
# reproj_r ## cropped and reprojected
# Map2 ## cropped, reprojected and/or resampled
map_to_save = Map2


## Check var_name and save final raster
var_name
resolution = "30m"
projection = "UTM"
writeRaster(map_to_save, paste0("./rasters/crop_rasters/cha_SS/30m/", var_name, "_", 
            projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff", overwrite=T)




### 5. BIOCLIMATICS ----

## Set names
var_name = "bioclim"

## Load data
dir <- "./rasters/raw_rasters/worldclim2/"
list <- list.files(path=dir, pattern =".tif", full.names=TRUE)
list ## 19 variables

## Use raster stack
bio_stack <- stack(list)
names(bio_stack)
projection(bio_stack) <- LL_prj


## Use PCA to selected variables
## Extract variables from rasters
names(bio_stack)
EVars <- raster::extract(bio_stack, coord_LongLat)
summary(EVars)

## Matrix to dataframe 
class(EVars)
EVars <- as.data.frame(EVars)

## Rename rows
head(coordinates)
rownames(EVars) <- coordinates[ ,"ID"]
head(EVars)

## Check dataframe
str(EVars)
head(EVars)
ncol(EVars) #19

### PCA 19 bioclimatics
### Run PCA
PC <- princomp(scale(EVars))
summary(PC)
head(PC) 
loadings(PC)
plot(PC)
biplot(PC)

##Correlate EVars with PC
var1 <- which.max(abs(cor(EVars, PC$scores[, 1])))
var2 <- which.max(abs(cor(EVars, PC$scores[, 2]))) 
var3 <- which.max(abs(cor(EVars, PC$scores[, 3]))) 
var4 <- which.max(abs(cor(EVars, PC$scores[, 4]))) 

sel <- c(var1, var2, var3, var4)

##Select most informative EVars
names(EVars)
Vars <- EVars[, sel]
head(Vars)
str(Vars)
names(Vars) ## "wc2.1_30s_bio_10" "wc2.1_30s_bio_15" "wc2.1_30s_bio_14" "wc2.1_30s_bio_3"

## Chart correlation
chart.Correlation(Vars, histogram=TRUE, pch=19) ##bio18 highly correlated

## Filter rasters selected
selected_rasters <- bio_stack[[sel[c(1:2,4)]]]
class(selected_rasters)
names(selected_rasters)


## Crop raster stack
rC <- crop(selected_rasters, E_LongLat)

## Plot cropped raster and check points
plot(rC[[1]], axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord, pch= 19, col="black", cex = 1.2)

## Reproject raster if necessary
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
reproj_r <- projectRaster(rC, crs = UTM_prj, method = 'bilinear', bylayer=TRUE, progress='text', snap="out") ## set new projection and reproject
crs(reproj_r)
plot(reproj_r[[1]])
points(coord_UTM)
summary(reproj_r)

## Resample to the same extent
resamp_r <- resample(reproj_r, res_1km, method = 'bilinear', bylayer=TRUE, progress='text', snap="out")
summary(resamp_r)

## Choose map to save:
map_to_save = resamp_r

## Save maps
var_name
resolution = "1km"
projection = "UTM"

writeRaster(map_to_save, 
            paste0("./rasters/crop_rasters/cha_SS/1km/", 
                   var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff",  
            bylayer = T, suffix='names', overwrite=T)


## End


## 6. ENVIREM ----

## Set names
var_name = "envirem"

## Load data
dir <- "./rasters/raw_rasters/envirem2/"
list <- list.files(path=dir, pattern =".bil", full.names=TRUE)
list ## 12 variables

## Use raster stack
env_stack <- stack(list)
names(env_stack)
projection(env_stack) <- LL_prj

## Use PCA to selected variables
## Extract variables from rasters
names(bio_stack)
EVars <- raster::extract(env_stack, coord_LongLat)
summary(EVars)

## Matrix to dataframe 
class(EVars)
EVars <- as.data.frame(EVars)

## Rename rows
rownames(EVars) <- coordinates[ ,"ID"]
head(EVars)

## Check dataframe
str(EVars)
head(EVars)
ncol(EVars) # 12

### PCA 12 Envirem
### Run PCA
PC <- princomp(scale(EVars))
summary(PC)
head(PC) 
loadings(PC)
plot(PC)
biplot(PC)

##Correlate EVars with PC
var1 <- which.max(abs(cor(EVars, PC$scores[, 1])))
var2 <- which.max(abs(cor(EVars, PC$scores[, 2]))) 
var3 <- which.max(abs(cor(EVars, PC$scores[, 3]))) 

sel <- c(var1, var2, var3)

##Select most informative EVars
names(EVars)
Vars <- EVars[, sel]
head(Vars)
str(Vars)
names(Vars) ## "current_30arcsec_PETDriestQuarter" "current_30arcsec_topoWet" "current_30arcsec_PETseasonality"

## Chart correlation
chart.Correlation(Vars, histogram=TRUE, pch=19)

## Filter rasters selected
selected_rasters <- env_stack[[sel]]
class(selected_rasters)
names(selected_rasters)

## Crop raster stack
rC <- crop(selected_rasters, E_LongLat)

## Plot cropped raster and check points
plot(rC[[1]], axes = T, xlab="Longitude", ylab="Latitude", main = paste0(var_name))
points(coord, pch= 19, col="black", cex = 1.2)

## Reproject raster if necessary
## method = 'bilinear' if working with continuous data
## method = 'ngb' if working with classes/categories
reproj_r <- projectRaster(rC, crs = UTM_prj, method = 'bilinear', bylayer=TRUE, progress='text', snap="out") ## set new projection and reproject
crs(reproj_r)
plot(reproj_r[[1]])
points(coord_UTM)

## Replace NA if necessary
reproj_r[is.na(reproj_r[])] <- 0.1 # meaning 10x less conductance in those areas
summary(reproj_r) ## No NA's anymore

## Resample to the same extent
resamp_r <- resample(reproj_r, res_1km, method = 'bilinear', bylayer=TRUE, progress='text', snap="out")
resamp_r
summary(resamp_r)

## Choose map to save:
map_to_save = resamp_r

## Save maps
var_name
resolution = "1km"
projection = "UTM"

writeRaster(map_to_save, 
            paste0("./rasters/crop_rasters/cha_SS/1km/", 
                   var_name, "_", projection, "_", resolution,"_", project_name, ".tif"), 
            format = "GTiff",  
            bylayer = T, suffix='names', overwrite=T)




## 7. CHECK ALL SAVED RASTERS ----

## 1km
dir <- "./rasters/crop_rasters/cha_SS/1km/"
final.list <- list.files(path=dir, pattern =".tif", full.names=TRUE)
length(final.list) ## 10

final_stack <- stack(final.list) ## ok
names(final_stack)

## 30m
dir <- "./rasters/crop_rasters/cha_SS/30m/"
final.list <- list.files(path=dir, pattern =".tif", full.names=TRUE)
length(final.list) ## 3

final_stack <- stack(final.list) ## ok
names(final_stack)


## END OF SCRIPT ----

### Bioclimate variables
#Alt = Altitude  
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter
