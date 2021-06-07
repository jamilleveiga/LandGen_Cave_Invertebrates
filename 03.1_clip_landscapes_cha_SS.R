###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## CLIPPING LANDSCAPES USING AN ELLIPSE ----

## Clean Global Environment 
rm(list=ls())

## Load packages
library(purrr)      # mapping
library(suppoRt)   # devtools::install_github("mhesselbarth/helpeR")
library(maptools)   # read shape files
library(PBSmapping) # make polygons
library(raster)     # read raster files and clipping
library(spatstat)   # for random walker simulation
library(rgdal)      # SIG settings
library(rasterVis)  # Visualize rasters in 3D

## Load source
source("source/source_gsm_geodiv.R")

## Set projection settings before running the script
A=make_EPSG()
LL_prj <- as.character(subset(A, A$code=="4326")[3])  ## decimal degrees
UTM_prj <- as.character(subset(A, A$code=="5383")[3]) ## utm

## Set project name
project_name = "charinus_SS"


## INPUT DATA FOR CLIPPING LANDSCAPES ----

## COORDINATES ----

## Load cave points
coords_cave <- read.csv("./coords/cha/coords_charinus_27ind_OneSampleByCave.csv", h=T)
head(coords_cave)
nrow(coords_cave) ## 27
unique(coords_cave$local_ID)

## Set cave coords projection
sp <- coords_cave[, 3:4]
coordinates(sp) <- sp
projection(sp) <- crs(UTM_prj)
plot(sp, axes=T)

## Setting numeric id of sampling points
id <- as.numeric(rownames(sp@data))
id

## Create distance matrix 
dist.m <- as.matrix(dist(sp@coords, method = "euclidean", 
                         diag = F, upper = T))

to_disk <- FALSE
rasterOptions(todisk = to_disk)


## Load layers
dir <- "./rasters/crop_rasters/cha_SS/30m/"
final.list <- list.files(path=dir, pattern =".tif", full.names=TRUE)
length(final.list) ## 3


## 1. ELEVATION ----
## Set var name
var_name = "elevation"

## Check list and and choose position to load raster
final.list
r <- final.list[1]

## Load raster
input_layer <-  raster(r)
crs(input_layer)
hist(input_layer)

## Clip ellipsoide landscapes
## Caution: when building ellipsoid landscapes, NAs are not allowed!!! So here we considered only one sample by cave.
## Loop for all sampling points

clipped <- clip_landscapes(id, sp, input_layer)
clipped

clippings_flatten <- purrr::flatten(clipped)
length(clippings_flatten) ## 351


## Check plots
plot(clippings_flatten[[10]])
plot(clippings_flatten[[50]])
points(sp)

## Check name and save
var_name
suppoRt::save_rds(clippings_flatten, 
                  filename =  paste0("clippings_flatten_", var_name, ".rds"), 
                  path =  "./rasters/clipped_landscapes/cha_SS/", 
                  overwrite = T)

## End


## 2. FORESTCOVER ----
## Set var name
var_name = "forestcover"

## Check list and and choose position to load raster
final.list
r <- final.list[2]

## Load raster
input_layer <-  raster(r)
crs(input_layer)
hist(input_layer)

## Clip ellipsoide landscapes
## Caution: when building ellipsoid landscapes, NAs are not allowed!!! So here we considered only one sample by cave.
## Loop for all sampling points
clipped <- clip_landscapes(id, sp, input_layer)
clipped

clippings_flatten <- purrr::flatten(clipped)
length(clippings_flatten) ## 351

plot(clippings_flatten[[17]])
plot(clippings_flatten[[50]])
points(sp)

## Check name and save
var_name
suppoRt::save_rds(clippings_flatten, 
                  filename =  paste0("clippings_flatten_", var_name, ".rds"), 
                  path =  "./rasters/clipped_landscapes/cha_SS/", 
                  overwrite = T)



## 3. ROUGHNESS ----
## Set var name
var_name = "roughness"

## Check list and and choose position to load raster
final.list
r <- final.list[3]

## Load raster
input_layer <-  raster(r)
crs(input_layer)
hist(input_layer)

## Clip ellipsoide landscapes
## Caution: when building ellipsoid landscapes, NAs are not allowed!!! So here we considered only one sample by cave.
## Loop for all sampling points
clipped <- clip_landscapes(id, sp, input_layer)
clipped

clippings_flatten <- purrr::flatten(clipped)
length(clippings_flatten) ## 351


## Check plots
plot(clippings_flatten[[10]])
plot(clippings_flatten[[50]])
points(sp)

## Check name and save
var_name
suppoRt::save_rds(clippings_flatten, 
                  filename =  paste0("clippings_flatten_", var_name, ".rds"), 
                  path =  "./rasters/clipped_landscapes/cha_SS/", 
                  overwrite = T)



### Plot 3D ----

## Load clipped rasters
clippings1 <- readRDS("./rasters/clipped_landscapes/pha_SS/clippings_flatten_elevation.rds")
clippings2 <- readRDS("./rasters/clipped_landscapes/pha_SS/clippings_flatten_roughness.rds")

## Check plots
plot(clippings1[[17]])
plot(clippings2[[17]])

## Visualize ellipses in 3D
elip1 <- clippings1[[17]]
plot(elip1)

elip2 <- clippings2[[17]]
plot(elip2)

## 3D visualization
rasterVis::plot3D(elip1, 
                  adjust=T, 
                  useLegend = T,
                  zfac=1,
                  rev=T)

rgl::open3d()

rasterVis::plot3D(elip2, 
                  adjust=T, 
                  useLegend = T,
                  zfac=1,
                  rev=T)


## END OF SCRIPT ----
