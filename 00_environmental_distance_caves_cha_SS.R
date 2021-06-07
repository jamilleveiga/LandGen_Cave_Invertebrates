###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


## LANDSCAPE METRICS ----

## Clean Global Environment 
rm(list=ls())

## Load packages
library(tidyverse)
library(reshape2)
library(dplyr)
library(ecodist)

## Load source
source("./source/source_ecodist.R")

## Set project name
project_name = "charinus_SS"

## LOAD LANDSCAPE METRICS AND CAVE METADATA ----
land_metrics <- read.table("./5-LandGenAnalysis/landscapemetrics/landscape_metrics.txt")
head(land_metrics)

meta_caves <- read.csv("./5-LandGenAnalysis/landscapemetrics/Cavedata_Charinus_FINAL.csv", sep=";", h=T)
head(meta_caves)


## 1. COORDINATES ----
## Load cave points
coords_cave <- read.csv("./coords/cha/coords_charinus_27ind_OneSampleByCave.csv", h=T)
head(coords_cave)
nrow(coords_cave) ## 27
coords_cave$local_ID

## Remove duplicated
pontos_unicos = coords_cave[!duplicated(coords_cave[c("longitude","latitude")]), ]

## Numero de pontos unicos:
length(pontos_unicos[, 1]) #27


## 2. FILTER LANDSCAPE METRICS ----
fil1_land_metrics <- land_metrics[land_metrics$Caverna %in% pontos_unicos$local_ID, ]
head(fil1_land_metrics)
nrow(fil1_land_metrics) ## 27

## Check order of IDs
identical(as.character(pontos_unicos$local_ID), as.character(fil1_land_metrics$Caverna))  ## False

## Arrange
fil1_land_metrics <- fil1_land_metrics %>%
                     mutate(Caves =  factor(Caverna, levels = as.factor(pontos_unicos$local_ID)))%>%
                     arrange(Caves); head(fil1_land_metrics)

## Check order of IDs
identical(as.character(pontos_unicos$local_ID), as.character(fil1_land_metrics$Caves)) ## True 



## 3. FILTER CAVE METADATA ----
fil3_meta_caves <- meta_caves[meta_caves$Cavidade %in% pontos_unicos$local_ID, ]
head(fil3_meta_caves)
nrow(fil3_meta_caves) ## 27

## Check order of IDs
identical(as.character(pontos_unicos$local_ID), as.character(fil3_meta_caves$Cavidade))  ## False

## Arrange
fil3_meta_caves <- fil3_meta_caves %>%
                   mutate(Caves =  factor(Cavidade, levels = as.factor(pontos_unicos$local_ID)))%>%
                   arrange(Caves); head(fil3_meta_caves)

## Check order of IDs
identical(as.character(pontos_unicos$local_ID), as.character(fil3_meta_caves$Caves)) ## True 


## 4. CALCULATE ENVIRONMENTAL DISTANCES ----

## FOREST B500 2019 ----
var_name = "B500_Forest_pland_2019"

env_d_B500_Forest <- as.matrix(dist(fil1_land_metrics$B500_Forest_pland_2019, method = "euclidian"))
env_d_B500_Forest[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B500_Forest = env_d_B500_Forest %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B500_Forest = env_d_B500_Forest  %>% 
  filter(ID1!=ID2)

class(env_d_B500_Forest)   
head(env_d_B500_Forest)
hist(env_d_B500_Forest[,3])
nrow(env_d_B500_Forest) ## 702


##  MINING B500 2019 ----
var_name = "B500_Minning_pland_2019"

env_d_B500_Mining <- as.matrix(dist(fil1_land_metrics$B500_Minning_pland_2019, method = "euclidian"))
env_d_B500_Mining[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B500_Mining = env_d_B500_Mining %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B500_Mining = env_d_B500_Mining  %>% 
  filter(ID1!=ID2)

class(env_d_B500_Mining)   
head(env_d_B500_Mining)
hist(env_d_B500_Mining[,3])
nrow(env_d_B500_Mining) ## 702




##  CANGA B500 2019 ----
var_name = "B500_Canga_pland_2019"

env_d_B500_Canga <- as.matrix(dist(fil1_land_metrics$B500_Canga_pland_2019, method = "euclidian"))
env_d_B500_Canga[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B500_Canga = env_d_B500_Canga %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B500_Canga = env_d_B500_Canga  %>% 
  filter(ID1!=ID2)

class(env_d_B500_Canga)   
head(env_d_B500_Canga)
hist(env_d_B500_Canga[,3])
nrow(env_d_B500_Canga) ## 702




## ALL VEGETATION B500 2019 ----
var_name = "B500_Allvegetation_pland_2019"

env_d_B500_Allvegetation <- as.matrix(dist(fil1_land_metrics$B500_Allvegetation_pland_2019, method = "euclidian"))
env_d_B500_Allvegetation[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B500_Allvegetation = env_d_B500_Allvegetation %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B500_Allvegetation = env_d_B500_Allvegetation  %>% 
  filter(ID1!=ID2)

class(env_d_B500_Allvegetation)   
head(env_d_B500_Allvegetation)
hist(env_d_B500_Allvegetation[,3])
nrow(env_d_B500_Allvegetation) ## 702




## FOREST B1000 2019 ----
var_name = "B1000_Forest_pland_2019"

env_d_B1000_Forest <- as.matrix(dist(fil1_land_metrics$B1000_Forest_pland_2019, method = "euclidian"))
env_d_B1000_Forest[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B1000_Forest = env_d_B1000_Forest %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B1000_Forest = env_d_B1000_Forest  %>% 
  filter(ID1!=ID2)

class(env_d_B1000_Forest) 
head(env_d_B1000_Forest)
hist(env_d_B1000_Forest[,3])
nrow(env_d_B1000_Forest) ## 702




##  MINING B1000 2019 ----
var_name = "B1000_Minning_pland_2019"

env_d_B1000_Mining <- as.matrix(dist(fil1_land_metrics$B1000_Minning_pland_2019, method = "euclidian"))
env_d_B1000_Mining[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B1000_Mining = env_d_B1000_Mining %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B1000_Mining = env_d_B1000_Mining  %>% 
  filter(ID1!=ID2)

class(env_d_B1000_Mining) 
head(env_d_B1000_Mining)
hist(env_d_B1000_Mining[,3])
nrow(env_d_B1000_Mining) ## 702


##  CANGA B1000 2019 ----
var_name = "B1000_Canga_pland_2019"

env_d_B1000_Canga <- as.matrix(dist(fil1_land_metrics$B1000_Canga_pland_2019, method = "euclidian"))
env_d_B1000_Canga[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B1000_Canga = env_d_B1000_Canga %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B1000_Canga = env_d_B1000_Canga  %>% 
  filter(ID1!=ID2)

class(env_d_B1000_Canga) 
head(env_d_B1000_Canga)
hist(env_d_B1000_Canga[,3])
nrow(env_d_B1000_Canga) ## 702


## ALL VEGETATION B1000 2019 ----
var_name = "B1000_Allvegetation_pland_2019"

env_d_B1000_Allvegetation <- as.matrix(dist(fil1_land_metrics$B1000_Allvegetation_pland_2019, method = "euclidian"))
env_d_B1000_Allvegetation[1:10,1:10]

## Convert distance matrix into data frame for analyses
env_d_B1000_Allvegetation = env_d_B1000_Allvegetation %>%
  melt %>%
  na.omit %>%
  setNames(c("ID1", "ID2", paste0("d_", var_name)))

env_d_B1000_Allvegetation = env_d_B1000_Allvegetation  %>% 
  filter(ID1!=ID2)

class(env_d_B1000_Allvegetation)   
head(env_d_B1000_Allvegetation)
hist(env_d_B1000_Allvegetation[,3])
nrow(env_d_B1000_Allvegetation) ## 702


## CAVE HORIZONTAL PROJECTION ----
## Soma da projecao horizontal entre cavernas
var_name = "hp"

env_d_hp = pairedsum(as.matrix(fil3_meta_caves$PH))
class(env_d_hp)
dim(env_d_hp)

#Converting distance matrix into data frame for analyses:
env_d_hp = env_d_hp %>%
  melt %>%
  na.omit %>%
  #arrange(., Var1) %>%
  setNames(c("ID1", "ID2", "Var", paste0("sum_", var_name)))

env_d_hp = env_d_hp  %>% 
  select(-Var) %>%
  filter(ID1!=ID2)

class(env_d_hp)   
head(env_d_hp)
summary(env_d_hp)
hist(env_d_hp[,3])
nrow(env_d_hp) ## 702


## 5. BIND DISTANCES, APPEND CAVE IDs AND SAVE ----

## Bind dataframes
final_data <- bind_cols(env_d_hp,
                        env_d_B1000_Allvegetation, 
                        env_d_B1000_Canga, env_d_B1000_Forest, env_d_B1000_Mining,
                        env_d_B500_Allvegetation,
                        env_d_B500_Canga, env_d_B500_Forest, env_d_B500_Mining)

head(final_data)



## Select variables and rename
names(final_data)
final <- final_data %>% dplyr::select (ID1...1, 
                                       ID2...2,
                                       sum_hp,
                                       d_B1000_Allvegetation_pland_2019,
                                       d_B1000_Canga_pland_2019,
                                       d_B1000_Forest_pland_2019,
                                       d_B1000_Minning_pland_2019,
                                       d_B500_Allvegetation_pland_2019,
                                       d_B500_Canga_pland_2019,
                                       d_B500_Forest_pland_2019,
                                       d_B500_Minning_pland_2019) %>%
  dplyr::rename(Cave1 = ID1...1,
                Cave2 = ID2...2,
                d_B1000_Mining2019 = d_B1000_Minning_pland_2019,
                d_B500_Mining2019 = d_B500_Minning_pland_2019)




str(final)


## Prepare meta based on original coordinates file used to compute distances
head(coords_cave)
meta <- pontos_unicos %>% dplyr::select(sample_name, local_ID) # make new df with sample and cave ids

## Cave IDs
id <- as.character(meta$local_ID)
unique(id)

## Replace number IDs by correct sample names
final$Cave1 <- plyr::mapvalues(final$Cave1, from = sort(unique(final$Cave1)), to = id)
final$Cave2 <- plyr::mapvalues(final$Cave2, from = sort(unique(final$Cave2)), to = id)

final <- final %>% filter(Cave1!= Cave2)

head(final)
nrow(final) ## 702

## Save
project_name
write.csv(final, 
          paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/Caves_Metrics_Dist_", project_name, ".csv"), 
          row.names = F)



## END OF SCRIPT ----


