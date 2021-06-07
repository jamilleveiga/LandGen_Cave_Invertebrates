###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## Calculate surface metrics

## Clean Global Environment 
rm(list=ls())

# load libraries
library(tidyverse)
library(geodiv)
library(plyr)

## Load source
source("source/source_gsm_geodiv.R")


## Set names
project_name = "charinus_SS"

# load sampling points
## Load cave points
coords_cave <- read.csv("./coords/cha/coords_charinus_27ind_OneSampleByCave.csv", h=T)
head(coords_cave)
nrow(coords_cave) ## 27
unique(coords_cave$local_ID)

## Cave IDs
cave_id <- as.character(coords_cave$local_ID) # setting unique id of sampling points
unique(cave_id)

## Split cave IDs in two sets
cave_id1 <- cave_id[1:26]
cave_id2 <- cave_id[2:27]

#### Load data ####

## 1. ELEVATION ----
var_name = "elevation"

# # load the clippings
clippings <- readRDS("./rasters/clipped_landscapes/cha_SS/clippings_flatten_elevation.rds")


# Calculate metrics locally but overall printing progress
gsm_calc <- calculate_gsm(clippings)

## Save RData
var_name
save(gsm_calc, file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))
load(file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))


# Bind to one dataframe
gsm_data <- dplyr::bind_rows(gsm_calc)
head(gsm_data)
nrow(gsm_data)

## Prepare full dataframe
gsm_new <- gsm_data %>%
            separate(sample, c("A", "Cave1", "Cave2")) %>% 
            dplyr::select(., -"A")

head(gsm_new)


## Check length -- should match
length(unique(gsm_new$Cave1))
length(unique(gsm_new$Cave2))

length(cave_id1)
length(cave_id2)

## Append cave ID
head(gsm_new)
gsm_new$caveID1 <- mapvalues(gsm_new$Cave1, from = unique(gsm_new$Cave1), to = cave_id1)
gsm_new$caveID2 <- mapvalues(gsm_new$Cave2, from = unique(gsm_new$Cave2), to = cave_id2)


## Check, select and rename
head(gsm_new)
gsm_final <- gsm_new %>% 
  dplyr::select(caveID1, caveID2, sa, s10, ssk, sdr, srwi, sbi) %>% 
  dplyr::rename(., sa_elev = sa, 
                s10z_elev = s10, 
                ssk_elev = ssk,
                sdr_elev = sdr, 
                srwi_elev = srwi, 
                sbi_elev = sbi)



gsm_final
head(gsm_final)
nrow(gsm_final) ## 351

## Save data
var_name
write.csv(gsm_final, paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_data_", var_name, "_", project_name, ".csv"), row.names = F)

## End


## 2. FORESTCOVER ----    
## Set var name
var_name = "forestcover"

# # load the clippings
clippings <- readRDS("./rasters/clipped_landscapes/cha_SS/clippings_flatten_forestcover.rds")

# Calculate metrics locally but overall printing progress
gsm_calc <- calculate_gsm(clippings)


## Save RData
var_name
save(gsm_calc, file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))
load(file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))

# Bind to one dataframe
gsm_data <- dplyr::bind_rows(gsm_calc)
head(gsm_data)
nrow(gsm_data)

## Prepare full dataframe
gsm_new <- gsm_data %>%
  separate(sample, c("A", "Cave1", "Cave2")) %>% 
  dplyr::select(., -"A")

head(gsm_new)


## Check length -- should match
length(unique(gsm_new$Cave1))
length(unique(gsm_new$Cave2))

length(cave_id1)
length(cave_id2)

## Append cave IDs
head(gsm_new)
gsm_new$caveID1 <- mapvalues(gsm_new$Cave1, from = unique(gsm_new$Cave1), to = cave_id1)
gsm_new$caveID2 <- mapvalues(gsm_new$Cave2, from = unique(gsm_new$Cave2), to = cave_id2)


## Check, select and rename
head(gsm_new)
gsm_final <- gsm_new %>% 
  dplyr::select(caveID1, caveID2, sa, s10, ssk, sdr, srwi, sbi) %>% 
  dplyr::rename(., sa_fcov = sa, 
                s10z_fcov = s10, 
                ssk_fcov = ssk,
                sdr_fcov = sdr, 
                srwi_fcov = srwi, 
                sbi_fcov = sbi)


gsm_final
head(gsm_final)
nrow(gsm_final) ## 351

## Save data
var_name
write.csv(gsm_final, paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_data_", var_name, "_", project_name, ".csv"), row.names = F)

## End



## 3. ROUGHNESS ----
## Set var name
var_name = "roughness"

# # load the clippings
clippings <- readRDS("./rasters/clipped_landscapes/cha_SS/clippings_flatten_roughness.rds")

# Calculate metrics locally but overall printing progress
gsm_calc <- calculate_gsm(clippings)


## Save RData
var_name
save(gsm_calc, file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))
load(file = paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_", var_name, ".RData"))


# Bind to one dataframe
gsm_data <- dplyr::bind_rows(gsm_calc)
head(gsm_data)
nrow(gsm_data)


## Prepare full dataframe
gsm_new <- gsm_data %>%
  separate(sample, c("A", "Cave1", "Cave2")) %>% 
  dplyr::select(., -"A")

head(gsm_new)


## Check length -- should match
length(unique(gsm_new$Cave1))
length(unique(gsm_new$Cave2))

length(cave_id1)
length(cave_id2)

## Append cave IDs
head(gsm_new)
gsm_new$caveID1 <- mapvalues(gsm_new$Cave1, from = unique(gsm_new$Cave1), to = cave_id1)
gsm_new$caveID2 <- mapvalues(gsm_new$Cave2, from = unique(gsm_new$Cave2), to = cave_id2)


## Check, select and rename variables 
head(gsm_new)
gsm_final <- gsm_new %>% 
  dplyr::select(caveID1, caveID2, sa, s10, ssk, sdr, srwi, sbi) %>%
  dplyr::rename(., sa_rough = sa, 
                s10z_rough = s10, 
                ssk_rough = ssk,
                sdr_rough = sdr,
                srwi_rough = srwi, 
                sbi_rough = sbi)


gsm_final
head(gsm_final)
nrow(gsm_final) ##  351

## Save data
var_name
write.csv(gsm_final, paste0("./5-LandGenAnalysis/results/surf_metrics/cha_SS/gsm_data_", var_name, "_", project_name, ".csv"), row.names = F)



## END OF SCRIPT ----


  

  
