###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## BIND SURFACE METRICS ----


## Clean Global Environment 
rm(list=ls())

## Load libraries
library(tidyverse)

## Set project name
project_name = "charinus_SS"


## 1. LOAD GENETIC DISTANCE ----
gen_dist <- read.csv("./5-LandGenAnalysis/results/gen_dist/cha/charinus_relatedness_27samples_OneSampleByCave.csv", sep = ",", h=T)
head(gen_dist)
nrow(gen_dist) ## 378
  
## 2. FILTER PAIRS OF CAVES ----
gen_dist <- gen_dist %>% 
            filter(ID1!= ID2) %>%
            filter(Cave1!= Cave2)

nrow(gen_dist) ## 351

## 4. SAVE FINAL DATAFRAME ----
## List all surface metrics
dir <- "./5-LandGenAnalysis/results/surf_metrics/cha_SS/"
list <- list.files(path=dir, pattern =".csv", full.names=TRUE)
list

## Load all surface metrics
sm_elev <- read.csv(list[1], sep = ",", h=T)
sm_fcover <- read.csv(list[2], sep = ",", h=T)
sm_rough <- read.csv(list[3], sep = ",", h=T)
head(sm_rough)
nrow(sm_rough) ## 351


## Merge dataframes by cave IDs
df1 <- merge(gen_dist, sm_elev, by.x = c("Cave1","Cave2"), by.y= c("caveID1","caveID2"))
df2 <- merge(df1, sm_fcover, by = c("Cave1","Cave2"), by.y= c("caveID1","caveID2"))
df_final <-merge(df2, sm_rough, by =  c("Cave1","Cave2"), by.y= c("caveID1","caveID2"))
head(df_final)
nrow(df_final) 

## Check length
length(unique(df_final$Cave1))
length(unique(df_final$Cave2))

setdiff(as.factor(df_final$Cave1), as.factor(df_final$Cave2)) ## "S11C_0159"
setdiff(as.factor(df_final$Cave2), as.factor(df_final$Cave1)) ## "S11B_0073"


## Save final dataframe
write.csv(df_final, 
          paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/final_SurfMetrics_", project_name, ".csv"), 
          row.names = F)


## END OF SCRIPT ----

