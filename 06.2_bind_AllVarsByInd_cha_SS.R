###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## BIND ALL DISTANCES - DATAFRAME ----

## Clean Global Environment 
rm(list=ls())

## Load libraries
library(tidyverse)

## Set project name
project_name = "charinus_SS"

## Load data
## 1. GENETIC DISTANCE DATA ----
gendist <- read.csv("./5-LandGenAnalysis/results/gen_dist/cha/charinus_relatedness_93samples.csv", sep = ",", h=T)
head(gendist)
nrow(gendist) ## 4371
length(unique(gendist$Cave1)) ## 27
length(unique(gendist$Cave2)) ## 27


gendist.ind <- gendist %>% 
    filter(ID1!= ID2) %>% 
    filter(Cave1!= Cave2)

nrow(gendist.ind)
head(gendist.ind)
length(unique(gendist.ind$Cave1)) ## 27
length(unique(gendist.ind$Cave2)) ## 27

## 2. RESISTANCE AND TOPODIST DATA ----

## Load resistance distance between individuals
res.ind <- read.csv(paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/final_ResDist_", project_name, ".csv"), sep= ",", h=T)
res.ind <- res.ind %>% filter(Cave1!=Cave2)
head(res.ind)
nrow(res.ind)
res <- res.ind[ ,5:15]

## Load topographic distance between individuals
topo.ind <- read.csv(paste0("./5-LandGenAnalysis/results/topo_dist/cha_SS/topdist_8_", project_name, ".csv"), sep= ",", h=T)
topo.ind <- topo.ind %>% filter(Cave1!=Cave2)
head(topo.ind)
nrow(topo.ind) 
topo <- topo.ind[,c(5,6)]


## Merge datasets
df1 <- gendist.ind %>% bind_cols(., res, topo) 

head(df1)
nrow(df1) 
names(df1)

## Nest resistance distances to bind with surface metrics
group1 <- df1 %>%
  filter(Cave1!= Cave2) %>%
  unite("Buffer_ID", Cave1:Cave2, sep = "_", remove = FALSE) %>%
  group_by(Buffer_ID) %>%
  nest(); nrow(group1) ## 543

head(group1)

group2 <- df1 %>%
  filter(Cave1!= Cave2) %>%
  unite("Buffer_ID", Cave2:Cave1, sep = "_", remove = FALSE) %>%
  group_by(Buffer_ID) %>%
  nest(); nrow(group2) ## 542

head(group2)

groups <- bind_rows(group1, group2) %>% 
                   unnest(cols = data) %>%
                   nest()
head(groups)
nrow(groups) ## 702

## 4. GRADIENT SURFACE METRICS DATA ----
gsm_final <- read.csv(paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/final_SurfMetrics_", project_name, ".csv"), sep=",", h=T)
head(gsm_final)
nrow(gsm_final) ## 351
ncol(gsm_final) ## 23

gsm <- as_tibble(gsm_final[ , -c(3:5)] %>% 
                   unite("Buffer_ID", Cave1:Cave2, sep = "_", remove = FALSE) %>%
                   dplyr::select(-c(Cave1, Cave2)))

head(gsm)
tail(gsm)
nrow(gsm) ## 351


## Filter Buffer_IDs from nested resistance distance data based on gsm
fil1 <- groups[groups$Buffer_ID %in% gsm$Buffer_ID, ]
nrow(fil1) ## 351
head(fil1)


## Join resistance distances and surface metrics
tbl1 <- inner_join(fil1, gsm,
                  by = "Buffer_ID")
str(tbl1)
head(tbl1)
nrow(tbl1) ## 351


## 5. ENVIRONMENTAL DISTANCE METRICS DATA ----
dist <- read.csv(paste0("./5-LandGenAnalysis/results/final_dataframe/cha_SS/Caves_Metrics_Dist_", project_name, ".csv"), h=T)
head(dist)
nrow(dist) ## 702

dist_metrics <- as_tibble(dist %>% 
                   unite("Buffer_ID", Cave1:Cave2, sep = "_", remove = FALSE) %>%
                   unite("Buffer_ID", Cave2:Cave1, sep = "_", remove = FALSE)) %>%
                   dplyr::select(-c(Cave1, Cave2))

head(dist_metrics)
nrow(dist_metrics) ## 702


## Filter Buffer_IDs from env distance metrics based on gradient surface metrics
nrow(tbl1) ## 351

fil2 <- dist_metrics[dist_metrics$Buffer_ID %in% tbl1$Buffer_ID, ]; nrow(fil2) ## 351
head(fil2)



## Join full dataset
final_tbl <- inner_join(tbl1, fil2,
                        by = "Buffer_ID")

### Unnest data
unnested <- final_tbl %>% 
  unnest(cols = data) 

unnested

nrow(unnested) ## 4059
head(unnested)
str(unnested)

## Check length caves
length(unique(unnested$Cave1)) ## 27
length(unique(unnested$Cave2)) ## 27

## Check length Individuals
length(unique(unnested$ID1)) ## 92
length(unique(unnested$ID2)) ## 92

## Check cols and rows
nrow(unnested) ## 4059
ncol(unnested) ## 46
names(unnested)

## Final to save
final2save <-  unnested

## Save full dataframe
project_name
write.csv(final2save, "./5-LandGenAnalysis/results/final_dataframe/cha_SS/dataframe2model_resDist+surfMetrics.csv", row.names =F)


## END OF SCRIPT ----

