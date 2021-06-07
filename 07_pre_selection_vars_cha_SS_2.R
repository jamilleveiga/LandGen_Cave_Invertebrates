###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


## PRE-SELECTION OF PREDICTORS USING MODEL SELECTION ----

## Clean Global Environment 
rm(list=ls())

## Load packages
library(corMLPE)
library(nlme)
library(MuMIn)
library(usdm)
library(tidyverse)
library(snow)
library(parallel)


## Set project name
project_name = "charinus_SS"

## LOAD DATA ----
DT <- read.csv("./5-LandGenAnalysis/results/final_dataframe/cha_SS/dataframe2model_resDist+surfMetrics.csv", head=TRUE)
nrow(DT) ## 4059
ncol(DT) ## 46
names(DT)

## Check with histograms and plots
hist(DT$RELATEDNESS_AJK_Yang, breaks=100)
hist(DT$geoDist, breaks=100)
hist(DT$cd_geo_distance, breaks=100)

plot(cd_geo_distance ~ geoDist, data = DT)

## SCALE VARIABLES ----
## Scale independent variables before analysis

names(DT)
variables <- as.data.frame(scale(DT[, -c(1:6)]))
names(variables)

DT.final <- cbind(DT[ , c(1:6)], variables)
head(DT.final)
names(DT.final)

## Check histograms
hist(DT.final$RELATEDNESS_AJK_Yang, breaks=100)
hist(DT.final$geoDist, breaks=100)
hist(DT.final$cd_bio15, breaks=100)

## SET PREDICTORS BY GROUP ----
names(DT.final)

## Isolation by resistance
bioclim <- DT.final[ , c(7:9, 15:17)];  colnames(bioclim)
env_res <- DT.final[ , c(10,11,14)]; colnames(env_res)

## Islation by gradients
elev    <- DT.final[ , c(20:25)];       colnames(elev)
fcover  <- DT.final[ , c(26:31)];       colnames(fcover)
rough   <- DT.final[ , c(32:37)];       colnames(rough)

## Isolation by environment
env_diss <- DT.final[ ,c(39:46)]; colnames(env_diss)
geo_dist <- DT.final[ ,c(12,18,19)]; colnames(geo_dist)


## SET CLUSTER ----
## Run parallelized dredge
## by Ben Bolker in: https://stackoverflow.com/questions/30261473/mumin-pdredge-error-with-glmer

## Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

## clusterExport assigns the values on the master R process of the variables named in varlist to variables of the same names in the global environment (aka ‘workspace’) of each node. 
## The environment on the master from which variables are exported defaults to the global environment.
clusterExport(cluster,"DT.final")

## clusterEvalQ evaluates a literal expression on each cluster node. 
## It is a parallel version of evalq, and is a convenience function invoking clusterCall
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))



## 1. GROUP 1: BIOCLIMATICS ----

## BUILD FULL MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(bioclim), collapse = " + ")))

Fullmodel <- lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)


## RUN PDREDGE 
## Specify the number of predictor variables and including the max.r function
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                           m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_bioclim_cha_SS.RData")
  
### LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_bioclim_cha_SS.RData")

## BEST MODELS
nrow(Allmodels) ## 7
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 1.1. SELECTED VARS GROUP 1 ---- 
var1 = "cd_topoWet"


## 2. GROUP 2: ENVIRONMENTAL RESISTANCE ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(env_res), collapse = " + ")))

Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")

## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## ALL VARIABLE PLOTS
plot(DT.final$cd_elevation, RES) ; abline(0,0, col="red")
plot(DT.final$cd_forestcover, RES) ; abline(0,0, col="red")
plot(DT.final$cd_roughness, RES) ; abline(0,0, col="red") 

## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_ENVRES_cha_SS_2.RData")

### LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_ENVRES_cha_SS.RData")

## BEST MODELS
nrow(Allmodels) ## 4
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 2.1. SELECTED VARS GROUP 2 ----
var2 = "cd_elevation"

## end


## 3. GROUP 3: ELEVATION SURFACE METRICS ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(elev), collapse = " + ")))

Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")

## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)


## ALL VARIABLE PLOTS
plot(DT.final$sa_elev, RES) ; abline(0,0, col="red")
plot(DT.final$s10z_elev, RES) ; abline(0,0, col="red")
plot(DT.final$ssk_elev, RES) ; abline(0,0, col="red")
plot(DT.final$sdr_elev, RES) ; abline(0,0, col="red") 
plot(DT.final$srwi_elev, RES) ; abline(0,0, col="red")
plot(DT.final$sbi_elev, RES) ; abline(0,0, col="red")


## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_gsm_elev_cha_SS.RData")

### LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_gsm_elev_cha_SS.RData")

## BEST MODELS
nrow(Allmodels) ## 7
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 3.1. SELECTED VARS GROUP 3 ----
var3 = "s10z_elev"


## end


## 4. GROUP 4: VEGETATION COVER SURFACE METRICS ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(fcover), collapse = " + ")))
Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)


## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)


## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_gsm_fcov_cha_SS.RData")

## LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_gsm_fcov_cha_SS.RData")

## Final model selection table
nrow(Allmodels) ## 7
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 4.1. SELECTED VARS GROUP 4 ----
var4 = "ssk_fcov"

## end


## 5. GROUP 5: ROUGHNESS SURFACE METRICS ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(rough), collapse = " + ")))

Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)


## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_gsm_rough_cha_SS.RData")


## LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_gsm_rough_cha_SS.RData")

## Final model selection table
nrow(Allmodels) ## 7
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 5.1. SELECTED VARS GROUP 5 ----
var5 = "s10z_rough"

## end



## 6. GROUP 6: EPIGEAN ENVIRONMENT ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(env_diss), collapse = " + ")))

Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)


## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_envdiss_cha_SS.RData")

## LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_envdiss_cha_SS.RData")

## Final model selection table
nrow(Allmodels) ## 9
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 6.1. SELECTED VARS GROUP 6 ----
var6 = "d_B1000_Canga_pland_2019"

## end


## 7. GROUP 7: GEOGRAPHIC DISTANCE ----

## BUILD MODEL
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(geo_dist), collapse = " + ")))

Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## ALL VARIABLE PLOTS
plot(DT.final$geoDist, RES) ; abline(0,0, col="red")
plot(DT.final$cd_geo_distance, RES) ; abline(0,0, col="red")

## RUN PDREDGE
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_GEODIST_cha_SS.RData")

## LOAD ALL MODELS
load(file="./5-LandGenAnalysis/RData/Allmodels_GEODIST_cha_SS.RData")

## Final model selection table
nrow(Allmodels) ## 3
BM <- model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])

## 7.1. SELECTED VARS GROUP 7 ----
var7 = "topodist"

## end


## 8. SUBSET SELECTED VARIABLES ----
str(DT.final)
sel <- c(var1, var2, var3, var4, var5, var6, var7)
save(sel, file = "./5-LandGenAnalysis/RData/sel_cha_SS_2.RDATA")
load(file = "./5-LandGenAnalysis/RData/sel_cha_SS_2.RDATA")

sel

## Bind selected predictors with geographic and topographic distance
names(DT.final)
DT.save <- bind_cols(DT.final[  , c(1:6,38)], DT.final[  , sel])
head(DT.save)
ncol(DT.save) ## 14
names(DT.save)

## 9. SAVE FINAL DATAFRAME ----
project_name
write.csv(DT.save, "./5-LandGenAnalysis/results/final_dataframe/cha_SS/DT_final_MLPE_FULL_2.csv", row.names = F)


## END OF SCRIPT ----


