###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## RUN MLPE ----    
     
## Clean Global Environment 
rm(list=ls())

## Load packages
library(corMLPE)
library(nlme)
library(MuMIn)
library(vegan)
library(usdm)
library(tidyverse)
library(stringr)
library(snow) 
library(parallel) 
library(effects)
library(RColorBrewer)
library(corrplot)


## Load source
source("./source/source_mlpe_models.R")


## Set project name
project_name = "charinus_SS"
  

## 1. LOAD DATA ----
DT.final <- read.csv("./5-LandGenAnalysis/results/final_dataframe/cha_SS/DT_final_MLPE_FULL_2.csv", h=T)
names(DT.final)
nrow(DT.final) ## 4059
ncol(DT.final) ## 14
names(DT.final)

hist(DT.final$RELATEDNESS_AJK_Yang)

## SET PRDICTORS ----
predictors <- DT.final[ , -c(1:6)]
names(predictors)


## SET CLUSTER ----
## Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

## Use clusterExport to send data to global environment (aka ‘workspace’) of each node. 
clusterExport(cluster,"DT.final")

## Load packages in the cluster
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))


## 2. RUN FULL MODEL ----

## Set formula
form <- as.formula(paste("RELATEDNESS_AJK_Yang ~ ", 
                         paste0(colnames(predictors), collapse = " + ")))

## Build full model
Fullmodel <- nlme::lme(form,
                       random = ~ 1|Buffer_ID,
                       correlation = corMLPE(form = ~ ID1 + ID2),
                       data = DT.final, method = "ML")


## Check model
summary(Fullmodel)
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## Correlations of estimates of model parameters ----
Fullmodel
corrFixed = summary(Fullmodel)$corFixed

## Rename vars
rownames(corrFixed) <- c("Intercept", 
                         "Sum of cave length",
                         "Topo Wetness resist",
                         "Elevation resist",
                         "Elevation s10z",
                         "Veget cover ssk", 
                         "Roughness s10z",
                         "Canga cover diss",
                         "Topographic distance")

colnames(corrFixed) <- c("Intercept", 
                         "Sum of cave length",
                         "Topo Wetness resist",
                         "Elevation resist",
                         "Elevation s10z",
                         "Veget cover ssk", 
                         "Roughness s10z",
                         "Canga cover diss",
                         "Topographic distance")
class(corrFixed)
dim(corrFixed)
corrFixed[2:9, 2:9]

## Corplot

## Save pdf 
pdf("./5-LandGenAnalysis/results/MLPE/cha_SS/Correlogram_cha_SS_2.pdf", onefile = T)

cols <- brewer.pal(10, "PuOr")

corrplot(corrFixed[2:9, 2:9], 
         #method="ellipse", 
         type = "upper", 
         outline = T, 
         diag = T, 
         tl.col = "black",
         col = cols)

dev.off()




## 3. RUN MODEL SELECTION WITH PARALLELIZED DREDGE ---- 
max.r(Fullmodel) ## 0.9545954

## Specify the number of predictor variables and including the max.r function
nrow(DT.final) ## 4059
options(na.action = na.fail)

## Run pdredge
start_time <- Sys.time()
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 8), extra = c(max.r), cluster)
end_time <- Sys.time()

## Compute run time
start_time1 = start_time
end_time1 = end_time
diff1 = end_time - start_time ## 1.859914 hours

## Save
save(Allmodels, file="./5-LandGenAnalysis/RData/Allmodels_lmeCaves_cha_SS_FULL_pdredge_2.RData")

## 3.1. LOAD ALLMODELS ----
load(file = "./5-LandGenAnalysis/RData/Allmodels_lmeCaves_cha_SS_FULL_pdredge_2.RData")


## 4. RETRIEVE MODELS WITH NOT COLLINEAR PREDICTORS ----
nrow(Allmodels) ## 256

## Retrieve Not Collinear Models, with max.r <=0.6

## Get Not Collinear Models using cluster too!
NCM <- get.models(Allmodels, subset = max.r <= 0.6, cluster)

  
## Save
save(NCM, file="./5-LandGenAnalysis/RData/NCM_lmeCaves_cha_SS_FULL_pdredge_2.RData")

## 4.1. LOAD NOT COLLINEAR MODELS ----
load(file="./5-LandGenAnalysis/RData/NCM_lmeCaves_cha_SS_FULL_pdredge_2.RData")
length(NCM) ##77


## 5. FINAL MODEL SELECTION TABLE ----
BM <- model.sel(NCM, rank=AIC)
nrow(BM) ## 59

## Check best models
BM[BM$delta <= 2, ]

## Save
save(BM, file="./5-LandGenAnalysis/RData/BM_lmeCaves_cha_SS_FULL_pdredge_2.RData")

## 5.1. LOAD RANKED BEST MODELS ----
load(file="./5-LandGenAnalysis/RData/BM_lmeCaves_cha_SS_FULL_pdredge_2.RData")

## Check best models
BM[BM$delta <= 2, ]

## Save FULL model selection table as dataframe ---- 
df_model_sel <- as.data.frame(BM)
head(df_model_sel)

## Save full tabl
write.csv(df_model_sel, "./5-LandGenAnalysis/files/originalwt_Best_cha_SS_full_table_2.csv", row.names = F)


## Original weight of best models
originalwt_Best <- df_model_sel[df_model_sel$delta <=2, ] ## #originalwt_Best <- df_model_sel[df_model_sel$delta !="NA",]
nrow(originalwt_Best) ## 2 models
write.csv(originalwt_Best, "./5-LandGenAnalysis/files/originalwt_Best_cha_SS_2.csv", row.names = F)


## Rescaled weight of best models
rescaledwt_Best <- BM[BM$delta <=2, ]
nrow(rescaledwt_Best) ## 2 models
write.csv(rescaledwt_Best, "./5-LandGenAnalysis/files/rescaledwt_Best_cha_SS_2.csv", row.names = F)


## MODEL AVERAGE ML ----
MAv <- model.avg(BM, subset = delta <=2)
summary(MAv)
confint(MAv) ## Confidence Intervals for model-averaged coefficients
as.data.frame(MAv$coefficients) ## Mean coefficients


## Save
save(MAv, file="./5-LandGenAnalysis/RData/ModelAverage_charinus_SS_2.RData")
write.csv(as.data.frame(confint(MAv)), row.names=TRUE, file= "./5-LandGenAnalysis/results/MLPE/cha_SS/BestModels_CI_coefficients_2.csv")
write.csv(as.data.frame(MAv$coefficients), row.names=TRUE, file="./5-LandGenAnalysis/results/MLPE/cha_SS/BestModels_mean_coefficients_2.csv")

## Save dataframe to plot estimates
df1 <- as.data.frame(MAv$coefficients[1, ])
rownames(df1) <- c("Intercept",
                   "Elevation resist", 
                   "Canga cover 1km buffer", 
                   "Elevation s10z", 
                   "Roughness s10z",
                   "Vegetation cover ssk", 
                   "Sum cave length")

df1$term <- rownames(df1)
colnames(df1) <- c("estimate", "term")
df2 <- bind_cols(df1$term, df1$estimate, as.data.frame(confint(MAv)[,1]), as.data.frame(confint(MAv)[,2]))
colnames(df2) <- c("term","estimate", "conf.low", "conf.high" )
rownames(df2) <- NULL
head(df2)

## Save
write.csv(df2, "./5-LandGenAnalysis/files/model_avg_df_cha_SS_2_2.csv", row.names = F)


## 6. MODEL VALIDATION ----
## AUTOCORRELATION OF RESIDUALS AND RESIDUALS X PREDICTORS

## Get top models separately and validate

## Top Model 1 ----
TopM1 <- get.models(BM, subset = 1)[[1]]
summary(TopM1)
intervals(TopM1)

## Check autocorrelation of residuals
RES <- residuals(TopM1, type="normalized")
FIT <- fitted(TopM1)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## Plot predictors x residuals
plot(DT.final$d_B1000_Canga_pland_2019, RES) ; abline(0,0, col="red")
plot(DT.final$s10z_elev, RES) ; abline(0,0, col="red")
plot(DT.final$s10z_rough, RES) ; abline(0,0, col="red")
plot(DT.final$ssk_fcov, RES) ; abline(0,0, col="red")


## Top Model 2 ----
TopM2 <- get.models(BM, subset = 2)[[1]]
summary(TopM2)
intervals(TopM2)

## Check autocorrelation of residuals
RES <- residuals(TopM2, type="normalized")
FIT <- fitted(TopM2)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## Plot predictors x residuals
plot(DT.final$d_B1000_Canga_pland_2019, RES) ; abline(0,0, col="red")
plot(DT.final$s10z_elev, RES) ; abline(0,0, col="red")
plot(DT.final$s10z_rough, RES) ; abline(0,0, col="red")
plot(DT.final$ssk_fcov, RES) ; abline(0,0, col="red")
plot(DT.final$sum_hp, RES) ; abline(0,0, col="red")

## 7. TEST SIGNIFICANCE OF PREDICTORS - LRT ----
## Run likelihood ratio tests

## TopM1
LRT <- drop1(TopM1, test="Chisq")
save(LRT, file = "./5-LandGenAnalysis//RData/LRT_TopM1_cha_SS_2.RData")
load(file = "./5-LandGenAnalysis//RData/LRT_TopM1_cha_SS_2.RData")

LRT

## Save likelihood ratio test
write.csv(LRT, file="./5-LandGenAnalysis/results/MLPE/cha_SS/LRT_Model1_FULL_pdredge_2.csv", row.names = TRUE)


## TopM2
LRT <- drop1(TopM2, test="Chisq")
save(LRT, file = "./5-LandGenAnalysis//RData/LRT_TopM2_cha_SS_2.RData")
load(file = "./5-LandGenAnalysis//RData/LRT_TopM2_cha_SS_2.RData")

LRT

## Save likelihood ratio test
write.csv(LRT, file="./5-LandGenAnalysis/results/MLPE/cha_SS/LRT_Model2_FULL_pdredge_2.csv", row.names = TRUE)



## 8. REFIT BEST MODELS USING REML ----
## Get accurate estimates

## TopM1
TopM1

Refit_TopM1 <- nlme::lme(RELATEDNESS_AJK_Yang ~ 
                           cd_elevation + d_B1000_Canga_pland_2019 +
                           s10z_elev + s10z_rough + ssk_fcov,
               random = ~ 1|Buffer_ID,
               correlation = corMLPE(form = ~ ID1 + ID2), 
               data = DT.final, method = "REML")

summary(Refit_TopM1)

## Save
save(Refit_TopM1, file="./5-LandGenAnalysis/results/MLPE/cha_SS/refit_REML_TopM1_2.RData")
load(file="./5-LandGenAnalysis/results/MLPE/cha_SS/refit_REML_TopM1_2.RData")

## Plot All effects (predictors) TopM1 ----
ae_TopM1 <- allEffects(Refit_TopM1)
plot(ae_TopM1, rescale.axis=F)
plot(ae_TopM1,type="response")

## Save
pdf("./5-LandGenAnalysis/results/MLPE/cha_SS/TopM1_AllEffectsPlot_cha_SS_2.pdf",
    width = 10,
    height = 7.5)
plot(ae_TopM1,type="response")
dev.off()


## TopM2
TopM2

Refit_TopM2 <- nlme::lme(RELATEDNESS_AJK_Yang ~ 
                         cd_elevation + d_B1000_Canga_pland_2019 +      
                         s10z_elev + s10z_rough + ssk_fcov + sum_hp,
                         random = ~ 1|Buffer_ID,
                         correlation = corMLPE(form = ~ ID1 + ID2), 
                         data = DT.final, method = "REML")

summary(Refit_TopM2)

## Save
save(Refit_TopM2, file="./5-LandGenAnalysis/results/MLPE/cha_SS/refit_REML_TopM2_2.RData")
load(file="./5-LandGenAnalysis/results/MLPE/cha_SS/refit_REML_TopM2_2.RData")

## Plot All effects (predictors) TopM2 ----
ae_TopM2 <- allEffects(Refit_TopM2)
plot(ae_TopM2, rescale.axis=F)
plot(ae_TopM2,type="response")

## Save
pdf("./5-LandGenAnalysis/results/MLPE/cha_SS/TopM2_AllEffectsPlot_cha_SS_2.pdf",
    width = 10,
    height = 7.5)
plot(ae_TopM2,type="response")
dev.off()


## MODEL AVERAGE REML ----
MAv_reml <- model.avg(Refit_TopM1, Refit_TopM2)
summary(MAv_reml)
confint(MAv_reml) ## Confidence Intervals for model-averaged coefficients
as.data.frame(MAv_reml$coefficients) ## Mean coefficients


## Save
save(MAv_reml, file="./5-LandGenAnalysis/RData/ModelAverage_REML_charinus_SS_2.RData")
write.csv(as.data.frame(confint(MAv_reml)), row.names=TRUE, file= "./5-LandGenAnalysis/results/MLPE/cha_SS/BestModels_CI_coefficients_REML_2.csv")
write.csv(as.data.frame(MAv_reml$coefficients), row.names=TRUE, file="./5-LandGenAnalysis/results/MLPE/cha_SS/BestModels_mean_coefficients_REML_2.csv")


## Save dataframe to plot estimates
df1 <- as.data.frame(MAv_reml$coefficients[1, ])
rownames(df1) <- c("Intercept",
                   "Elevation resist",
                   "Canga cover 1km buffer",
                   "Elevation s10z", 
                   "Roughness s10z","Vegetation cover ssk", 
                   "Sum cave length")

df1$term <- rownames(df1)
colnames(df1) <- c("estimate", "term")
df2 <- bind_cols(df1$term, df1$estimate, as.data.frame(confint(MAv_reml)[,1]), as.data.frame(confint(MAv_reml)[,2]))
colnames(df2) <- c("term","estimate", "conf.low", "conf.high" )
rownames(df2) <- NULL
head(df2)


## Save
write.csv(df2, "./5-LandGenAnalysis/files/model_avg_df_REML_cha_SS_2.csv", row.names = F)

## Compare model average ML x REML ----
summary(MAv)
summary(MAv_reml)

## END OF SCRIPT ----

