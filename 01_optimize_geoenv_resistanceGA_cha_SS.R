###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jamilleveiga/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## OPTIMIZE RESISTANCE SURFACES USING resistanceGA ----

## Clean Global Environment 
rm(list=ls())


## Install packages
#if(!("devtools" %in% list.files(.libPaths()))) {
#  install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
#} 

#devtools::install_github("wpeterman/ResistanceGA", build_vignettes = F)


## Load packages
library(ResistanceGA) # Carregar o pacote resitanceGA
library(GeNetIt)
library(reshape2)
library(tidyverse)

## Set project name
project_name = "charinus_SS"

## 1. LOAD DATA ----
gd <- read.csv(file = "./5-LandGenAnalysis/results/gen_dist/cha/charinus_GenDistMatrix_93samples.csv", header = T)
head(gd)
nrow(gd) ## 93
ncol(gd) ## 93
class(gd)

matrix.gen <- as.matrix(gd)
matrix.gen[1:7,1:7]

dim(matrix.gen)
nrow(matrix.gen)
ncol(matrix.gen)

## 2. LOAD COORDINATES ----
Sample.coord <- read.csv("./coords/cha/coords_charinus_93ind.csv", header=TRUE)
head(Sample.coord)
Sample.coord <- SpatialPoints(Sample.coord[,c(3,4)])
length(Sample.coord) ## 92

## 3. LOAD RASTER TO OPTIMIZE ----
landuse <- raster("./5-LandGenAnalysis/results/resistGA/cha_SS/asc/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
plot(landuse)
plot(Sample.coord, add=T)
res(landuse) # checar a resolucao do raster
extent(landuse) # checar a extensao do raster

## Check classes
unique(landuse[]) ## 70   8  12  13   6   2   3   7   5 100
sort(unique(landuse[]))

## CLASSES DO MOSAICO -- RASTER DE GEOAMBIENTES
## [1]   Buritizal
## [2]   Campo Brejoso ## OK
## [3]   Campo Graminoso ## OK
## [4]   Estrutura Relativa a Mineracao
## [5]   Lagoa ## OK
## [6]   Lajedo ## OK
## [7]   Mata Alta ## OK
## [8]   Mata Baixa ## OK
## [9]   Pastagem
## [10]  Samambaial
## [11]  Solo Exposto
## [12]  Vegetação Rupestre Aberta ## OK
## [13]  Vegetação Rupestre Arbustiva ## OK
## [70]  Floresta ## 70
## [100] Mineracao ## 100



## Check plot
Ricker.plot <- Plot.trans(PARM = c(5.46, 2456), ## Parametros para transformar a superficie continua - scale e maximum
                          Resistance = landuse, #
                          transformation= list(NA,"A"))



## RUN 1 (seed = 123) ----

## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res1/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 123,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run1.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run1.RData")

SS_RESULTS.gdist$CategoricalResults

## RUN 2 (seed = 456) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res2/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 456,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run2.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run2.RData")



## RUN 3 (seed = 789) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res3/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 789,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run3.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run3.RData")




## RUN 4 (seed = 321) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res4/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 321,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run4.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run4.RData")



## RUN 5 (seed = 654) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res5/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 654,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run5.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run5.RData")



## RUN 6 (seed = 987) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res6/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 987,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run6.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run6.RData")




## RUN 7 (seed = 213) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res7/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 213,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run7.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run7.RData")


## RUN 8 (seed = 546) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res8/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 546,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run8.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run8.RData")




## RUN 9 (seed = 879) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res9/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 879,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run9.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run9.RData")




## RUN 10 (seed = 101) ----
## 4. SET PATHS TO READ FILES AND SAVE RESULTS
read.dir   <- "./5-LandGenAnalysis/results/resistGA/cha_SS/asc/"
result.dir <- "./5-LandGenAnalysis/results/resistGA/cha_SS/res10/"

## 5. RESISTANCEGA SETTINGS

## Response variables
CS.response <- lower(matrix.gen)
class(CS.response)

## Genetic Algorithm inputs
GA.inputs <- GA.prep(ASCII.dir = read.dir,
                     seed = 101,
                     Results.dir = result.dir,
                     select.trans = list(NA,"A"),
                     parallel = 4)


## gdistance inputs
gdist.inputs <- gdist.prep(length(Sample.coord),
                            directions = 16,
                            samples = Sample.coord,
                            response = CS.response,
                            method = 'commuteDistance') ## Optimize using commute distance


## Results
SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                              GA.inputs = GA.inputs)


SS_RESULTS.gdist
save(SS_RESULTS.gdist , file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run10.RData")
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run10.RData")


## END OF OPTIMIZATION ----
            


## Load runs ----
load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run1.RData")
run1 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run2.RData")
run2 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run3.RData")
run3 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run4.RData")
run4 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run5.RData")
run5 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run6.RData")
run6 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run7.RData")
run7 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run8.RData")
run8 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run9.RData")
run9 = SS_RESULTS.gdist

load(file = "./5-LandGenAnalysis/results/resistGA/cha_SS/Results/results_cha_SS_run10.RData")
run10 = SS_RESULTS.gdist

## Compare IBD and null model with optimized surface
run1$AICc ## Distance 
run2$AICc ## Distance
run3$AICc ## Distance
run4$AICc ## Distance
run5$AICc ## Distance
run6$AICc ## Distance
run7$AICc ## optimized better than resistance
run8$AICc ## Distance
run9$AICc ## Distance
run10$AICc ## Distance

## Build dataframe w/ optimized values in each run

df = rbind(run1$CategoricalResults,
run2$CategoricalResults,
run3$CategoricalResults,
run4$CategoricalResults,
run5$CategoricalResults,
run6$CategoricalResults,
run7$CategoricalResults,
run8$CategoricalResults,
run9$CategoricalResults,
run10$CategoricalResults)


#df = run7$CategoricalResult

## Plots ----
names(df)

## Check classes
unique(landuse[]) ## 70   8  12  13   6   2   3   7   5 100
sort(unique(landuse[]))


## CLASSES DO MOSAICO -- RASTER DE GEOAMBIENTES
## [1]   Buritizal
## [2]   Campo Brejoso ## OK
## [3]   Campo Graminoso ## OK
## [4]   Estrutura Relativa a Mineracao
## [5]   Lagoa ## OK
## [6]   Lajedo ## OK
## [7]   Mata Alta ## OK
## [8]   Mata Baixa ## OK
## [9]   Pastagem
## [10]  Samambaial
## [11]  Solo Exposto
## [12]  Vegetação Rupestre Aberta ## OK
## [13]  Vegetação Rupestre Arbustiva ## OK
## [70]  Floresta ## 70
## [100] Mineracao ## 100


colnames(df) <- c("Surface", "obj.func_LL", "k", "AIC",
                  "AICc", "R2m", "R2c",  "LL",
                  "Campo Brejoso", 
                  "Campo Graminoso", 
                  "Lagoa", "Lajedo",
                  "Mata Alta", "Mata Baixa", 
                  "Vegetação Rupestre Aberta", "Vegetação Rupestre Arbustiva",
                  "Floresta", "Mineração")

## Save dotplot of first run

## Some data wrangling
library(tidyverse)

head(df)
df1 <- df %>% 
       summarise(camp_brej = mean(`Campo Brejoso`),
                 camp_gram = mean (`Campo Graminoso`),
                 lagoa = mean(Lagoa),
                 laj = mean(Lajedo),
                 mt_alt = mean (`Mata Alta`),
                 mt_bx = mean (`Mata Baixa`),
                 veg_abert = mean (`Vegetação Rupestre Aberta`),
                 veg_arb = mean (`Vegetação Rupestre Arbustiva`),
                 floresta = mean (Floresta),
                 minera = mean (Mineração))


head(df1)

df1$ID = rownames(df1)

mdata1 <- melt(df1, id=c("ID"))
head(mdata1)
levels(mdata1$variable)
sort(mdata1$value)


## Plot
optmized_values_mean <- 
  mdata1 %>%
  mutate(name = fct_reorder(variable, value)) %>%
  ggplot(aes(x=name, y=value)) +
  geom_point(color = "black", size = 5) +
  coord_flip() + ylim(1, 2000) +
  ylab("Value") + 
  xlab(NULL) + 
  theme(axis.text=element_text(size=12, color = "black", face = "bold"),
                          axis.title.y = element_text(size=15, color = "black", face = "bold"),
                          axis.title.x = element_text(size=15, color = "black", face = "bold"))


optmized_values_mean

## Save
pdf("./5-LandGenAnalysis/results/resistGA/cha_SS/optimized_values_cha_SS_mean.pdf")
optmized_values_mean
dev.off()



## MEAN VALUES OF OPTIMIZED RESISTANCES (10 RUNS) ----


### If boxplot corresponding to all ten runs is needed:

## Save boxplot correspondent to variation in 10 runs
head(df[, 9:18])
boxplot(df[, 9:18])

df1 = df[, 9:18]
df1$ID = rownames(df1)


mdata2 <- melt(df1, id=c("ID"))
head(mdata2)
levels(mdata2$variable)
max(mdata2$value) # 4809.004
min(mdata2$value) #1


## Boxplot 10 runs ----
optmized_values <- mdata2 %>%
  mutate(name = fct_reorder(variable, value)) %>%
  ggplot(aes(x=name, y=value)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", shape=20, size=6, color="blue", fill="blue")+
  coord_flip() + ylim(1, 5050) +
  ylab("Value") + 
  xlab(NULL) + 
  theme_minimal() + theme(axis.text=element_text(size=12, color = "black", face = "bold"),
                          axis.title.y = element_text(size=15, color = "black", face = "bold"),
                          axis.title.x = element_text(size=15, color = "black", face = "bold"))


optmized_values


## Save boxplot ----
pdf("./5-LandGenAnalysis/results/resistGA/cha_SS/optimized_values_cha_SS.pdf")
optmized_values
dev.off()


## Final raster with mean values ----

## Load rasters
r1 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res1/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r2 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res2/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r3 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res3/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r4 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res4/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r5 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res5/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r6 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res6/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r7 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res7/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r8 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res8/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r9 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res9/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")
r10 = raster("./5-LandGenAnalysis/results/resistGA/cha_SS/res10/Results/ASC_geoenvironment_UTM_250m_charinus_SS.asc")

## Stack rasters -----
## Build an optimized resistance surface with mean values correspondent to all runs
list_files <- c(r1,r2, r3, r4, r5, r6, r7, r8, r9, r10)
raster_stack <- stack(list_files)

## Calc mean values -----
raster_mean_values <- calc(raster_stack, fun = mean)
hist(raster_mean_values[])

## Check
unique(sort(raster_mean_values[]))
plot(raster_mean_values)

## Save raster ----
writeRaster(raster_mean_values, 
            filename = "./5-LandGenAnalysis/results/resistGA/cha_SS/mean_optimized_resist_values_geoenvrionment.tif", 
            format = "GTiff", overwrite=T)


## END OF SCRIPT ----