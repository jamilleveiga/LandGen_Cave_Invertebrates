###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

## EFFECTIVE POPULATION SIZE ----

#to install strataG:
#library(devtools)
#has_devel()
#devtools::install_github('ericarcher/strataG', build_vignettes = TRUE, force = T)

## REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

## LOAD THE PACKAGES:
library(r2vcftools)
library(adegenet)
library(vcfR)
library(dartR)
library(poppr)
library(strataG)

##6. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("source/functions_LanGen.R")
#source("source/r2vcftools_utils_2019.R")

## Set project name
project_name = "charinus"

## Load VCF file
snpsR = read.vcfR(paste0("vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"), verbose = T)
snpsR
# 82 samples; 1297 CHROMs; 2,291 variants; Object size: 17.5 Mb; 0 percent missing data

## Defining populations using the metafile:
snps =  vcfLink(paste0("vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"), overwriteID=T)
VCFsummary(snps) ## 82 individuals and 2291 SNPs.
population = as.factor(snps@meta$PopID_snmf)

##  Convert VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

## Add strata (pops) into Genind
snps_genind@pop = population

## Convert Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

## Convert Gtypes to GENEPOP to run in NeEstimator:
setwd("./3-GenDiversity/cha/NeEstimator/")
genepopWrite(snps_gtypes, paste0(project_name, "_genepop_82ind"))

setwd("./3-GenDiversity/cha/NeEstimator/Ne_82ind/")
genepopWrite(snps_gtypes, paste0(project_name, "_genepop_82ind"))

## END OF SCRIPT ----



## RUN NeEstimator

## 1. Download software from: http://www.molecularfisherieslaboratory.com.au/neestimator-software/
## 2. Activate NeEstimator
## 3. Run the software using "sudo java -jar ./NeEstimator2x1.jar"
## 4. Choose "Random mating" and Linkage Disequelibrium as method

