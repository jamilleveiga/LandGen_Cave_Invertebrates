###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

## POPULATION STRUCTURE ----

## REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


## INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline OR CREATE ONE FILE WITH THE FUNCTIONS. FUNCTIONS ARE IN THE NUMBER 9.


## GOALS FOR THIS STEP:
#A. ESTIMATE THE NUMBER OF GENETIC CLUSTERS IN THE DATASET. 
#B. CALCULATE GENETIC DIVERSITY IN POPULATIONS/CLUSTERS


## REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


## INSTALL AND LOAD THE PACKAGES ----

## LOAD SOURCE WITH FUNCTIONS TO BE USED ON THIS STEP.
source("source/r2vcftools_utils_2019.R")

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, SNPRelate, qvalue, poppr, mmod, dartR, vegan, ade4, vcfR, adegenet, seqinr, tess3r, ggplot2, tidyverse, reshape2) 

##6. CHOOSE A FOLDER FOR ANALYSIS. THE FOLDER "vcf" MUST BE THERE! 
#A. In RStudio go to SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. in RStudio tool bar or use the shortcut CTRL+SHIFT+H

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
dir_names = c("./2-Structure/1.snmf/cha/Results_Metafiles", 
              "./2-Structure/1.snmf/cha/Results_snmf_runs", 
              "./2-Structure/1.snmf/cha/Results_Diversity")

for (i in 1:length(dir_names)){
  if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
  else (print (paste0(dir_names[i], " has already been created. Be careful with overwritting")))}


## sNMF ANLYSIS ----

###1. CHOOSE A NAME FOR THE PROJECT. THE SAME ONE FROM FILTERING STEP:
#A. Project name:
project_name = "charinus"

###5. POPULATION ASSIGNMENT ANALYSIS USING sNMF WITHOUT FST OUTLIERS:
#Now we will carry out population assignment again, but using only neutral loci dataset.
#A. Load neutral .vcf file
snps_fil_ldF_neutral <- vcfLink(paste0("./vcf/cha/Neutral/rel_cutoff/", project_name, "_CutOffRelatedness_neutral_snps.vcf"), overwriteID = T) 
VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.
nrow(snps_fil_ldF_neutral@meta)

##B. Convert to geno object.You need to specify the output file. It will automatically subset the vcf file and assign it as a new object
snps_fil_ldF_neutral <- Geno(snps_fil_ldF_neutral, output.file = paste0("vcf/cha/Neutral/final/", project_name, "_filtered_ld_hw_neutral.geno"))
VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.

#B. Choose the 6 values of alpha. Alpha is the value of the regularization parameter (by default: 10). The results can depend on the value of this parameter, especially for small data sets. Less than 10,000 Snps we use values from 1000 to 10000. More than 10,000 SNP between 1 to 2,000.
# >= 10,000 SNPs
#alpha_values = c(1, 10, 100, 500, 1000, 2000)
# < 10,000 SNPs
alpha_values = c(1000, 2000, 4000, 6000, 8000, 10000)

#C. Create folders for alpha values and copy .geno object in each folder:
for (i in alpha_values){
  path = paste0("./2-Structure/1.snmf/cha/Results_snmf_runs/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("./vcf/cha/Neutral/final/", project_name, "_filtered_ld_hw_neutral.geno"), path )
}

#D. Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # numeber of replication by K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores

#E. Run SNMF (LEA) using different alpha.
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./2-Structure/1.snmf/cha/Results_snmf_runs/Alpha", i,"n/", project_name, "_filtered_ld_hw_neutral.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU, seed=123)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}



#F. To load the SNMF projects in a new R session (after quitting R), use:  project = load.snmfProject("Alpha1000//pilocarpus_filtered_ld_hw.snmfProject") ##This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./2-Structure/1.snmf/cha/Results_snmf_runs/Alpha", i,"n/", project_name, "_filtered_ld_hw_neutral.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop, "n"), pro_snmf)
}

#G. summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
summary(project5n)
summary(project6n)

#H. View Cross-Entropy plots
PlotK(project1n) ## 1
PlotK(project2n) ## 1
PlotK(project3n) ## 1
PlotK(project4n) ## 1
PlotK(project5n) ## 1
PlotK(project6n) ## 1

plot(project1n, lwd = 5, col = "red", pch=1, main = "alpha1000")  #1
plot(project2n, lwd = 5, col = "red", pch=1, main = "alpha2000")  #1
plot(project3n, lwd = 5, col = "red", pch=1, main = "alpha4000")  #1
plot(project4n, lwd = 5, col = "red", pch=1, main = "alpha6000")  #1
plot(project5n, lwd = 5, col = "red", pch=1, main = "alpha8000")  #1
plot(project6n, lwd = 5, col = "red", pch=1, main = "alpha10000") #1

#I. Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./2-Structure/1.snmf/cha/Results_snmf_runs/project",  i, "n.pdf"), onefile = F)
  path = paste0("./2-Structure/1.snmf/cha/Results_snmf_runs/Alpha", i,"n/", project_name, "_filtered_ld_hw_neutral.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
  }


#J. Select optimal K value
optimal_K = 1

#K. Select best run (lowest cross-entropy)
Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n, p4=project4n, p5=project5n, p6=project6n)
## [1] "Best run is: project = 1, run = 1" "Best run is: project = 2, run = 5"
## [3] "Best run is: project = 3, run = 9" "Best run is: project = 4, run = 1"
## [5] "Best run is: project = 5, run = 5" "Best run is: project = 6, run = 9"

project= project1n #edit this
run = 1 #edit this

#L. Barplot replace best run information
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
order = barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
write.table(order[[1]], "./2-Structure/1.snmf/cha/Results_Metafiles/snmf_bestK_pipegraph_order_neutral.txt")

pdf("./2-Structure/1.snmf/cha/Results_snmf_runs/bestk_BarChart_filtered_neutral.pdf", onefile = T) #replace best run information
my.colors <- rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors,
              xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)
dev.off()

## Add admixture coefficient and replace the population ID to vcf file
Qmat <- as.data.frame(Q(project, run=run, K=optimal_K))
head(Qmat)
nrow(Qmat)

columns = c() 

for (i in 1:ncol(Qmat)){
  columns[i] = paste0("Adx_Coeff_", i)
}

colnames(Qmat) = columns
head(Qmat)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:optimal_K){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat[i]
  }
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

popIds = apply(Qmat, 1, which.max)
snps_fil_ldF_neutral@meta$PopID_snmf <- popIds
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)


#O. Save new vcf file with sNMF results and pop ID
Save(snps_fil_ldF_neutral, paste0("./vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"))
write.csv(snps_fil_ldF_neutral@meta, paste0("./2-Structure/1.snmf/cha/Results_Metafiles/ancestry_coef_LEA_", project_name, ".csv"), quote = F)
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_snmf), file= paste0("./2-Structure/1.snmf/cha/Results_Metafiles/", project_name, "_neutral_popIDs_LEA.csv")) #save only pop ID from sNMF

VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.

## Reload VCF
snps_fil_ldF_neutral <- vcfLink(paste0("vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"), overwriteID = T) 
VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.
names(snps_fil_ldF_neutral@meta)

## Population size
length( snps_fil_ldF_neutral@meta$PopID_snmf[ snps_fil_ldF_neutral@meta$PopID_snmf==1]) ## 82

## Check cave IDs
caves1 <-  snps_fil_ldF_neutral@meta[snps_fil_ldF_neutral@meta$PopID_snmf==1,]
sort(unique(caves1$local_ID)) # 26 caves 

## END OF SCRIPT ----

