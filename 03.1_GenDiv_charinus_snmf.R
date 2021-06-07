###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS TUTORIAL IS ORGANIZED IN DIFFERENT "Steps" AND INSIDE EACH STEP THERE ARE "Actions" DENOTED BY NUMBERS (#1) AND INSIDE EACH ACTION COULD EXIST "Operations" INDICATED BY LETTES (#A.)
#D. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


##2. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline


##3. GOALS FOR THIS STEP:
#A. CALCULATE GENETIC DIVERSITY INTRA AND INTER POPULATIONS/CLUSTERS
#A. CALCULATE GENETIC DISTANCE AMONG POPULATIONS/CLUSTERS AND INDIVIDUALS
#B. CALCULATE TAJIMA D FOR POPULATION EXPANSION


##4. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H


##5. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##6. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
#source("functions_LanGen.R")
source("source/r2vcftools_utils_2019.R")

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vegan, ecodist, vcfR, adegenet, poppr, mmod, reshape2, ggplot2, dartR)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
dir_names =  c("./3-GenDiversity/cha/Results_Diversity", 
               "./3-GenDiversity/cha/Results_Distance",
               "./3-GenDiversity/cha/Results_TajimaD")

for (i in 1:length(dir_names)){
  if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
  else (print (paste0(dir_names[i], " has already been created. Be careful with overwritting")))}


#------------------------------------------------------------------------------
#                        1. Loading Files
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "charinus"

#B. Choose the number of clusters according to step 2. In my case TESS3. 

###2.1. LOAD VCF FILES AND GEOGRAPHICAL INFORMATIONS: 
#A. Load neutral .vcf file with geographical information and genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"), overwriteID=T)

VCFsummary(snps_neutral) # 82 individuals and 2291 SNPs.
names(snps_neutral@meta) # check col names in metafile

#B. Number of cluster and method:
optimal_K = 1
method = "SNMF"

#B. Position of samples by populion by DAPC approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_snmf))){
  pop = which(snps_neutral@meta$PopID_snmf == i)
  assign(paste0("pop_SNMF_", i), pop)
}


#------------------------------------------------------------------------------
#                            2. Genetic Diversity
#------------------------------------------------------------------------------
###2.1. ESTIMATE GENETIC DIVERSITY:
#A. By species/taxon (all samples together)
Overall = GenDiv(snps_neutral)
write.csv(Overall, file=paste0("./3-GenDiversity/cha/Results_Diversity/Diversity_Overall_neutral_", project_name, "_", method, ".csv"))

#B. By individual
ind = GenDiv_IND(snps_neutral)
write.csv(ind, file=paste0("./3-GenDiversity/cha/Results_Diversity/Diversity_Individual_neutral_", project_name, "_", method, ".csv"))

#C. By cluster. ONLY IF YOUR K >= 2
#define a list of positions for each cluster
clusters = list(pop_SNMF_1)

#subset
for (i in 1:optimal_K){
  UNIND = snps_neutral@sample_id[clusters[[i]]]
  pop = Subset(snps_neutral, samples=UNIND)
  assign(paste0("pop_", i), pop)
}

#genetic diversity by cluster
pops = list(pop_1)

for (j in 1:optimal_K){
  cluster_div = GenDiv(pops[[j]])
  write.csv(cluster_div, file=paste0("./3-GenDiversity/cha/Results_Diversity/Diversity_Cluster_neutral_", project_name, "_POP", j, "_", method, ".csv"))
}


#------------------------------------------------------------------------------
#                              6. Tajima D
#------------------------------------------------------------------------------
###6.1. FILTER TO A SINGLE SNP PER CONTIG
#A. Load the VCF file:
snps_neutral = vcfLink(paste0("vcf/cha/PopStructure/snmf/", project_name, "_filtered_ld_hw_neutral_LEA.vcf"), overwriteID=T)
VCFsummary(snps_neutral) # 82 individuals and 2291 SNPs.

#B. A single SNP per contig. This thins SNPs to a given distance in bp from one another. Setting the distance higher than the length of the contig ensures that you'll have a single SNP per contig.
snps_thin = Filter(snps_neutral, filterOptions(thin=300)) 
VCFsummary(snps_thin) # 82 individuals and 1297 SNPs.

###6.2. SUBSET BY GENETIC CLUSTERS/POPULATIONS. SAME AS #3.C
#define a list of positions for each cluster
clusters = list(pop_SNMF_1)

#subset
for (i in 1:optimal_K){
  UNIND = snps_thin@sample_id[clusters[[i]]]
  pop = Subset(snps_thin, samples=UNIND)
  assign(paste0("popD_", i), pop)
}

# check
VCFsummary(popD_1) # 82 individuals and 1297 SNPs.
unique(popD_1@meta$local_ID)


###6.3. ESTIMATE TAJIMA'S D, BIAS-CORRECTED FOR MAF
set.seed(123)
tajd_p1 = TajimaD(popD_1, nboot=10000, maf=0.05, use_vcftools_D=FALSE)

###6.4. SAVE AND LOAD TAJIMA'S D RESULTS
#A. Save
save(tajd_p1, file = "./3-GenDiversity/cha/Results_TajimaD/tajd_p1_SNMF.Rdata")

#B. Load 
load("./3-GenDiversity/cha/Results_TajimaD/tajd_p1_SNMF.Rdata")

## Check results
tajd_p1$results # ## 'maf' may not agree with data: minimum count in data was 6 while minimum count allowed by maf filter was 5
as.data.frame(tajd_p1$results)
write.csv(as.data.frame(tajd_p1$results), file=paste0("./3-GenDiversity/cha/Results_TajimaD/TajimaD_Results_", project_name, "_", method, ".csv"))

###6.5. PLOT OBSERVED TAJIMA'D AGAINST THE NULL DISTRIBUTION
## The null distribution (histogram) is shown next to the observed Tajima's D value (red line)

#A. Create the theme for ggplot graphs
theme_pca = theme(legend.text = element_text(face = "italic",
                                             colour="black",
                                             family = "Helvetica",
                                             size = rel(1)), 
                  axis.title = element_text(colour="black",
                                            family = "Helvetica",
                                            size = rel(1.2)), 
                  axis.text = element_text(family = "Helvetica",
                                           colour = "black",
                                           size = rel(1)), 
                  axis.line = element_line(size = 1,colour = "black"), 
                  axis.ticks = element_line(colour="black",size = rel(1)),
                  
                  panel.grid.minor = element_blank(), 
                  panel.background = element_rect(fill = "whitesmoke"), 
                  panel.grid.major = element_line(colour="black",size = rel(0.2), linetype = "dotted"),
                  legend.key = element_blank(), 
                  legend.title = element_text(colour = "black",
                                              size = rel(1.5),
                                              family = "Helvetica"), 
                  plot.title = element_text(colour = "black",
                                            face = "bold",
                                            hjust = 0.5, #alingment
                                            size = rel(1.7),
                                            family = "Helvetica"))

#B.POP1
pdf(paste0("./3-GenDiversity/cha/Results_TajimaD/TajimaD_POP1_SNMF.pdf"), onefile = F)
ggplot(data.frame(x=tajd_p1$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + 
  geom_vline(xintercept=mean(tajd_p1$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + 
  geom_vline(xintercept=0) + theme_pca + ggtitle("Simulations for POP1 - SNMF") +
  labs(y= "Frequency", x = "Tajima's D") 

dev.off()


### END OF SCRIPT ----
