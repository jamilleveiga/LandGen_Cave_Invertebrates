###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

###############################################################################
############################### PRE-ANALYSIS ##################################
###############################################################################


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

## LOAD SOURCE WITH FUNCTIONS TO BE USED ON THIS STEP.
source("source/r2vcftools_utils_2019.R")

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, SNPRelate, qvalue, poppr, mmod, dartR, vegan, ade4, vcfR, adegenet, seqinr, tess3r, ggplot2, tidyverse, reshape2) 


##6. CHOOSE A FOLDER FOR ANALYSIS. THE FOLDER "vcf" MUST BE THERE! 
#A. In RStudio go to SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. in RStudio tool bar or use the shortcut CTRL+SHIFT+H

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
dir_names = c("./2-Structure/3.DAPC/cha/Results_Metafiles", 
              "./2-Structure/3.DAPC/cha/Results_DAPC", 
              "./2-Structure/3.DAPC/cha/Results_Diversity",
              "./2-Structure/3.DAPC/cha/Results_PCA")

for (i in 1:length(dir_names)){
  if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
  else (print (paste0(dir_names[i], " has already been created. Be careful with overwritting")))}


#######################################################################################
##################################### ANALYSES ########################################
#######################################################################################

##################################
############   DAPC  #############
##################################
#Here we will use other population clustering methods to assign individual
#to different population. This method is DAPC

#############################
#PCA is its ability to identify genetic structures in very large datasets within negligible computational time, and the absence of any assumption about the underlying population genetic model. PCA aims to summarize the overall variability among individuals, which includes both the divergence between groups (i.e., structured genetic variability), and the variation occurring within groups (‘random’ genetic variability).
#################################
#DA Discriminant Analysis focus on between-group variability, while neglecting within-group variation. this method also allows for a probabilistic assignment of individuals to each group, as in Bayesian clustering methods. the method requires the number of variables (alleles) to be less than the number of observations (individuals). This condition is generally not fulfilled in Single Nucleotide Polymorphism (SNP) or re-sequencing datasets. Second, it is hampered by correlations between variables, which necessarily occur in allele frequencies due to the constant-row sum constraint [i.e., compositional data]. Uncorrelated variables will be even more blatant in the presence of linkage disequilibrium
##################################
#DAPC relies on data transformation using PCA as a prior step to DA, which ensures that variables submitted to DA are perfectly uncorrelated, and that their number is less than that of analysed individuals (1:10 proportion). Without implying a necessary loss of genetic information, this transformation allows DA to be applied to any genetic data.
##################################
#K-means relies on the same model as DA to partition genetic variation into a between-group and a within-group component, and attempts to find groups that minimize the latter. we use Bayesian Information Criterion (BIC) to assess the best supported model, and therefore the number and nature of clusters.

###5. POPULATION ASSIGNMENT ANALYSIS USING DAPC:


###1. CHOOSE A NAME FOR THE PROJECT. THE SAME ONE FROM FILTERING STEP:
#A. Project name:
project_name = "charinus"

#A. Load neutral .vcf file without FST outliers
# LOAD VCF WITH ONLY NEUTRAL SNPS
snps_fil_ldF_neutral <- vcfLink(paste0("./vcf/cha/Neutral/rel_cutoff/", project_name, "_CutOffRelatedness_neutral_snps.vcf"), overwriteID = T) 
VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.

# LOAD SAME DATASET USING read.vcfR
vcf <- read.vcfR(paste0("./vcf/cha/Neutral/rel_cutoff/", project_name, "_CutOffRelatedness_neutral_snps.vcf"), verbose = FALSE)
vcf
#82 samples; 1297 CHROMs; 2,291 variants

#B. Convert "VCF" to "GENLIGHT"
my_genind_pop <- vcfR2genlight(vcf) # It's used for alleles numbers counted by individuals
my_genind_pop #total missing data 18.4 % %

#B. Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input

#C. Perform a PCA to choose the number of PC in the DAPC:
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE)

#D. % of PC variation
pc_pca = as.data.frame(pca_input$eig)
pc_pca[,2] = (pc_pca/sum(pc_pca))*100

#D1. Rule of 100% of variance:
index_100 = length(rownames(input@tab))-1
index_100 # number of PC to reach to 100% of explained variance

#D2. Rule of at least 95% of variance:
index_95 = length(which(cumsum(pc_pca[,2]) <= 95))
index_95

#D3. Rule of at least 70% of variance:
index_70 = length(which(cumsum(pc_pca[,2]) <= 70))
index_70 

#D4. Rule of at least 50% of variance:
index_70 = length(which(cumsum(pc_pca[,2]) <= 50))
index_70 

#D4. Rule of minimum variance per PCs:
variance_pc = 100/(nrow(input@tab)-1)
variance_pc #PCs that increase the explained variance bellow this threshould will be removed

#calculate number of PCs
index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
index_min

#E. Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust=40)
index_100 #For me, 100% of variation is more coherent.
index_95  
index_70  
index_min


## Index_100 ----
index_100
#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(123) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size ## 73  9

#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/Bestk_DAPC_index_100.pdf", onefile = T)
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Best K")
dev.off()

## Index_95 ----
index_95
#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(123) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size ## 73  9

#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/Bestk_DAPC_index_95.pdf", onefile = T)
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Best K")
dev.off()


## Index_70 ----
index_70
#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(123) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size ## 73  9

#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/Bestk_DAPC_index_70.pdf", onefile = T)
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Best K")
dev.off()


## Index_min ----
index_min

#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(123) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size ## 73  9

#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/Bestk_DAPC_index_min.pdf", onefile = T)
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Best K")
dev.off()

  
#G. ATTENTION, if your dataset is k = 1 go to PCA analysis (Action #4)! Graphs for DAPC need datasets with K >= 2

#H. Choose the best PCs number to recover correctly the clusters. The input should be a Genind object without missing data. Maximum number of PCs is number of individuals -1. Replication by deafaul is 30. Save automatically the graph
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/Best_PCs_Number_DAPC_index100.pdf", onefile = T)
number_PCs = xvalDapc(tab(input, NA.method = "mean"), grp$grp, scale = T, n.pca.max = (nrow(input@tab)-1))
dev.off()

#I. Verify the number of PCs and DA used and summary of DAPC
number_PCs$DAPC$n.pca
number_PCs$DAPC$n.da
summary(number_PCs$DAPC)

#J. Verify scatter plot and define the group colors and names:
# color and names
col_dapc= c('yellow', "darkslategrey")
legend_dapc = c("POP1", "POP2")

#plot graph
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = col_dapc, txt.leg = legend_dapc,
        clabel = FALSE, posi.leg = "topright", scree.pca = TRUE, pch=19:20,
        posi.pca = "topleft", posi.da = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

#K. Save scatterplots. Please edit the position of caption
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/scatter_DAPC.pdf", onefile = T)

# color and names
col_dapc= c('yellow', "darkslategrey")
legend_dapc = c("POP1", "POP2")

#plot graph
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = col_dapc, txt.leg = legend_dapc,
        clabel = FALSE, posi.leg = "topright", scree.pca = TRUE, pch=19:20,
        posi.pca = "topleft", posi.da = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

dev.off()

#L. Save a STRUCTURE-like graph 
pdf("./2-Structure/3.DAPC/cha/Results_DAPC/barchart_DAPC.pdf", onefile = T)
compoplot(number_PCs$DAPC, col=col_dapc, legend=TRUE, txt.leg=legend_dapc, cleg=.8)
dev.off()

#M. Add DAPC clusters ID in the VCF file
snps_fil_ldF_neutral@meta$PopID_DAPC = as.character(grp$grp)
head(snps_fil_ldF_neutral@meta)

#N. Add DAPC posterior probabilities in the VCF file:
Qmat_2 <- as.data.frame(number_PCs$DAPC$posterior)
head(Qmat_2)
columns = c() 
for (i in 1:ncol(Qmat_2)){
  columns[i] = paste0("DAPC_Posterior_", i)
}
colnames(Qmat_2) = columns
head(Qmat_2)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:length(grp$size)){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat_2[i]
}

head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)


# Check individuals and maximum Adx_Coeff and add to metafile:
popIds = apply(Qmat_2, 1, which.max)
snps_fil_ldF_neutral@meta$PopID_DAPC <- popIds
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

## Population size
length(snps_fil_ldF_neutral@meta$PopID_DAPC[snps_fil_ldF_neutral@meta$PopID_DAPC==1]) ## 72
length(snps_fil_ldF_neutral@meta$PopID_DAPC[snps_fil_ldF_neutral@meta$PopID_DAPC==2]) ## 20

## Check cave IDs
caves1 <- snps_fil_ldF_neutral@meta[snps_fil_ldF_neutral@meta$PopID_DAPC==1,]
length(sort(unique(caves1$local_ID))) # 23 caves

caves2 <- snps_fil_ldF_neutral@meta[snps_fil_ldF_neutral@meta$PopID_DAPC==2,]
length(sort(unique(caves2$local_ID))) ## 3 caves

## Veen Diagram
library(VennDiagram)
library(tidyverse)
library(seqinr)

A = unique(caves1$local_ID)
B = unique(caves2$local_ID)

## Creating a Venn Diagram ----
w1 <- venn.diagram(list(pop1=B, pop2=A),
                   lty = c("blank", "blank"),
                   fill = c("blue", "yellow"),
                   alpha = c(0.5, 0.5), cat.cex = 1.2, cex= 1.5,  cat.pos = 0, 
                   filename=NULL )

## The default plot
grid.newpage()
grid.draw(w1)

intersect(A, B) ## "S11D_0013"
caves2[caves2$local_ID=="S11D_0013",]
caves1[caves1$local_ID=="S11D_0013",]


#O. Save new vcf file and pop ID file
Save(snps_fil_ldF_neutral, paste0("vcf/cha/PopStructure/dapc/", project_name, "_filtered_ld_hw_neutral_POP_DAPC.vcf"))
write.csv(snps_fil_ldF_neutral@meta, paste0("./2-Structure/3.DAPC/cha/Results_Metafiles/ancestry_coef_DAPC_", project_name, ".csv"))
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_DAPC), file= paste0("./2-Structure/3.DAPC/cha/Results_Metafiles/", project_name , "_neutral_popIDs_DAPC.csv")) #edit column position by the number of coefficients.

## END OF SCRIPT ----

