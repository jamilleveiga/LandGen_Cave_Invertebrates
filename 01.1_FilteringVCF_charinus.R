###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


### GOAL: VCF FILTERING - PRE-ANALYSIS ----

##2. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER POPULATIONS STEP IN STACKS ORGANIZED BY POPULATIONS/GROUPS OR A REGULAR ".VCF" OR ".GVCF" FROM ANOTHER SOFTWARE.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline OR CREATE ONE FILE WITH THE FUNCTIONS. FUNCTIONS ARE IN THE NUMBER 9.
#C. THE FILE ".CSV" WITH THE GEOGRAPHICAL COORDENATES OF GENETIC SAMPLES IN DECIMALS DEGREES PREFERENCIALLY AND ONE COLUMN WITH LOCALITIES ID. 


##3. GOALS FOR THIS STEP:
#A. FILTER RAW SNPS DATASET INTO NEUTRAL SNPS DATASET. 

##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##  5. INSTALL AND LOAD THE PACKAGES
#A. install the packages automatically
if("remotes" %in% rownames(installed.packages()) == FALSE){install.packages("remotes")
} else {print (paste0("'remotes' has already been installed in library"))}
if("BiocManager" %in% rownames(installed.packages()) == FALSE){install.packages("BiocManager")
} else {print (paste0("'BiocManager' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("devtools" %in% rownames(installed.packages()) == FALSE){install.packages("devtools")
} else {print (paste0("'devtools' has already been installed in library"))}

if("r2vcftools" %in% rownames(installed.packages()) == FALSE){remotes::install_github("nspope/r2vcftools")
} else {print (paste0("'r2vcftools' has already been installed in library"))}
if("LEA" %in% rownames(installed.packages()) == FALSE){BiocManager::install("LEA")
} else {print (paste0("'LEA' has already been installed in library"))}
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")
} else {print (paste0("'tidyverse' has already been installed in library"))}
  
#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, tidyverse, vcfR, dartR, poppr, adegenet, 
               LEA, pcadapt, qvalue, psych, tess3r, seqinr)

##6. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H
#B. CREATE A FOLDER NAMED "vcf" INSIDE YOUR WORKING DIRECTORY AND COPY THE .vcf OR .gvcf THERE. "vcf/pilocarpus.vcf"
#C. CREATE A FOLDER NAMED "coords" INSIDE YOUR WORKING DIRECTORY AND COPY THERE, THE FILE .csv WITH GEOGRAPHICAL INFORMATION OF THE SAMPLES

##7. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("source/r2vcftools_utils_2019.R")

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
dir_names = c("./1-Filtering/cha/",
              "./1-Filtering/cha/ld", 
              "./1-Filtering/cha/Results_Metafiles", 
              "./1-Filtering/cha/Results_snmf",
              "./1-Filtering/cha/Results_tess3",
              "./1-Filtering/cha/Results_pcadapt",
              "./1-Filtering/cha/FST-Outliers/",
              "./1-Filtering/cha/Results_PCA/",
              "./vcf/cha")

for (i in 1:length(dir_names)){
  if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
  else (print (paste0(dir_names[i], " has already been created. Be careful with overwritting")))}

#######################################################################################
##################################### ANALYSES ########################################
#######################################################################################


###1. CHOOSE A NAME FOR THE PROJECT AND THE PATH FOR VCF AND COORDS FILES:

#A. Project name:
project_name = "charinus"

#B. vcf file path:
vcf_file = "raw_vcf/cha/cavespider95.vcf"

#C. Coordenates file path and load it:
coord_file = "coords/cha/coords_gen_utm_charinus95.csv"
coords = read.csv(coord_file)
head(coords)
tail(coords)

#D. Define the columns positions of localities ID, sample ID, longitude, and latitude. Add more columns if you need
sample_ID = 1
local = 2
longitude = 3
latitude = 4

###2. LOAD THE VCF FILE
#A. VCF
snps_raw <-  vcfLink(vcf_file, overwriteID=T)
snps_raw

#B. Remove loci that are not SNPs
snps <- Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2), indels="remove")
VCFsummary(snps) # 96 individuals and 18055 SNPs.

#C. Basic stats sample size (SZ)
raw_SZ = capture.output(VCFsummary(snps_raw)) 
filtered_SZ = capture.output(VCFsummary(snps)) 
raw_SZ # 96 individuals and 21062 SNPs.
filtered_SZ # 96 individuals and 18055 SNPs

## Basic statistics
snps@sample_id
snps@site_id
snps@meta
Chrom(snps) ##Chromosome, possitions, and IDs

###3. REMOVE INDIVIDUALS IF IT IS NECESSARY
#A. Remove one individual
#snps@meta$sample_name
#UNIND <- snps@sample_id[snps@sample_id != "9202_rep_sorted"]
#snps_unind <- Subset(snps, samples=UNIND)
#snps_unind@meta
#VCFsummary(snps_unind) ## 290 individuals and 35214 SNPs.

#B. Remove more than one individual. In the case Granito locality (expect for 9521) and individuals (9592, 9541, 9561) and with coefficients <50% (9586 9599, 9507, 9519, 9511, 9508):
snps@meta$sample_name
remove_ind = "CH_245_rep_sorted"
UNIND <- snps@sample_id[!snps@sample_id %in% remove_ind]
snps_unind <- Subset(snps, samples=UNIND)

#verify sample size
length(snps@meta$sample_name)-length(remove_ind)
length(snps_unind@meta$sample_name)
VCFsummary(snps_unind) ## 95 individuals and 18055 SNPs.

#C. If you do not need to remove individuals:
#snps_unind = snps
#VCFsummary(snps_unind)

###4. CHECK THE QUALITY OF DATA
#A. See GenotypeMatrix
genotypes <- GenotypeMatrix(snps_unind) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing; otherwise, gives the number of derived alleles in the genotype -- ie. a 0 and 2 are both homozygotes

#B. Look at depth, quality, HWE, HE, allele frequencies, and Pi
site.depth <- Query(snps_unind, type="site-mean-depth") #avaliar coertura
summary(site.depth$MEAN_DEPTH) #Mean = 59.43/ Median = 45.45
hist(site.depth$MEAN_DEPTH, breaks=30) 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <100]) ### >20 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <1000]) ### <800

quality <- Query(snps_unind, type="site-quality")
summary(quality$QUAL) # Mean = 45.00/ Median = 46.66
hist(quality$QUAL) ## >30

PI <- Query(snps_unind, type="site-pi")
mean(PI$PI) ## Mean nucleotide divergency per-site ### 0.237
hist(PI$PI)

HWE <- Query(snps_unind, type="hardy")
summary(HWE$P_HWE) #Mean = 0.168 / Median = 0.0003339
hist(HWE$P_HWE)
hist(HWE$P_HWE[HWE$P_HWE<0.0001])

HE <- Query(snps_unind, type="het")
hist(HE$O.HOM) ## O.HOM, E.HOM, N_SITES, F
hist(HE$E.HOM) 
hist(HE$N_SITES) 
hist(HE$F) ## Inbreeding coefficient

freq <- Query(snps_unind, type="freq2")
hist(freq$X.FREQ.)
hist(freq$X.FREQ.1)
head(freq)

#C. Missing per locus
Missing <- apply(GenotypeMatrix(snps_unind), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 66.32%
hist(Missing)

#D. Missing per individual
Missing_ind <- apply(GenotypeMatrix(snps_unind),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 66.55%
hist(Missing_ind)

###5. ADD COORD INFORMATION AND MISSING DATA AMOUNT TO A VCF FILE:
#A. Save as data frame
Missingind.df = as.data.frame(Missing_ind)
#row.names(Missingind.df) = str_remove(row.names(Missingind.df), pattern="_sorted") #remove string "_sorted" from samples id if you need it
Missingind.df$ID = row.names(Missingind.df)
Missingind.df

#B. Remove the individuals in the coord files that you excluded in step #3 
missing_coord = coords[coords$vcf_ID %in% Missingind.df$ID, ]
head(missing_coord)
tail(missing_coord)
length(missing_coord[,1]) ## 95

#C. When coords are in UTM you need to transform to lon/lat in decimals. Use this code:
#C1. remove NAs if they are present.
#missing_coord = missing_coord[!is.na(missing_coord[,longitude]),]

#C2. Add a datum
#My UTM data are in SAD69 UTM zone 22S. Search in internet the EPSG code for that area. For example https://epsg.io/ 
#In the case is EPSG:32722.
#coord_gen_pts = missing_coord[,c(longitude, latitude)] #select lat and long
#coordinates(coord_gen_pts) <- coord_gen_pts
#projection(coord_gen_pts) = crs("+init=epsg:32722") #define utm datum only for coord cols

#C3. Convert:
#coords_gen_dec = spTransform(coord_gen_pts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"))
#create a data.frame with converted coords:
#x = missing_coord[, c(local, sample_ID)]
#y = as.data.frame(coords_gen_dec)
#y = y[,c(longitude, latitude)]
#coords_gen = cbind(x, y)
#colnames(coords_gen) = c("local_ID", "samples", "long", "lat")
#coords_dec = coords_gen
#head(coords_dec)
#tail(coords_dec)
#missing_coord = coords_dec

#D. Sort the coord file according to the missing file:
missing_coord_sorted = missing_coord[order(match( missing_coord$vcf_ID, Missingind.df$ID)),]
head(missing_coord_sorted)

#E. Check the differences. Should be "character (0)"
setdiff(Missingind.df$ID, missing_coord_sorted[,sample_ID])

#F. Check if identical samples were selected. Should be "TRUE"
identical(as.character(Missingind.df$ID), as.character(missing_coord_sorted[,sample_ID]))

#G. Create a  data frame merging missing and locality information
missing_local = cbind(missing_coord_sorted[,c(local, longitude, latitude)], Missingind.df)
head(missing_local)
tail(missing_local)

#H. Add coords and missing information to vcf file. Be careful to the order
#verify order
identical(as.character(missing_local$ID), as.character(snps_unind@meta$sample_name))

#add the information
snps_unind@meta$local_ID = missing_local[,1] #local ID
snps_unind@meta$longitude = missing_local[,2] #longitude
snps_unind@meta$latitude = missing_local[,3] #latitude
snps_unind@meta$missing_PC = missing_local[,4] #% of missing data
snps_unind@meta$ind_ID = missing_local[,5] #sample ID

#verify
head(snps_unind@meta) #verify
tail(snps_unind@meta) #verify

#I. Save missing results and coord information:
#write.csv(snps_unind@meta, paste0("./1-Filtering/cha/Results_Metafiles/missing_data_coords", project_name, ".csv"))

###6. REMOVE INDIVIDUALS WITH LARGE AMOUNT OF MISSING GENOTYPES

# define threshold for individual missing data tolerance - shows individuals that should be removed from ID data
#Missing per locus
Missing <- apply(GenotypeMatrix(snps_unind), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 66.32
Missing.df <- as.data.frame(Missing)
Missing.df$ID <- row.names(Missing.df)

# define threshold for missing data tolerance per locus
## Remove loci with large amounts of missing (set to max 30% missing)
locustokeep <- Missing.df[Missing.df$Missing <= 30,] ##
snps_lowlocusmiss <- Subset(snps_unind, sites=locustokeep$ID)
snps_lowlocusmiss@meta
VCFsummary(snps_unind) # 95 individuals and 18055 SNPs.
VCFsummary(snps_lowlocusmiss) # 95 individuals and 17310 SNPs.

#Missing per individual
Missing_ind <- apply(GenotypeMatrix(snps_lowlocusmiss), 1, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing_ind) ##Max missing = 66.724
hist(Missing_ind)
Missingind.df <- as.data.frame(Missing_ind)
Missingind.df$ID <- row.names(Missingind.df)

#define threshold for individual missing data tolerance 
Missingind.df[Missingind.df$Missing_ind>50, ] 

# These are the individuals that need to be excluded 
#Missing_ind    ID
#CH_255_sorted    66.72444 CH_255_sorted
#CH_256_sorted    62.29347 CH_256_sorted

##Remove individuals with large amounts of missing (set to max 50% missing)
## Set value
value = 50

indtokeep <- Missingind.df[Missingind.df$Missing_ind<=value, ]
length(indtokeep$ID) # 93
snps_ind_pc  <- Subset(snps_lowlocusmiss, samples=indtokeep$ID) 
snps_ind_pc@meta
VCFsummary(snps_unind) ## 95 individuals and 18055 SNPs.
VCFsummary(snps_lowlocusmiss) # 95 individuals and 17310 SNPs.
VCFsummary(snps_ind_pc) # 93 individuals and 17310 SNPs.

removed_missing_SZ = capture.output(VCFsummary(snps_ind_pc)) 
removed_missing_SZ ## 93 individuals and 17310 SNPs.

###################################################################################
################# FILTER NEUTRAL LOCI by depth, quality, maf ######################
############ and removing high LD loci within same and between contigs ############
###################################################################################

###11. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE:
#A. from dataset with low missing by individuals

snps_fil_hwe_low <- Filter(snps_ind_pc, filterOptions(minQ=30, maf=0.03, 
                                          min.meanDP=20, max.meanDP=800)) 

VCFsummary(snps_fil_hwe_low) # 93 individuals and 15088 SNPs.

#B. from dataset with all individuals
#snps_fil_hwe_high <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.8, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
#VCFsummary(snps_fil_hwe_high)

#C. choose one dataset:
snps_fil_hwe = snps_fil_hwe_low
neutralFilter_SZ = capture.output(VCFsummary(snps_fil_hwe))

#D. Define R² value:
r2 = 0.8

###12. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:
## If you have more than one population, use the code below. This code identifies SNPs deviating from HW equilibrium within each population, and then removes those SNPs that are in desequilibrium in all populations.

#A. remove snps with R² value
ld_within2 <- Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=r2))
head(ld_within2)
hist(ld_within2$R.2)
#save(ld_within2, file = "1-Filtering/cha/ld/ld_within2.RData")
load("1-Filtering/cha/ld/ld_within2.RData")

#B. Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within2$ID1
nold_snps <- snps_fil_hwe@site_id[!(snps_fil_hwe@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil_hwe, sites=nold_snps) # Keep snps that are not in LD.
neutralLDWithin_SZ = capture.output(VCFsummary(snps_fil_ld))
neutralLDWithin_SZ ## 93 individuals and 10468 SNPs.

###13. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) BETWEEN CONTIGS:
#A. remove snps with R²
ld_between <- Linkage(snps_fil_ld, type="interchrom-geno-r2", linkageOptions(min.r2=r2)) 
head(ld_between)
hist(ld_between$R.2)
#save(ld_between, file = "1-Filtering/cha/ld/ld_between.RData")
load("1-Filtering/cha/ld/ld_between.RData")

#B. Select one set of the correlated snps (ID1 or ID2)
ld2_snps <- ld_between$ID1
nold2_snps <- snps_fil_ld@site_id[!(snps_fil_ld@site_id %in% ld2_snps)]
snps_fil_ld2 <- Subset(snps_fil_ld, sites=nold2_snps) # Keep snps that are not in LD.
neutralLDBetween_SZ = capture.output(VCFsummary(snps_fil_ld2)) 
neutralLDBetween_SZ ## 93 individuals and 4835 SNPs.

## Filter by HWE
snps_fil_ldF <- Filter(snps_fil_ld2, filterOptions(hwe=0.0001)) 
VCFsummary(snps_fil_ldF) # 93 individuals and 3093 SNPs.

neutralFilter_HWE_SZ2 = capture.output(VCFsummary(snps_fil_ldF))
neutralFilter_HWE_SZ2 # 93 individuals and 3093 SNPs.

###14. VERIFY THE QUALITY OF THE FILTERED DATASET:
#A. Quality indexes:
site.depth2 <- Query(snps_fil_ldF, "site-mean-depth")
quality2 <- Query(snps_fil_ldF, "site-quality")
HWE2 <- Query(snps_fil_ldF, type="hardy")

hist(site.depth2$MEAN_DEPTH, breaks=30)
hist(site.depth2$MEAN_DEPTH[site.depth2$MEAN_DEPTH <200])
hist(quality2$QUAL)
hist(HWE2$P_HWE[HWE2$P_HWE<200])

#B. Real value of missing data by individuals:
Missing_ind <- apply(GenotypeMatrix(snps_fil_ldF),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 43.55%
hist(Missing_ind)

#C. Missing per locus:
Missing <- apply(GenotypeMatrix(snps_fil_ldF), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 30.11%
hist(Missing)

VCFsummary(snps_unind) # 95 individuals and 18055 SNPs.
VCFsummary(snps_fil_ldF) ## 93 individuals and 3093 SNPs.

###15. SAVE THE .VCF FILE:
Save(snps_fil_ldF, paste0("vcf/cha/Neutral/parcial/", project_name,"_filtered_ld_hwe.vcf"))


#######################################################################################
################### FILTER OUTLIER SNPS BY FST OUTLIER TESTS  #########################
#######################################################################################

#------------------------------------------------------------------------------
#                   Filtering Neutral Loci - TESS3 
#------------------------------------------------------------------------------
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html

###7.1. INPUT FILES FOR TESS:
#A. Load the .VCF file with only neutral SNPs:
project_name = "charinus"
snps <-  vcfLink(paste0("vcf/cha/Neutral/parcial/", project_name,"_filtered_ld_hwe.vcf"), overwriteID=T)
VCFsummary(snps) ## 93 individuals and 3093 SNPs.

#B. Create a Genotype matrix
genotypes = GenotypeMatrix(snps) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing;
class(genotypes)
genotypes = replace(genotypes, genotypes == -1, NA)

#C. Create a Matrix with long and lat 
coordinates = snps@meta[,4:5]
class(coordinates)
coordinates = data.matrix(coordinates, rownames.force = NA)
class(coordinates)

#verify the coords
plot(coordinates, pch = 19, cex = .5, xlab = "Longitude", ylab = "Latitude")


###7.2. RUNNING THE TESS3R FUNCTION:
#A. Costumize values of lambda
lambda_values = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)

#B.Run tess for all lambda values

### LAMBDA = 0
set.seed(123)
TESS3_Lambda_0 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 0,
                   method = "projected.ls", max.iteration = 5000, rep = 10,
                   ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_0, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_0.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_0, pch = 19, col = "blue", type="b",lwd=2,
       crossvalid = T, crossentropy = T,
       xlab = "Number of ancestral populations", cex.lab=1.5,
       ylab = "Cross-validation score", main="lambda = 0")
dev.off()


### LAMBDA = 0.25
set.seed(123)
TESS3_Lambda_0.25 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 0.25,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_0.25, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.25.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_0.25.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_0.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 0.25")
dev.off()


### LAMBDA = 0.5
set.seed(123)
TESS3_Lambda_0.5 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 0.5,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_0.5, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.5.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_0.5.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_0.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 0.5")
dev.off()


### LAMBDA = 0.75
set.seed(123)
TESS3_Lambda_0.75 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 0.75,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_0.75, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.75.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_0.75.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_0.75, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 0.75")
dev.off()

### LAMBDA = 1
set.seed(123)
TESS3_Lambda_1 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 1,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_1, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_1.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_1, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 1")
dev.off()


### LAMBDA = 1.25
set.seed(123)
TESS3_Lambda_1.25 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 1.25,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_1.25, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.25.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_1.25.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_1.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 1.25")
dev.off()


### LAMBDA = 1.5
set.seed(123)
TESS3_Lambda_1.5 = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = 1.5,
                       method = "projected.ls", max.iteration = 5000, rep = 10,
                       ploidy = 2, openMP.core.num = 4)

save(TESS3_Lambda_1.5, file = paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.5.RData"))

pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_Lambda_1.5.pdf"), onefile =F)
plot.new()
plot(TESS3_Lambda_1.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main="lambda = 1.5")
dev.off()



### Load lambda projects
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.25.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.5.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_0.75.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.25.RData"))
load(paste0("./1-Filtering/cha/Results_tess3/TESS3_Lambda_1.5.RData"))


#C. Plot results for all lambda values:
pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_PlotK_CrossValidation_AND_CrossEntropy.pdf"), onefile = T)


#Lambda 0
plot(TESS3_Lambda_0, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 0")  

#Lambda 0.25
plot(TESS3_Lambda_0.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 0.25")  

#Lambda 0.5
plot(TESS3_Lambda_0.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 0.5")  

#Lambda 0.75
plot(TESS3_Lambda_0.75, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 0.75")  

#Lambda 1
plot(TESS3_Lambda_1, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 1")   

#Lambda 1.25
plot(TESS3_Lambda_1.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 1.25")  

#Lambda 1.5
plot(TESS3_Lambda_1.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score", main = "lambda = 1.5")   

dev.off()


#D. Choose best K:
optimal_K = 6

a <- min(TESS3_Lambda_0[[optimal_K]]$crossvalid.crossentropy) 
b <- min(TESS3_Lambda_0.25[[optimal_K]]$crossvalid.crossentropy) 
c <- min(TESS3_Lambda_0.75[[optimal_K]]$crossvalid.crossentropy) 
d <- min(TESS3_Lambda_0.5[[optimal_K]]$crossvalid.crossentropy)
e <- min(TESS3_Lambda_1[[optimal_K]]$crossvalid.crossentropy) 
f <- min(TESS3_Lambda_1.25[[optimal_K]]$crossvalid.crossentropy)  
g <- min(TESS3_Lambda_1.5[[optimal_K]]$crossvalid.crossentropy)

list <- c(a, b, c, d, e, f, g)
min <- min(list); min == list


# The best lambda value has lower cross entropy values
#Best = Lambda_1

best_lambda = TESS3_Lambda_1

# Plot
plot(best_lambda, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-valid & Cross-entropy")

## K = 6 (segunda quebra de tendência)

## Save plot
pdf(paste0("./1-Filtering/cha/Results_tess3/TESS3_CrossValidation_AND_CrossEntropy_BestK_filtering.pdf"), onefile = T)

plot(best_lambda, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-valid & Cross-entropy")

dev.off()
## FST OUTLIER DETECTION
## Use the best lambda project to retrieve results
res = Gettess3res(best_lambda, K=6)

#A. Compute the FST statistics using best run
FST = res$Fst

#B. Compute the GIF - genomic inflation factor
gif = res$gif; gif #  2.077798

#C. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(res$Q)[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/gif, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

#D. Test different lambda values and plot histogram of p-values
adj.p.values2 = pchisq(z.scores^2/1.05, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
adj.gif = 1.05
hist(adj.p.values2, col = "green")

## Test different FDR
C_fst_0.1 <- candidates(alpha=0.1, adj.p.values2); length(C_fst_0.1) 
ManPlot(adj.p.values2, C_fst_0.1, paste0("10% FDR - Removing ",  length(C_fst_0.1), " SNPs"))

C_fst_0.05 <- candidates(alpha=0.05, adj.p.values2); length(C_fst_0.05) 
ManPlot(adj.p.values2, C_fst_0.05, paste0("5% FDR - Removing ",  length(C_fst_0.05), " SNPs"))

## Save vcf for comparison (Venn Diagram)
## Subset and and save FST outlier for later comparison

snps_fil_ldF_candidate <- Subset(snps, sites = C_fst_0.1)
VCFsummary(snps_fil_ldF_candidate) # 93 individuals and 1066 SNPs.
Save(snps_fil_ldF_candidate, file="./1-Filtering/cha/FST-Outliers/candidates_fst_tess_0.1_FDR.vcf")

## Subset and and save FST outlier for later comparison
snps_fil_ldF_candidate <- Subset(snps, sites = C_fst_0.05)
VCFsummary(snps_fil_ldF_candidate) # 93 individuals and 802 SNPs.
Save(snps_fil_ldF_candidate, file="./1-Filtering/cha/FST-Outliers/candidates_fst_tess_0.05_FDR.vcf")

#E. Save Results
pdf("./1-Filtering/cha/Results_tess3/Charinus_adj_p_values_Manhattan.pdf", onefile = T)

hist(adj.p.values, col = "red", main=paste0("Genomic Inflation Factor = ", gif))

hist(adj.p.values2, col = "green", main=paste0("Adj. Genomic Inflation Factor = ", adj.gif))

ManPlot(adj.p.values2, C_fst_0.1, paste0("10% FDR - Removing ",  length(C_fst_0.1), " SNPs"))

ManPlot(adj.p.values2, C_fst_0.05, paste0("5% FDR - Removing ",  length(C_fst_0.05), " SNPs"))

dev.off()


## Choose one FDR
alpha=0.05

C_fst <- candidates(alpha, adj.p.values2); length(C_fst) # 802

ManPlot(adj.p.values2, C_fst, paste0("Fst Outlier Detected"," ", 
                                     length(C_fst)," ", "SNPs"))

###7.4. EXCLUDING FST OUTLIERS

## Subset and and save FST outlier for later comparison
snps_fil_ldF_candidate <- Subset(snps, sites = C_fst)

snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_neutral <- Subset(snps, sites=N_snp)
VCFsummary(snps_neutral) # 93 individuals and 2291 SNPs.
length(N_snp)
length(snps@site_id)-length(candidates) 
length(snps_neutral@site_id)

neutral_after_fst_tess = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_tess # 93 individuals and 2291 SNPs."

#B. Save vcf
Save(snps_neutral, paste0("vcf/cha/Neutral/tess/", project_name, "_filtered_neutral_TESS_FST-OUtlier.vcf"))
Save(snps_neutral, paste0("vcf/cha/Neutral/final/", project_name, "_filtered_neutral_TESS_final.vcf"))
Save(snps_neutral, paste0("vcf/cha/PopStructure/AllSamples/", project_name, "_filtered_neutral_TESS_final.vcf"))


###############################################################################
############# SELECT ONE SAMPLE BY CAVE BASED ON LOWER MISSING VALUES #########
###############################################################################

snps_fil_ldF_neutral <- vcfLink(paste0("vcf/cha/Neutral/final/", project_name, "_filtered_neutral_TESS_final.vcf"), overwriteID=T)
VCFsummary(snps_fil_ldF_neutral) ## 93 individuals and 2291 SNPs.

#Dataframe of Missing per individual
Missing_ind <- apply(GenotypeMatrix(snps_fil_ldF_neutral), 1, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing_ind) ##Max missing =  43.91 
hist(Missing_ind)
Missingind.df <- as.data.frame(Missing_ind)
Missingind.df$ID <- row.names(Missingind.df)
head(Missingind.df)

### Prepare subset with individuals with lower missing values
caves <- snps_fil_ldF_neutral@meta[ ,c(2,3,4,5)]
head(caves)

### Check IDs
identical(rownames(Missingind.df$ID), as.character(caves$sample_name)) ## check if the order of Sample are the same in both datasets

### Sort IDs and check
caves_sorted <- caves[order(match(caves[,1], Missingind.df$ID)),]
head(caves_sorted)
identical(as.character(caves_sorted$sample_name), as.character(Missingind.df$ID)) # TRUE

## Add caves
Missingind.df$caves <- caves$local_ID
head(Missingind.df)
length(unique(Missingind.df$caves)) ## 27 levels

## Plot to check
boxplot(Missingind.df$Missing_ind ~ Missingind.df$caves, ylab="Missing/ind (%)", xlab =NULL, las=2)

##keep one individual per cave with lower missing values
library(dplyr)
indtokeep <- Missingind.df %>% 
  group_by(caves) %>% 
  filter(Missing_ind==min(Missing_ind)) %>% 
  ungroup() %>%
  distinct(caves, .keep_all = TRUE)


head(indtokeep)

length(indtokeep$caves) # 27 
length(indtokeep$ID) # 27
boxplot(indtokeep$Missing_ind ~ indtokeep$caves, ylab="Missing/ind (%)", xlab =NULL, las=2)

### Subset
snps_caves <- Subset(snps_fil_ldF_neutral, samples=indtokeep$ID) 
snps_caves@meta

#Saves filtered SNPs
Save(snps_caves, paste0("vcf/cha/PopStructure/OneSampleByCave/", project_name, "_filtered_ld_hw_neutral_OneSampleByCave.vcf"))
VCFsummary(snps_caves) # 27 individuals and 2291 SNPs.



####################################################################################
###### FILTER ADAPTIVE LOCI by depth, quality, maf, snps with missing data #########
################ and removing high LD loci within same contig ######################
####################################################################################

project_name = "charinus"

###7. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE:
#A. from dataset with low missing by individuals
VCFsummary(snps_ind_pc) # 93 individuals and 17310 SNPs.

snps_fil_low = Filter(snps_ind_pc, filterOptions(minQ=30, maf=0.03, min.meanDP=20, max.meanDP=800)) 
VCFsummary(snps_fil_low) # 93 individuals and 15088 SNPs.

#B. from dataset with all individuals
#snps_fil_high <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.8, maf=0.05, min.meanDP=20, max.meanDP=200)) 
#VCFsummary(snps_fil_high)

#C. choose one dataset:
snps_fil = snps_fil_low
adaptativeFilter_SZ =  capture.output(VCFsummary(snps_fil))

###8. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:

# Define R² value:
r2 = 0.8

#A. remove snps with R² value
ld_within <- Linkage(snps_fil, type="geno-r2", linkageOptions(min.r2= r2))
head(ld_within)
hist(ld_within$R.2)
save(ld_within, file = "1-Filtering/cha/ld/ld_within.RData")
load("1-Filtering/cha/ld/ld_within.RData")

#B. Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within$ID1
nold_snps <- snps_fil@site_id[!(snps_fil@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil, sites=nold_snps) # Keep snps that are not in LD.
adaptativeLD_SZ = capture.output(VCFsummary(snps_fil_ld))
adaptativeLD_SZ # 93 individuals and 10468 SNPs.

###9. VERIFY THE REAL MISSING DATA VALUE IN THE DATASET:
#A. Missing per individual:
Missing_ind <- apply(GenotypeMatrix(snps_fil_ld),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 44.908%
hist(Missing_ind)

#B. Missing per locus:
Missing <- apply(GenotypeMatrix(snps_fil_ld), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 30.108%
hist(Missing)
VCFsummary(snps_fil_ld) ## 93 individuals and 10468 SNPs.

###10. SAVE THE .VCF FILE FILTERED WITH ADAPTATIVE SNPS
Save(snps_fil_ld, paste0("vcf/cha/Adaptive/", project_name,"_filtered_within_ld_adaptive.vcf"))


###16. SAVE THE FILTERING RESULTS BY DATASETS:
#A. Combine all results
datasets = c(raw_SZ, filtered_SZ, removed_missing_SZ, 
             adaptativeFilter_SZ, 
             adaptativeLD_SZ, 
             neutralFilter_SZ, 
             neutralLDWithin_SZ, 
             neutralLDBetween_SZ, 
             neutralFilter_HWE_SZ2,
             neutral_after_fst_tess)


#B. Create a matrix to save the results. Let the last line for FST outliers filtering step.
results_snps = matrix("NA", 12, 2)
colnames(results_snps) = c("Individuals", "SNPs")
rownames(results_snps) = c("Raw", "Biallelic SNPs", "Low Missing Data per Individual", 
                           "Quality Filter for Adaptive SNPs", 
                           "Linkage Desequilibrium Filter for Adaptive SNPs (Within Contigs)", 
                           "Quality Filter for Neutral SNPs", 
                           "Linkage Desequilibrium Filter for Neutral SNPs (Within Contigs)", 
                           "Linkage Desequilibrium Filter for Neutral SNPs (Between Contigs)", 
                           "Hardy-Weinberg Equilirbirum Filter for Neutral SNPs (HWE)",
                           "Neutral SNPs Tess (FST outliers)")

#C. Save the results in matrix
for (i in 1:length(datasets)){
  text = scan(text = datasets[i], what = "")
  results_snps[i,] = c(text[1], text[4])
}

#D. Verify the results
as.data.frame(results_snps)

#E. Save as .csv
write.csv(as.data.frame(results_snps), paste0("./1-Filtering/cha/Results_Metafiles/Results_SNPs_datasets_", project_name, ".csv"))

## END

