###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

## THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3



## GOAL: DETERMINE RELATEDNESS CUTOFF AND FILTER OUT HIGHLY RELATED INDIVIDUALS ----

## REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

## Load packages ----
library(r2vcftools)
library(anchors)
library(dplyr)
library(ggplot2)
library(gridExtra)

## Load source
source("source/r2vcftools_utils_2019.R")

## Set project name
project_name = "charinus"

## Set paths
vcf_file <- "./vcf/cha/PopStructure/charinus_filtered_neutral_TESS_final.vcf"

# Load VCF ----
neutralsnps <- vcfLink(vcf_file, overwriteID = TRUE) 
VCFsummary(neutralsnps) # 93 individuals and 2291 SNPs.
names(neutralsnps@meta)
#genotypes <- GenotypeMatrix(neutralsnps)
#genotypes[1:10, 1:10]

#write.table(genotypes, "./vcf/cha/PopStructure/charinus_filtered_neutral_TESS_final_genotypematrix.txt",
#            append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)

#### 1. Calculate Relatedness (Yang) ----
REL <- Relatedness(neutralsnps, type = "yang", verbose = TRUE)
length(unique(REL$INDV1)) ## 93
length(unique(REL$INDV2)) ## 93
max(REL$RELATEDNESS_AJK) ## 1.4471

## Check self relatadness
REL_self <- REL %>% filter(INDV1 == INDV2)
head(REL_self)
tail(REL_self)

nrow(REL_self) # 93
hist(REL_self$RELATEDNESS_AJK, main = "Self-relatedness \n Charinus")
min(REL_self$RELATEDNESS_AJK) ## 0.758738
max(REL_self$RELATEDNESS_AJK) ## 1.4471


## Filter out self relatedness
REL <- REL %>% filter(INDV1 != INDV2)
nrow(REL) # 4278
colnames(REL) <- c("ID1","ID2","RELATEDNESS_AJK_Yang")

min(REL$RELATEDNESS_AJK_Yang) ## -0.160505
max(REL$RELATEDNESS_AJK_Yang) ##  0.560098

length(unique(REL$ID1)) ## 92
length(unique(REL$ID2)) ## 92

## ID differences between columns before filtering
dif_ID1_ID2 <- setdiff(unique(REL$ID1), unique(REL$ID2)); dif_ID1_ID2 ## "CH_167_sorted"
dif_ID2_ID1 <- setdiff(unique(REL$ID2), unique(REL$ID1)); dif_ID2_ID1 ## "CH_270_sorted"


## Add cave IDS using the inner_join function from dplyr
head(neutralsnps@meta)
meta <- neutralsnps@meta %>% dplyr::select(sample_name, local_ID) # make new df with sample and cave ids
REL1 <- REL %>% left_join(meta, by=c("ID1" = "sample_name")) %>% # merge dfs by sample name (ID1)
  left_join(meta, by=c("ID2" = "sample_name")) %>% # merge dfs by sample name (ID2)
  dplyr::rename(Cave1 = "local_ID.x",
                Cave2 = "local_ID.y")

## Individuals from the same cave
RELw <- REL1 %>% 
        filter(Cave1 == Cave2) %>% 
        mutate(Group = "Within caves")

## Individuals from different caves
RELb <- REL1 %>% 
        filter(Cave1 != Cave2) %>% 
        mutate(Group = "Between caves")

## Bind information and check later
REL_full <- bind_rows(RELw, RELb)

head(REL_full)
nrow(REL_full) ## 4278


## Check distribution of full dataset

REL_full %>% ggplot(aes(x=RELATEDNESS_AJK_Yang)) + geom_histogram(bins=100) + 
  xlab("Mean relatedness (full dataset)")


## Histogram
h1 <- REL_full %>% ggplot(aes(x=RELATEDNESS_AJK_Yang)) + 
  geom_histogram(bins=100) + 
  xlab("Genetic relatedness") +
  ylab("Number of pairs of samples") +
  annotate("text", x=0.4, y=200, label= "Full dataset \n N caves = 27", size = 6, col="blue") + 
  annotate("text", x = 0.6, y = 320, label = "", size=10) +
  theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
        axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
  theme(axis.text = element_text(size = 18, face = "bold", color = "black"))

h1

## Save
tiff(filename = "./1-Filtering/cha/Results_cutoff/relatedness_distribution_FULL.TIF",
     units="in", width=12, height=10, 
     compression = "lzw", bg = "white", res = 300)
h1

dev.off()


## Boxplot
boxpl1 <- REL_full %>% 
          ggplot(aes(x=Group, y=RELATEDNESS_AJK_Yang)) + 
          geom_boxplot() + 
          ylab("Genetic relatedness") +
          theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
          axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
          theme(axis.text = element_text(size = 18, face = "bold", color = "black"))


## Save
tiff(filename = "./1-Filtering/cha/Results_cutoff/boxplot_relatedness_distribution_FULL.TIF",
     units="in", width=12, height=10, 
     compression = "lzw", bg = "white", res = 300)

boxpl1

dev.off()



## Check maximum values within and between caves before filtering

## Within caves
names(REL_full)
REL_w <- REL_full %>% filter(Group == "Within caves")
max(REL_w$RELATEDNESS_AJK_Yang) ## 0.560098

## Between caves
REL_b <- REL_full %>% filter(Group == "Between caves")
max(REL_b$RELATEDNESS_AJK_Yang) ## 0.480254


## Determine threshold value = lowest value in the dataset ----
value = min(REL_full$RELATEDNESS_AJK_Yang)
value ## -0.160505


## Filter by range and check subset (above threshold = high related individuals)
REL_fil2 <- REL_full %>% filter(RELATEDNESS_AJK_Yang > abs(value))
nrow(REL_fil2) ## 131

REL_fil2 %>% ggplot(aes(x=RELATEDNESS_AJK_Yang)) + geom_histogram(bins=100) +
  xlab("Genetic relatedness") +
  annotate("text", x=0.4, y=5, label= "Above threshold \n 0.1605", size = 5, col="blue") + theme_minimal()

max(REL_fil2$RELATEDNESS_AJK_Yang) ## 0.560098
head(REL_fil2)
tail(REL_fil2)


## Filter by range and check subset (above 2x threshold = high related individuals)
REL_fil3 <- REL_full %>% filter(RELATEDNESS_AJK_Yang > abs(2*value))
nrow(REL_fil3) ## 23

min(REL_fil3$RELATEDNESS_AJK_Yang) ## 0.3212
max(REL_fil3$RELATEDNESS_AJK_Yang) ## 0.560098

head(REL_fil3)
tail(REL_fil3)

sort(unique(REL_fil3$Cave1))
sort(unique(REL_fil3$Cave2))


## 2. Drop samples from full dataset ----
drop_samples <- REL_fil3$ID1
unrel_samples <- REL_full[!(REL_full$ID1 %in% drop_samples), ] ## all "unrelated" samples (low related individuals)
names(unrel_samples) 
nrow(unrel_samples) ##  3531
hist(unrel_samples$RELATEDNESS_AJK_Yang)

min(unrel_samples$RELATEDNESS_AJK_Yang) ## -0.160505
max(unrel_samples$RELATEDNESS_AJK_Yang) ##  0.315019

length(unique(unrel_samples$ID1)) ## 81
length(unique(unrel_samples$ID2)) ## 92

dif_ID1_ID2 <- setdiff(unique(unrel_samples$ID1), unique(unrel_samples$ID2)); dif_ID1_ID2 ## 
dif_ID2_ID1 <- setdiff(unique(unrel_samples$ID2), unique(unrel_samples$ID1)); dif_ID2_ID1 ## 

length(unique(unrel_samples$Cave1)) ## 26
length(unique(unrel_samples$Cave2)) ## 27


## Save histogram
h3 <- unrel_samples %>% ggplot(aes(x=RELATEDNESS_AJK_Yang)) + geom_histogram(bins=100) +
  xlab("Mean relatedness (filtered dataset)") +
  annotate("text", x=0.2, y=150, label= "Threshold = 0.3150 \n N caves = 26", size = 6, col="blue") + 
  xlab("Genetic relatedness") +
  ylab("Number of pairs of samples") +
  theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
        axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
  theme(axis.text = element_text(size = 18, face = "bold", color = "black"))

h3

tiff(filename = "./1-Filtering/cha/Results_cutoff/relatedness_distribution_below_2x_threshold.TIF",
     units="in", width=12, height=10, 
     compression = "lzw", bg = "white", res = 300)

h3

dev.off()


## Histogram (supplementary material) ----
h_supp_mat <- REL_full %>% ggplot(aes(x=RELATEDNESS_AJK_Yang)) + 
  geom_histogram(bins=100) + 
  xlab("Genetic relatedness") +
  ylab("Number of pairs of samples") +
  annotate("text", x=0.19, y=200, label= paste0("Cutoff = ", abs(2*value)), size = 6, col="red") + 
  annotate("text", x = 0.6, y = 320, label = "", size=10) +
  geom_vline(xintercept = abs(2*value), linetype="dashed", color="red", size = 1.5) +
  theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
        axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
  theme(axis.text = element_text(size = 18, face = "bold", color = "black"))

h_supp_mat

## Save
tiff(filename = "./1-Filtering/cha/Results_cutoff/relatedness_distribution_FULL_CUTOFF_SUPP_MATERIAL.TIF",
     units="in", width=7, height=5, 
     compression = "lzw", bg = "white", res = 300)
h_supp_mat

dev.off()

## Check differences in caves
dif_Cave_1_2 <- setdiff(unique(unrel_samples$Cave1), 
                        unique(unrel_samples$Cave2)); dif_Cave_1_2 ## 0

dif_Cave_2_1 <- setdiff(unique(unrel_samples$Cave2), 
                        unique(unrel_samples$Cave1)); dif_Cave_2_1 ## S11D_0078


## Rename set of unrelated samples
REL_fil4 = unrel_samples

## Count individuals within each cave
counts <- REL_fil4 %>% group_by(Cave1) %>%
  summarise(n_samples = length(unique(ID1))) %>%
  arrange(n_samples)


## 3. Create new dataset grouping within and between cave relatedness ----
names(REL_fil4)
max(REL_fil4$RELATEDNESS_AJK_Yang) ## 0.315019


## Plot relatedness dsitributions

h4 <- REL_fil4 %>% filter(Group=="Within caves") %>% 
      ggplot(aes(x=RELATEDNESS_AJK_Yang)) + geom_histogram() + 
      xlab("Genetic relatedness \n within caves")  +
      ylab("Number of pairs of samples") +
      theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
      axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
      theme(axis.text = element_text(size = 18, face = "bold", color = "black"))

h4

h5 <- REL_fil4 %>% filter(Group=="Between caves") %>% 
      ggplot(aes(x=RELATEDNESS_AJK_Yang)) + geom_histogram() + 
      xlab ("Genetic relatedness \n between caves")  +
      ylab("Number of pairs of samples") +
      theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
      axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
      theme(axis.text = element_text(size = 18, face = "bold", color = "black"))

h5

## 4. Test differences between groups (parametric and non-parametric tests) ---- 

t.test(RELATEDNESS_AJK_Yang ~ Group, data=REL_fil4) # parametric test
wilcox.test(RELATEDNESS_AJK_Yang ~ Group, data=REL_fil4) # non-parametric test

boxpl3 <- REL_fil4 %>% 
          ggplot(aes(x=Group, y=RELATEDNESS_AJK_Yang)) + 
          geom_boxplot() + 
          ylab("Genetic relatedness") +
          theme(axis.title.y = element_text(size=22, color = "black", face = "bold"),
          axis.title.x = element_text(size=22, color = "black", face = "bold")) + 
          theme(axis.text = element_text(size = 18, face = "bold", color = "black"))


boxpl3


## Save panel ----

grid.arrange(h1, boxpl1, h3, boxpl3, ncol=2)


tiff(filename = "./1-Filtering/cha/Results_cutoff/relatedness_PANEL.tif",
     width = 4000, height = 3000, units = "px", 
     compression = "lzw", bg = "white", res = 300)

grid.arrange(h1, boxpl1, h3,boxpl3, ncol=2)


dev.off()


## 5. Count individuals within each cave in filtered dataset ----
counts <- REL_fil4 %>% group_by(Cave1) %>%
  summarise(n_samples = length(unique(ID1))) %>%
  arrange(n_samples)

## View full tables
View(counts)

## 6. Check sample ID differences between raw and filtered dataset ----
dif1 <- setdiff(unique(REL$ID1), unique(REL_fil4$ID1)); dif1
dif2 <- setdiff(unique(REL$ID2), unique(REL_fil4$ID2)); dif2


## Check caves 
neutralsnps@meta[neutralsnps@meta$sample_name %in% dif1, ] 
neutralsnps@meta[neutralsnps@meta$sample_name %in% dif2, ]  


## 7. Filter droped sample IDs ----

## Check sample names
sort(neutralsnps@meta$sample_name)

## Use this sample IDs to filter highly related individuals from the VCF
sample_IDs <- neutralsnps@meta$sample_name[!(neutralsnps@meta$sample_name %in%
                                               dif1)]
sort(sample_IDs)

## Check caves
cave_IDs_1 <- unique(REL_fil4$Cave1); length(cave_IDs_1) ## 26
cave_IDs_2 <- unique(REL_fil4$Cave2); length(cave_IDs_2) ## 27

### 8. Subset and save filtered vcf ----
sub_snps <- Subset(neutralsnps, samples = sample_IDs)
VCFsummary(sub_snps) # 82 individuals and 2291 SNPs
Save(sub_snps, file= paste0("./vcf/cha/Neutral/rel_cutoff/", project_name, "_CutOffRelatedness_neutral_snps.vcf"))


## 9. Check relatedness in subset  ----
VCFsummary(sub_snps) # 82 individuals and 2291 SNPs. ## this dataset will be used to test population structure

REL <- Relatedness(sub_snps, type = "yang", verbose = TRUE)
max(REL$RELATEDNESS_AJK) ## 1.84376

colnames(REL) <- c("ID1","ID2","RELATEDNESS_AJK_Yang")
fil <- REL %>% filter(ID1 !=  ID2)
head(fil)
nrow(fil) ## 3321
hist(fil$RELATEDNESS_AJK)
max(fil$RELATEDNESS_AJK) ## 0.523807




## 10. Prepare inputs for Landscape Genetic Analysis ----
## Here we will use the full dataset, because highly related individuals among caves will be modelled anyway

# Load VCF
neutralsnps <- vcfLink(vcf_file, overwriteID = TRUE) 
VCFsummary(neutralsnps) # 93 individuals and 2291 SNPs.
names(neutralsnps@meta)

## A. Save sample name, location and coordinates ----
write.csv(neutralsnps@meta[ , 2:5], "./coords/cha/coords_charinus_93ind.csv", row.names = F)

#### B. Relatedness dataframe ----
REL <- Relatedness(neutralsnps, type = "yang", verbose = TRUE)
colnames(REL) <- c("ID1","ID2","RELATEDNESS_AJK_Yang")
head(REL)
length(unique(REL$ID1)) ## 93
length(unique(REL$ID2)) ## 93
max(REL$RELATEDNESS_AJK) ## 1.4471

## Add cave IDS using the inner_join function from dplyr
head(neutralsnps@meta)
meta <- neutralsnps@meta %>% dplyr::select(sample_name, local_ID) # make new df with sample and cave ids
REL_IDs <- REL %>% left_join(meta, by=c("ID1" = "sample_name")) %>% # merge dfs by sample name (ID1)
  left_join(meta, by=c("ID2" = "sample_name")) %>% # merge dfs by sample name (ID2)
  dplyr::rename(Cave1 = "local_ID.x",
                Cave2 = "local_ID.y")
head(REL_IDs)

## Save
write.csv(REL_IDs, "./5-LandGenAnalysis/results/gen_dist/cha/charinus_relatedness_93samples.csv", row.names = F)


## C. Relatedness matrix ----
nb_cha <- length(unique(REL_IDs$ID1))
nb_cha

rb_cha <- matrix(0, nb_cha, nb_cha)
rb_cha[lower.tri(rb_cha, diag=T)] <- REL_IDs$RELATEDNESS_AJK
rb_cha <- rb_cha + t(rb_cha) # to make symmetric
diag(rb_cha) <- diag(rb_cha)/2 # to correct diagonal for previous line
class(rb_cha)
dim(rb_cha)

## Save genetic distance matrix to run resistanceGA
write.csv(rb_cha, "./5-LandGenAnalysis/results/gen_dist/cha/charinus_GenDistMatrix_93samples.csv", row.names = F)


## D. Relatedness of One Sample By Cave ----

#A. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/cha/PopStructure/OneSampleByCave/", project_name, "_filtered_ld_hw_neutral_OneSampleByCave.vcf"), overwriteID=T)
VCFsummary(snps_neutral) # 27 individuals and 2291 SNPs.
head(snps_neutral@meta[, 2:5])

## Save coords
write.csv(snps_neutral@meta[, 2:5], "./coords/cha/coords_charinus_27ind_OneSampleByCave.csv", row.names = F)

## Relatedness dataframe
REL <- Relatedness(snps_neutral, type = "yang", verbose = TRUE)
colnames(REL) <- c("ID1","ID2","RELATEDNESS_AJK_Yang")
nrow(REL)
head(REL)
hist(REL$RELATEDNESS_AJK)

## Append cave IDs
head(snps_neutral@meta)
meta <- snps_neutral@meta %>% dplyr::select(sample_name, local_ID) # make new df with sample and cave ids
REL_sub <- REL %>% left_join(meta, by=c("ID1" = "sample_name")) %>% # merge dfs by sample name (ID1)
           left_join(meta, by=c("ID2" = "sample_name")) %>% # merge dfs by sample name (ID2)
           dplyr::rename(Cave1 = "local_ID.x",
                         Cave2 = "local_ID.y")
  
nrow(REL_sub) #378
head(REL_sub)

## Save
write.csv(REL_sub, "./5-LandGenAnalysis/results/gen_dist/cha/charinus_relatedness_27samples_OneSampleByCave.csv", row.names = F)


## END OF SCRIPT ----

