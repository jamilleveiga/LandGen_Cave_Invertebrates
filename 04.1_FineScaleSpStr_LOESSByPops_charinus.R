###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################

##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


## FINE SCALE SPATIAL STRUCTURE ----

## REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS TUTORIAL IS ORGANIZED IN DIFFERENT "Steps" AND INSIDE EACH STEP THERE ARE "Actions" DENOTED BY NUMBERS (#1) AND INSIDE EACH ACTION COULD EXIST "Operations" INDICATED BY LETTES (#A.)
#D. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


## INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline


## GOALS FOR THIS STEP:
#A. PERFORM FINE-SCALE SPATIAL GENETIC STRUCTURE BASED ON:
    #A1. LOESS MODEL
    #A2. sPCA MODEL

## CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

## REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

## LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("source/r2vcftools_utils_2019.R")

## Load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(ggplot2, r2vcftools, geosphere, raster, rgdal, gridExtra, dplyr, LEA,  
               vcfR, dartR, adegenet, splancs, usdm, adespatial, spdep, phuse, remotes, scales)


## CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
dir_names = (c("4-FineScaleSpatialStructure/LOESS/cha/Results_LOESS/POP1/",
               "4-FineScaleSpatialStructure/LOESS/cha/Results_LOESS/POP2/",
               "4-FineScaleSpatialStructure/LOESS/cha/temp/POP1/",
               "4-FineScaleSpatialStructure/LOESS/cha/temp/POP2/",
               "4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/",
               "4-FineScaleSpatialStructure/sPCA/cha/temp/"))

for (i in 1:length(dir_names)){
  if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
  else (print (paste0(dir_names[i], " has already been created. Be careful with overwritting")))}


## FIT LOESS ----

###1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "charinus"

###2. LOAD VCF FILES AND GEOGRAPHICAL INFORMATIONS: 
#A. Load neutral .vcf file with geographical information:
snps_fil_ldF_neutral = vcfLink(paste0("./vcf/cha/PopStructure/tess/", project_name, "_filtered_neutral_POP_TESS.vcf")
                               , overwriteID=T)

VCFsummary(snps_fil_ldF_neutral) # 82 individuals and 2291 SNPs.

REL <- Relatedness(snps_fil_ldF_neutral, type = "yang", verbose = TRUE)
nrow(REL)
head(REL)
hist(REL$RELATEDNESS_AJK)

###Transform Relatedness dataframe into matrix
nb_cha <- length(unique(REL$INDV1))
nb_cha

rb_cha <- matrix(0, nb_cha, nb_cha)
rb_cha[lower.tri(rb_cha, diag=T)] <- REL$RELATEDNESS_AJK
rb_cha <- rb_cha + t(rb_cha) # to make symmetric
diag(rb_cha) <- diag(rb_cha)/2 # to correct diagonal for previous line
rb_cha


## Geographic distances
coord <- snps_fil_ldF_neutral@meta[,4:5]
head(coord)
coordinates(coord) <- snps_fil_ldF_neutral@meta[,4:5]
projection(coord) <- crs("+init=epsg:32722") #define utm datum only for coord cols

#transform to decimal degrees
coord_cha <- spTransform(coord, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"))
projection(coord_cha)

sdist_cha  <- distm(coord_cha, fun=distGeo)
sdist_cha  <- sdist_cha /1000
nrow(sdist_cha)
ncol(sdist_cha)
max(sdist_cha)
dim(sdist_cha)
sdist_cha

ot <- mDistoLoess(sdist_cha, rb_cha, 999)
save(ot, file="./4-FineScaleSpatialStructure/LOESS/cha/temp/ot.RData")
load(file="./4-FineScaleSpatialStructure/LOESS/cha/temp/ot.RData")

dd <- data.frame(cov = sdist_cha[lower.tri(sdist_cha)], obs = attr(ot$disto, "obs"), lci = attr(ot$disto, "CI")[1,], uci = attr(ot$disto, "CI")[2,])
save(dd, file="./4-FineScaleSpatialStructure/LOESS/cha/temp/dd.RData")
load(file="./4-FineScaleSpatialStructure/LOESS/cha/temp/dd.RData")

### Plot figure
LP3_cha <- ggplot(dd, aes(x=cov,y=obs)) + geom_line(size=1.2) + 
  geom_ribbon(aes(ymax = uci, ymin = lci), alpha = 0.5) + 
  theme_bw() + 
  geom_hline(yintercept = mean(ot$disto), lty = 3) + 
  ylab("Mean relatedness (r) \npredicted by LOESS") + xlab("Distance (Km)") + 
  geom_rug(sides = "b", alpha = 0.02) +
  annotate("text", x = 8, y = 0.005, label = "6.1 Km", size=7, color = "red") +
  geom_vline(xintercept = 6.1, linetype="dashed", color="red") +
  theme(axis.title.y = element_text(size=20, color = "black", face = "bold"),
        axis.title.x = element_text(size=20, color = "black", face = "bold")) + 
  theme(axis.text = element_text(face = "bold", color = "black", size = 20))

LP3_cha


# Save figure in pdf
pdf(paste0("./4-FineScaleSpatialStructure/LOESS/cha/Results_LOESS/Loess_",
    project_name,"_red_line_CI.pdf"), 
    width = 7, height = 5)
LP3_cha
dev.off()



## sPCA ----

## Fine-scale population strucutre using sPCA to compare the allelic frequency of individuals to allelic frequencies of their neighbours
###7. LOAD INPUTS

#A1. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("./vcf/cha/PopStructure/tess/", project_name, "_filtered_neutral_POP_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) # 82 individuals and 2291 SNPs.

#A2. Create a Genind object based on neutral VCFR
vcf_neutral = read.vcfR(paste0("./vcf/cha/PopStructure/tess/", project_name, "_filtered_neutral_POP_TESS.vcf"), verbose = FALSE)
vcf_neutral # 82 samples; 1297 CHROMs; 2,291 variants

#B. Convert to GENELIGHT
gl = vcfR2genlight(vcf_neutral)

#D. Convert to GENEIND
genind = gl2gi(gl)
genind ## 82 individuals; 2,291 loci; 4,576 alleles; size: 2.6 Mb
  
#E. Create a PCA object
input_scaled = scaleGen (genind, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = 93)

#F. Geographic data, removing equal coordinates
coordinates = snps_neutral@meta[,4:5]

#Here we will run sPCA following the tutorial provided by Jombart.
#  http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf

#You can test two different connection networks:
#Delaunay triangulation (type=1) = more suited to uniform sampling. Other softwares use this one, for comparasion it is a good one.

#Gabriel graph (type=2). The Gabriel graph is a subgraph of the Delaunay triangulation. It can be found in linear time if the Delaunay triangulation is given (Matula & Sokal 1980). The Gabriel graph contains, as subgraphs, the Euclidean minimum spanning tree (type=4), the relative neighborhood graph (type=3), and the nearest neighbor graph (type=6). 

#Relative neighbours (type=3)

#Minimum spanning tree (type=4)

#Neighbourhood by distance (type=5) =  You need to inform minimum (d0) and maximum distance (d1) IN NUMBER OF INDIVIDUALS 0 AND 2 IS A PATTERN. Biologically-defined distance would be better in this case (home-range, dispersion for plants)

#K nearests neighbours (type=6) = You need to inform the number of neighbours per point (argument k).Shoul be lower than one-third of the of the number of data points

#Inverse distances (type=7) = recommended when clustering of samples is present. This is not a true neighbouring graph: all sites are neighbours, but the spatial weights are directly proportional to the inversed spatial distances. You need to inform the minimum distance between any two distinct points (dmin) and the exponent of the inverse distance matrix (a).

#G. Connection network - CN1-4 doesn't allow duplicated locations
CN5 = chooseCN(coordinates, type=5, d1 = 0, d2 = "dmin", plot=T)

#H. Then convert your CN to a listw using the nb2listw function.
cn = nb2listw(CN5) ## OK 

## RUN sPCA ----

# A. Run:
mySpca = multispati(pca_input, cn, nfposi = 10, nfnega = 10, scannf = F)
save(mySpca, file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/mySpca", project_name, ".RData"))
load(file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/mySpca", project_name, ".RData"))

# check number of PC in global and local structure
mySpca
plot(mySpca)

#check lag vector onto the principal axis dataframe with n rows and (nfposi + nfnega) columns
head(mySpca$li)

###9. TEST THE LOCAL AND GLOBAL STRUCTURE
# Lollipop plot interpretation: significant if lollipop is far from distribution of permutation values, not significant if lollipop falls along the distribution of permutation values.
# GLOBAL TEST: neighbouring individuals are more similar than expected
# LOCAL TEST: neighbouring individuals are more dissimilar than expected

#change the cn according to the sPCA chosen
myGtest = global.rtest(pca_input$tab, cn , nperm=9999)
myGtest
plot(myGtest) ## significant
save(myGtest, file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/GlobalTest", project_name, ".RData"))
load(file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/GlobalTest", project_name, ".RData"))


##change the cn according to the sPCA chosen. If no one Axis was retain you do need run this step
myLtest = local.rtest(pca_input$tab, cn, nperm=9999)
myLtest
plot(myLtest) ## significant
save(myLtest, file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/LocalTest", project_name, ".RData"))
load(file = paste0("./4-FineScaleSpatialStructure/sPCA/cha/temp/LocalTest", project_name, ".RData"))


###10. EXTRACT RESULTS TO DO THE INTERPOLATON ON QGIS
#verify lag vector onto the principal axesdata frame with n rows and (nfposi + nfnega) columns
head(mySpca$ls)
head(cbind(coordinates, mySpca$ls))

write.csv(cbind(coordinates, mySpca$ls),
paste0("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/", 
       project_name , "_sPCA_Axis_CN5.csv"))
  
###11. SAVE RESULTS
#A. Tables
results_sPCA = summary(mySpca)
write.csv(results_sPCA, paste0("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/POP1/", project_name, "_sPCA_Summary_CN5.csv"))

#B. Spatial and Inertial Components:
pdf("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Spatial_Inertia_CN5.pdf", onefile=F)
plot.new()
fgraph(mySpca)
dev.off()


#C. Barplot with Global and Local Structure in spectral color pattern
pdf("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Global_Local_Structure_Spectral_CN5.pdf", onefile =F)


barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="black")
dev.off()

#D. Barplot with Global and Local Structure in red color for PC significants

## Check number of PCs
length(mySpca$eig) ## 35

#Global Significant PCs:
n_posi = 2
n_neg = 1

#For Local PCs test diferent values to put color on graph, for no one PC let ir on 0. Usually you need to plus number of positive PC and Negative to put correctly colors:
n_eig = length(mySpca$eig)-3

# Adjust colors:
barplot(mySpca$eig, main="Eigenvalues of sPCA", 
        col= rep(c("red","grey", "blue"), c(n_posi, n_eig, n_neg)))

pdf("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Global_Local_Structure_Significant_CN5.pdf", onefile =F)
# Adjust colors:
barplot(mySpca$eig, main="Eigenvalues of sPCA", 
        col= rep(c("red","grey", "blue"), c(n_posi, n_eig, n_neg)))
dev.off()


#E. Colorplot of mySpca. It represents a cloud of points with colors corresponding to a combination of 
# 1,2 or 3 quantitative variables, assigned to RGB (Red, Green, Blue) channels. 
# For instance, this can be useful to represent up to 3 principal components in space.
#define the number of variables UP TO THREE (Positive + Negatives PCs)
variables = 3

pdf("4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Colorplot_Significant_PC_CN5.pdf", onefile =F)
data = as.data.frame(mySpca$ls[ , c(1,2,3)])
head(data)

colorplot(coordinates, 
          data, cex=3, 
          main="Colorplot of sPCA \n Neutral SNPs")
dev.off()



#E. Global Test of mySpca
pdf("./4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Global_Test_Significant_CN5.pdf", onefile =F)
plot(myGtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()

#E. Local Test of mySpca
pdf("./4-FineScaleSpatialStructure/sPCA/cha/Results_sPCA/Local_Test_Significant_CN5.pdf", onefile =F)
plot(myLtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()


## Make Interpolated Maps using QGIS 3.10     

#1. Install "Processing" plug-in in QGis
#2. Go to "Add a delimited text file" and open the file with spca axes and coordinates
#3. Save the opened file as a shapefile
#4. Interpolate each axis separately using the Interpolation tool in the Processing Toolbox:
#  Processing Toolbox > Interpolation > Interpolation IDW.
#5. At the Interpolation window: 
# a) select the shapefile with spca axes; 
# b) add the axis in the second box as a vertorial entry layer;
# c) select the "extension";
# d) Adjust the number of columns and rows;
# e) run the interpolation (the raster will appear in the main window);
#6. Repeat the previous steps to each axis
#7. Combine the three interpolated axis raster in a rgb raster using Raster > Micellaneous > Mosaic:
# a) Select the three interpolated axis raster;
# b) Select the option "put each file in a separate band";
# c) Run to build your mosaic.


## END OF SCRIPT ----
