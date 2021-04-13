# Landscape genomics of cave invertebrates

These suit of scripts were developed to perform the following analysis:

## A) Population structure

1. SNPs filtering
 1.1. Filter SNPs dataset
 1.2. Filter highly related individuals

2. Assessment of population structure
 2.1. Population structure assessment using sNMF
 2.2. Population structure assessment using DAPC

3. Assessment of genetic diversity
 3.1. Estimate Ho, He, pi, F and Tajima's D
 3.2. Prepare file to estimate effective population size (Ne)

4. Assessment of fine-scale spatial genetic structure


## B) Landscape genomic analysis

- Caculate environmental distances
- Prepare rasters

1. Optimize resistance distances using resistanceGA

2. Reclassify rasters

3. Assessment of gradient surface metrics
 3.1. Clip landscapes (a buffer including two locations)
 3.2. Calculate gradient surface metrics

4. Calculate resistance distances

5. Calculate topographic distances

6. Bind final data frame for analysis
 6.1. Bind surface metrics by caves
 6.2. Bind all variables

7. Pre-selection of variables

8. Fit MLPE models  

