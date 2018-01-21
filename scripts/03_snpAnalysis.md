---
title: "SNP Analysis"
author: "Steve Pederson"
date: "22 January, 2018"
output: 
  html_document: 
    fig_caption: yes
    fig_height: 8
    fig_width: 10
    keep_md: yes
    toc: yes
    toc_depth: 2
---




```r
library(tidyverse)
library(parallel)
library(pander)
library(scales)
library(reshape2)
library(readxl)
library(magrittr)
library(UpSetR)
library(grid)
library(qqman)
library(sp)
library(ggmap)
library(rgdal)
```



```r
nCores <- min(12, detectCores() - 1)
logit <- binomial()$linkfun
alpha <- 0.05
theme_set(theme_bw())
```


# Introduction

This analysis takes the filtered SNPs from the [previous section](02_snpFiltering) and analyses their population frequencies amongst the two populations.
The two populations are referred to below as Population `1` (*Gum Creek*, 1996) and Population `2` (*Oraparinna*, 2012).

## Outline of analysis

The analysis uses two simple models fully described by Lewis in [Genetic association studies: Design, analysis and interpretation](https://doi.org/10.1093/bib/3.2.146) to detect statistically significant differences between the two populations

1. A *full genotype* model based on **genotype counts**, which will detect changes in the structure of heterozygous or homozygous nature of genotypes
2. A *multiplicative model* based on **allele counts** in which the impact of a single copy of an allele will be tested.

- Prior to analysis, PCA was performed as a data exploration procedure
- Models were tested using Fisher's Exact Test to allow for low counts in contingency tables
- Multiple testing adjustments were performed using both *Bonferroni's method* and the *Benjamini-Hochberg method* to assess results in the context of both FWER and FDR control.

# Data Setup


```r
allData <- file.path("..", "data", "filteredSNPs.tsv.gz") %>%
  gzfile() %>%
  read_delim(delim = "\t") 
```


```r
popSizes <- allData %>%
  group_by(`Pop ID`) %>%
  summarise(N = max(N)) %>%
  mutate(Population = c("1996", "2012")[`Pop ID`],
         Description = c("Gum Creek", "Oraparinna")[`Pop ID`]) %>%
  ungroup() %>%
  dplyr::select(contains("Pop"), Description, N)
```

Data from the 20,336 SNPs selected for testing was loaded, and sample sizes formed as a simple table.
Using the observed major allele frequencies (`P`) and heterozygote frequencies (`Obs Het`), the allele frequency and genotype frequency tables were constructed.


```r
allCounts <- allData %>%
  mutate(A = round(P*2*N, 0), 
         B = round((1-P)*2*N, 0), 
         AB = round(`Obs Het`*N, 0), 
         AA = round((A - AB)/2, 0), 
         BB = round((B - AB)/2, 0)) %>%
  dplyr::select(snpID, `Pop ID`, A, B, AA, AB, BB)
```

The complete set of minor allele counts derived from the original genepop file was also loaded


```r
minorAlleleCounts <- file.path("..", "data", "minorAlleleCounts.tsv.gz") %>%
  gzfile() %>%
  read_tsv() %>%
  as.data.frame() %>%
  column_to_rownames("sampleID") %>%
  as.matrix()
```

Populations are defined by the sample names in this object were also defined.


```r
allelePops <- gsub("(gc|ora|TF|Y|pt|tf).+", "\\1", rownames(minorAlleleCounts))
allelePops[!allelePops %in% c("gc", "ora")] <- "tf"
allelePops <- as.factor(allelePops)
names(allelePops) <- rownames(minorAlleleCounts)
```

In addition, the metadata from 2012 sample collection was loaded.


```r
sampleMetadata <- file.path("..", "data", "Samples_Nov2012_GumCreek.xlsx") %>%
  read_excel(skip = 1) %>%
  mutate(ID = gsub("[Oo][Rr][Aa] ", "ora", ID))
```

GPS Locations were added to the metadata.


```r
gpsPoints <-file.path("..", "data", "GPS_Locations.xlsx") %>%
  read_excel() %>%
  dplyr::select(ID = Sample, ends_with("tude")) %>%
  mutate(ID = gsub("[Oo][Rr][Aa] ([0-9ABC]*)", "ora\\1", ID))
sampleMetadata %<>% left_join(gpsPoints, by = "ID")
```



# Principal Component Analysis

## Missing Values

In order to perform PCA, missing values must be either imputed or the entire sample/SNPs must be removed.


```r
missingBySample <- minorAlleleCounts %>%
  apply(MARGIN = 1, function(x){mean(is.na(x))}) 
```

8 samples were missing information for more than 1/3 of SNPs and for the purposes of PCA, these samples were removed from the dataset.



```r
missingByAllele <- minorAlleleCounts %>%
  extract(missingBySample <= 1/3,) %>%
  apply(MARGIN = 2, function(x){mean(is.na(x))}) 
```

Of these remaining samples 9261 alleles were identified in >95% of individuals, and this subset of the data was then chosen as the basis for PCA.


```r
dataForPCA <- minorAlleleCounts[missingBySample <= 1/3, missingByAllele < 0.05]
```

Missing values were then randomly sampled using the distribution of alleles in the remaining samples.
Population structure was not taken into account at this point so as to not artificially inflate correlations within populations.


```r
dataForPCA %<>%
  apply(MARGIN = 2, function(x){
    missing <- is.na(x)
    if (sum(missing) > 0) {
      x[missing] <- sample(x[!missing], sum(missing),replace = TRUE)
    }
    x
  })
```

After this process it was noted that there was no allelic differences amongst the remaining samples for 2 SNPs, implying that the observed variants at these loci were contained within the samples with large numbers of missing data.
These loci were also removed before performing PCA.


```r
dataForPCA %<>% extract(, apply(., MARGIN = 2, function(x){length(unique(x))}) > 1)
```


## PCA Including The Outgroup Samples


```r
set.seed(53)
outgroupPCA <- dataForPCA %>%
  prcomp(center = TRUE)
```

The two Gum Creek populations (1996 & 2012) clearly showed differences to the outgroup, however a strong "tail" was noted for some of the 2012 population along PC2.

The samples forming this tail (i.e. PC2 > 5 & PC3 < 0 ) were identified and the collection region for these samples was investigated, using the GPS collection point.
The vast majority were found to come from the central collection region, and this sub-population was added to the PCA plot.



```r
set.seed(143)
pcaForPlot <- outgroupPCA$x %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  as_data_frame() %>%
  dplyr::select(sampleID, PC1, PC2, PC3) %>%
  mutate(Population = allelePops[sampleID]) %>%
  rowwise() %>%
  mutate(Population = if_else(Population == "gc", "Gum Creek (1996)",
                              if_else(Population == "ora", "Oraparinna (2012)", "Turretfield (2010)"))) %>%
  left_join(sampleMetadata, by = c("sampleID" = "ID")) %>%
  ungroup() %>%
  mutate(Cluster = kmeans(cbind(PC1, PC2, PC3), 3)$cluster)
getRegions <- pcaForPlot %>% 
  filter(grepl("2012", Population)) %>% 
  droplevels() %>% 
  group_by(Cluster) %>% 
  summarise(maxY = max(Latitude)) %>%
  mutate(Region = c("Central", "Outer")[(maxY == max(maxY)) + 1])
pcaForPlot %<>% 
  left_join(getRegions) %>%
  mutate(Population = ifelse(grepl("Oraparinna", Population), 
                             paste("Oraparinna", Region, "(2012)"), Population))
```

![PCA for all samples including the outgroup and indicating the sample collection region for the 2012 samples.](03_snpAnalysis_files/figure-html/finalPCA-1.png)





```r
ausPolygon <- file.path("..", "60804_shp", "framework") %>%
  readOGR(layer = "aus25fgd_r") %>%
  subset(AREA > 6)
ausPts <- SpatialPoints(cbind(x = loc[1], y = loc[2]))
proj4string(ausPts) <- proj4string(ausPolygon)
fullPlot <- ggplot()+
  geom_polygon(data=ausPolygon, 
               aes(long, lat, group = group), fill = "white", colour = "black") + 
  geom_point(data = as.data.frame(ausPts), aes(x, y), size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
```


![Collection points for all 2012 samples with colours showing sub-populations defined by PCA analysis.](03_snpAnalysis_files/figure-html/plotWithInset-1.png)


![Zoomed-in view of the central region for 2012 samples with colours showing sub-populations defined by PCA analysis. The region considered to be the Central Region is shaded in red. Due to overlapping GPS points a small amount of jitter has been added to the x-axis.](03_snpAnalysis_files/figure-html/plotZoom-1.png)

# Phylogenetic Analysis

## Data Export For Generation of Phylogenies

The SNP information for the original 20336 alleles was then exported for phylogenetic analysis.
The same subset of 135 samples was used for this analysis as was used for the PCA analysis above.


```r
alleles <- allData %>%
  dplyr::select(snpID, contains("Nuc")) %>%
  distinct(snpID, .keep_all = TRUE) %>%
  rowwise() %>%
  mutate(AA = paste0(`P Nuc`, `P Nuc`),
         AB = paste0(`P Nuc`, `Q Nuc`),
         BB = paste0(`Q Nuc`, `Q Nuc`)) %>%
  ungroup %>%
  dplyr::select(snpID, AA, AB, BB) %>%
  as.data.frame() %>%
  column_to_rownames("snpID") %>%
  as.matrix()
```


```r
gzOut <- file.path("..", "data", "snpForPhylo.tsv.gz") %>% gzfile("w")
colnames(minorAlleleCounts) %>% 
  vapply(function(x){
    alleles[x,  minorAlleleCounts[,x] + 1]
    }, character(nrow(minorAlleleCounts))) %>%
  as_data_frame() %>%
  mutate(sampleID = rownames(minorAlleleCounts)) %>%
  filter(sampleID %in% pcaForPlot$sampleID) %>%
  mutate(Population = gsub("gc.+", "Gum Creek (1996)", sampleID),
         Population = gsub("(TF|Y|pt|tf).+", "Turretfield (2010)", Population),
         Population = ifelse(sampleID %in% filter(pcaForPlot, grepl("Central", Population))$sampleID,
                             "Oraparinna Central (2012)", Population),
         Population = gsub("ora.+", "Oraparinna Outer (2012)", Population),
         colour = ifelse(grepl("Turretfield", Population), rgb(0.5, 0.5, 1, 0.7), ""),
         colour = ifelse(grepl("Gum", Population), rgb(0, 0, 0, 0.7), colour),
         colour = ifelse(grepl("Central", Population), rgb(0.8, 0.2, 0.2, 0.9), colour),
         colour = ifelse(grepl("Outer", Population), rgb(0.2, 0.8, 0.2, 0.9), colour)) %>%
  dplyr::select(sampleID, Population, colour, everything()) %>%
  write_tsv(gzOut)
close(gzOut)
```


# Analysis

## Removal of SNPs Associated with Collection Region

The structure observed within the 2012 population in the PCA could possibly be explained by recent migration into this region.
As the samples collected in the outer regions appeared very similar to the 1996 population in the above plots, this would possibly indicate this a very recent event as the genetic influence of this has not spread through the wider area.
Although this may be due to other factors such as sampling bias, this structure was addressed by identifying SNPs which showed an association with the sub-populations identified by PCA analysis.
In this way, any candidate SNPs obtained below will be less impacted by this structure, and will be more reflective of the intended variable under study, as opposed to any internal structure of the 2012 population.


```r
oraRegions <- pcaForPlot %>% 
  filter(grepl("2012", Population)) %>% 
  rowwise() %>%
  mutate(yCentral = cut(Latitude, breaks = central["y",], include.lowest = TRUE), 
         xCentral = cut(Longitude, breaks = central["x",], include.lowest = TRUE),
         Central = (is.na(yCentral) + is.na(xCentral)) == 0) %>%
  dplyr::select(sampleID, Central) %>%
  as.data.frame() %>%
  column_to_rownames("sampleID") 
```


### Testing for Structure in 2012

This model tests:  
H<sub>0</sub>: No association between genotypes and collection region  
H<sub>A</sub>: An association exists between genotypes and collection region


```r
ora <- grep("ora", rownames(minorAlleleCounts), value = TRUE)
regionResults <- minorAlleleCounts %>%
  magrittr::extract(ora,) %>%
  apply(MARGIN = 2, function(x){
    mat <- table(oraRegions[ora, "Central"], x)
    if (ncol(mat) > 1){
      fisher.test(mat)$p.value
    }
    else{
      NA
    }
  })
```

A total of 1943 SNPs were detected as showing a significant association between genotype and the collection region.
Under H<sub>0</sub>, the number expected using &#945; = 0.05 would be 1016, and as this number was approximately double that expected, this was taken as evidence of this being a genuine point of concerning this data.

Type II errors were of principle concern in this instance, and as such every SNP with p < 0.05 in the above test was excluded from downstream analysis.


```r
regionSNPs <- names(which(regionResults < 0.05))
regionSNPs %>% writeLines(file.path("..", "results", "regionSNPs.txt"))
```

Under this additional filtering step, **the original set of 20336 SNPs will be reduced to 18393** for testing by genotype and allele frequency.


### Verification Of Removal

In order to verify that the removal of the above SNPs removed the undesired population structure from the 2012 population, the above PCA was repeated, excluding the SNPs marked for removal.
The previous structure noted in the data was no longer evident, and as such, these SNPs were marked for removal during analysis by genotype and allele frequency. 


```r
reducedPCA <- dataForPCA %>%
  magrittr::extract(, !colnames(.) %in% regionSNPs) %>%
  prcomp()
```

![Principal Component Analysis after removal of SNPs showing evidence of structure within the 2012 population.](03_snpAnalysis_files/figure-html/plotReducedPCA-1.png)


## Genotype Frequency Model

This model tests:  
H<sub>0</sub>: No association between genotypes and populations  
H<sub>A</sub>: An association exists between genotypes and populations


```r
genotypeResults <- allCounts %>%
  filter(!snpID %in% regionSNPs) %>%
  split(f = .$snpID) %>%
  mclapply(function(x){
    mat <- as.matrix(x[,c("AA", "AB",  "BB")])
    ft <- fisher.test(mat)
    data_frame(snpID = unique(x$snpID),
               p = ft$p.value)
  }, mc.cores = nCores) %>%
  bind_rows() %>%
  mutate(adjP = p.adjust(p, method = "bonferroni"),
         FDR = p.adjust(p, method = "fdr")) %>%
  arrange(p)
```


Under the full genotype model:

- 0 genotypes were detected as being significantly associated with the two populations when controlling the FWER at &#945; = 0.05
- 0 genotypes were detected as being significantly associated with the two populations when controlling the FDR at &#945; = 0.05
- If controlling the FDR at 10% however, a total of 29 genotypes were considered as potentially associated with the population structure 
- For the most highly ranked SNP (147965_18), the minor allele has been completely lost in the 2012 population


| snpID     |      Chr |          BP |         p |     FDR |
|:----------|---------:|------------:|----------:|--------:|
| 147965_18 |        7 | 131,862,327 | 6.312e-06 | 0.08885 |
| 158509_87 |        5 |  12,946,260 | 1.365e-05 | 0.08885 |
| 104906_36 |       13 | 125,904,892 | 1.811e-05 | 0.08885 |
| 98522_63  |       14 |  34,667,589 | 2.712e-05 | 0.08885 |
| 101831_18 |       14 |  92,793,012 | 2.741e-05 | 0.08885 |
| 72765_47  |       18 |  25,311,139 | 4.157e-05 | 0.08885 |
| 72766_17  |       18 |  25,311,176 | 4.157e-05 | 0.08885 |
| 72765_61  |       18 |  25,311,153 | 5.227e-05 | 0.08885 |
| 37345_16  | GL018754 |   1,365,061 | 5.665e-05 | 0.08885 |
| 157157_45 |        6 |  27,103,576 |  5.68e-05 | 0.08885 |
| 21896_11  | GL018881 |      24,230 | 6.265e-05 | 0.08885 |
| 157156_58 |        6 |  27,103,573 | 6.578e-05 | 0.08885 |
| 80772_39  |       17 |  69,838,525 | 8.448e-05 | 0.08885 |
| 185586_34 |        2 |  47,053,622 | 8.657e-05 | 0.08885 |
| 151791_54 |        7 |  38,250,131 | 9.501e-05 | 0.08885 |
| 151249_77 |        7 |   2,702,916 | 9.629e-05 | 0.08885 |
| 151251_13 |        7 |   2,702,919 | 9.629e-05 | 0.08885 |
| 233206_29 | GL018881 |      88,327 | 9.849e-05 | 0.08885 |
| 37350_8   | GL018754 |   1,395,711 |     1e-04 | 0.08885 |

Table: SNPs with raw p-values < 1e-04 when analysing by genotype. All SNPs were considered significant using an FDR < 0.1

![Manhattan plot showing results for all SNPs on chromosomes 1:21 for analysis by genotype. The horizontal line indicates the cutoff for an FDR of 10%](03_snpAnalysis_files/figure-html/manhattanGenotype-1.png)


## Allele Frequency Model

This model tests:  
H<sub>0</sub>: No association between allele frequencies and populations  
H<sub>A</sub>: An association exists between allele frequencies  and populations


```r
alleleResults <- allCounts %>%
  filter(!snpID %in% regionSNPs) %>%
  split(f = .$snpID) %>%
  mclapply(function(x){
    mat <- as.matrix(x[,c("A", "B")])
    ft <- fisher.test(mat)
    data_frame(snpID = unique(x$snpID),
               p = ft$p.value)
  }, mc.cores = nCores) %>%
  bind_rows() %>%
  mutate(adjP = p.adjust(p, method = "bonferroni"),
         FDR = p.adjust(p, method = "fdr")) %>%
  arrange(p)
```

Under this model:

- 2 SNP alleles were detected as being significantly associated with the two populations when controlling the FWER at &#945; = 0.05.
However, as these SNPs were within 21nt of each other, this may represent the same haplotype 
- 16 SNP alleles were detected as being significantly associated with the two populations when controlling the FDR at &#945; = 0.05
- extending the FDR to 10% yielded 31 SNP alleles



| snpID     |      Chr |          BP |         p |    adjP |     FDR |
|:----------|---------:|------------:|----------:|--------:|--------:|
| 167108_14 |        4 |  84,940,235 | 9.143e-07 | 0.01682 | 0.01075 |
| 167107_60 |        4 |  84,940,214 | 1.169e-06 |  0.0215 | 0.01075 |
| 101831_18 |       14 |  92,793,012 | 4.112e-06 | 0.07564 | 0.02059 |
| 50206_34  | GL018713 |     365,938 | 4.478e-06 | 0.08237 | 0.02059 |
| 98522_63  |       14 |  34,667,589 | 8.562e-06 |  0.1575 |  0.0235 |
| 21896_11  | GL018881 |      24,230 | 1.159e-05 |  0.2132 |  0.0235 |
| 147965_18 |        7 | 131,862,327 | 1.235e-05 |  0.2272 |  0.0235 |
| 233206_29 | GL018881 |      88,327 | 1.266e-05 |  0.2328 |  0.0235 |
| 151249_77 |        7 |   2,702,916 | 1.431e-05 |  0.2632 |  0.0235 |
| 151251_13 |        7 |   2,702,919 | 1.431e-05 |  0.2632 |  0.0235 |
| 80772_39  |       17 |  69,838,525 | 1.468e-05 |  0.2699 |  0.0235 |
| 53831_42  | GL018704 |   4,853,505 | 1.533e-05 |   0.282 |  0.0235 |
| 21896_47  | GL018881 |      24,194 | 2.106e-05 |  0.3874 |   0.029 |
| 72765_61  |       18 |  25,311,153 | 2.434e-05 |  0.4477 |   0.029 |
| 72765_47  |       18 |  25,311,139 | 2.523e-05 |   0.464 |   0.029 |
| 72766_17  |       18 |  25,311,176 | 2.523e-05 |   0.464 |   0.029 |

Table: SNPs considered as significant when analysing by genotype using an FDR cutoff of 0.05


![Manhattan plot showing results for all SNPs on chromosomes 1:21 when analysing by allele frequencies. The horizontal line indicates the cutoff for an FDR of 10%, with SNPs considered significant under the Bonferroni adjustnt shown in green.](03_snpAnalysis_files/figure-html/manhattanAllele-1.png)

## FLK Analysis

The FLK analysis performed separately was also included in order to compare the differing approaches.
It should be noted that this analytic approach is very closely related to analysis by allele frequency, however, instead of Fisher's Exact Test a Chi-squared model is utilised incorporating genetic distances to ccount for genetics drift.


```r
rmarkdown::render("S1_FLK.Rmd")
```



```r
flkResults <- file.path("..", "results", "flkResults.tsv") %>%
  read_tsv()
```

## Bayescan Analysis


```r
snpsForBayescan <- file.path("..", "data", "filteredSNPs.genepop.gz") %>%
  gzfile() %>%
  readLines(n = 2) %>%
  extract2(2) %>%
  strsplit(",") %>%
  extract2(1)
```

The analysis tool Bayescan was also run on the set of 20,336 SNPs exported as a `genepop` file after the previous filtering steps.
The program was run setting the FDR to 0.1.
Any SNPs detected above as tracking with the internal structure of the 2012 population were discarded from the results.


```r
bayes10 <- file.path("..", "results", "BayOut_FDR10.csv") %>%
  read_csv() %>%
  mutate(snpID = snpsForBayescan[SNP])%>%
  dplyr::select(snpID, everything(), -SNP) %>%
  filter(!snpID %in% regionSNPs)
```

## Comparison of Results


```r
fdr <- c(genotype = 0.1, allele = 0.05, flk = 0.05)
sigSNPs <- c(filter(genotypeResults,FDR < fdr["genotype"])$snpID,
             filter(alleleResults, FDR < fdr["allele"])$snpID,
             filter(flkResults, FDR < fdr["flk"])$snpID,
             bayes10$snpID) %>%
  unique %>%
  as.data.frame() %>%
  set_names("snpID") %>%
  mutate(genotype = snpID %in% filter(genotypeResults,FDR < fdr["genotype"])$snpID,
         allele = snpID %in% filter(alleleResults,FDR < fdr["allele"])$snpID,
         flk = snpID %in% filter(flkResults, FDR < fdr["flk"])$snpID,
         bayescan = snpID %in% bayes10$snpID)
```

No novel SNPs were detected by Bayescan or by allele frequency down to an FDR of 0.1 or 0.05 respectively.
As such, the Bayescan & allele frequency analyses were disregarded going forward.

![Overlap between the lists of SNPs considered as associated with the different populations under either of the two analytic approaches, using an FDR of 0.1 for 'Analysis by Genotype' and Bayescan, with an FDR of 0.05 for 'Analysis by Allele' and FLK.](03_snpAnalysis_files/figure-html/upSetSNPs-1.png)

### SNPs Associated With Populations Under Both Approaches

The list of SNPs detected as associated with the population structure under both FLK and analysis by genotype is given below.


| snpID     | Chr      |          BP | Change in log(OR) | P_1996 | P_2012 | Genotype_p |     FLK_p |
|:----------|:---------|------------:|------------------:|-------:|-------:|-----------:|----------:|
| 127156_20 | 10       |  14,373,457 |             1.291 |    0.5 | 0.7843 |  0.0001159 | 8.942e-05 |
| 104906_36 | 13       | 125,904,892 |            -1.492 | 0.9302 |   0.75 |  1.811e-05 | 4.988e-05 |
| 98522_63  | 14       |  34,667,589 |            -1.897 | 0.9375 | 0.6923 |  2.712e-05 | 6.373e-08 |
| 80772_39  | 17       |  69,838,525 |             1.584 | 0.5952 | 0.8776 |  8.448e-05 | 5.139e-05 |
| 72765_47  | 18       |  25,311,139 |            -1.417 | 0.8333 | 0.5481 |  4.157e-05 | 2.369e-06 |
| 72765_61  | 18       |  25,311,153 |            -1.413 | 0.8333 |  0.549 |  5.227e-05 | 2.541e-06 |
| 72766_17  | 18       |  25,311,176 |            -1.417 | 0.8333 | 0.5481 |  4.157e-05 | 2.369e-06 |
| 68946_78  | 19       |  32,825,478 |            -1.547 | 0.9062 | 0.6731 |  0.0001032 | 3.433e-06 |
| 204810_79 | GL018802 |     439,482 |            -1.836 | 0.9583 | 0.7857 |  0.0001362 | 6.111e-06 |
| 21896_47  | GL018881 |      24,194 |            0.9543 | 0.5426 | 0.7549 |  0.0001127 | 4.562e-05 |
| 21896_11  | GL018881 |      24,230 |            0.9158 | 0.5521 | 0.7549 |  6.265e-05 | 2.594e-05 |
| 233206_29 | GL018881 |      88,327 |             1.466 | 0.5357 | 0.8333 |  9.849e-05 | 3.363e-05 |

Table: Summary of changes in the major (P) allele between the two timepoints. Changes in the log Odds ratio of observing the major allele are given, along with estimated population frequencies. Results from testing by genotype or FLK are given as raw p-values. All SNPs were considered as differentially associated with the two populations under both analyses to an FDR of 10\% (genotype) or 5\% (FLK)


![Estimated genotype frequencies for SNPs found to be significant under both models. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/genotypesBothModels-1.png)


### SNPs Associated With Populations Under Analysis Using FLK Only



| Chr      | BP          |     snpID | Change in log(OR) | P_1996 | P_2012 |     FLK_p |
|:---------|:------------|----------:|------------------:|-------:|-------:|----------:|
| 1        | 20,066,891  | 195917_31 |           -0.9545 | 0.7442 | 0.5283 | 5.149e-05 |
| 1        | 88,290,774  | 200447_80 |            -1.354 |   0.86 | 0.6132 |  1.39e-05 |
| 10       | 14,373,501  | 127157_45 |             1.266 |    0.5 |   0.78 | 0.0001145 |
| 10       | 15,009,977  | 127237_50 |            -1.278 |  0.869 | 0.6489 | 6.437e-05 |
| 10       | 15,010,006  | 127237_21 |            -1.303 | 0.8333 | 0.5761 | 1.765e-05 |
| 11       | 36,239,538  |  123550_8 |            -1.427 | 0.9091 | 0.7059 | 3.233e-05 |
| 11       | 36,239,606  | 123551_30 |            -1.317 | 0.8864 | 0.6765 | 6.409e-05 |
| 13       | 15,938,245  | 107422_35 |            -2.351 | 0.9881 | 0.8878 | 6.571e-05 |
| 13       | 76,001,299  | 111958_69 |            -1.341 | 0.8958 | 0.6923 | 6.507e-05 |
| 13       | 128,331,083 | 105215_39 |            0.1309 | 0.6304 | 0.6604 | 5.535e-05 |
| 16       | 735,716     |  87819_88 |            -2.016 | 0.9762 | 0.8452 | 2.811e-05 |
| 16       | 69,850,541  |  87407_11 |            0.2036 |   0.62 | 0.6667 | 7.397e-05 |
| 18       | 68,517,919  | 208110_71 |            -1.304 | 0.8478 |  0.602 | 2.381e-05 |
| 18       | 68,517,942  | 208111_21 |            -1.312 | 0.8478 |    0.6 | 2.064e-05 |
| 2        | 93,561,889  | 214772_69 |            -2.607 | 0.9898 | 0.8774 | 9.255e-06 |
| 20       | 28,859,654  |  65474_28 |            -1.437 | 0.9149 | 0.7188 | 3.942e-05 |
| 21       | 14,643,101  |  62998_77 |            -2.474 | 0.9889 | 0.8824 | 2.511e-05 |
| 3        | 53,556,574  | 175501_56 |            -1.324 | 0.8796 | 0.6604 | 4.414e-05 |
| 4        | 84,940,214  | 167107_60 |             1.987 | 0.6354 | 0.9271 | 1.564e-05 |
| 4        | 84,940,235  | 167108_14 |             1.997 | 0.6383 | 0.9286 |  1.64e-05 |
| 7        | 8,307,536   | 154818_20 |            -1.427 | 0.9135 |  0.717 | 4.198e-05 |
| 9        | 89,862,459  | 211772_50 |            -1.807 | 0.9535 | 0.7708 | 4.259e-06 |
| 9        | 89,862,481  | 211772_28 |            -1.627 | 0.9432 | 0.7653 | 2.053e-05 |
| 9        | 89,862,497  | 211772_12 |            -1.596 | 0.9432 | 0.7708 | 3.351e-05 |
| 9        | 90,037,571  | 137949_40 |             -1.78 | 0.9667 | 0.8302 | 7.017e-05 |
| 9        | 90,037,630  | 137951_16 |            -1.803 | 0.9674 | 0.8302 | 5.925e-05 |
| GL018704 | 4,853,505   |  53831_42 |           -0.3921 | 0.7045 |  0.617 | 4.169e-06 |
| GL018713 | 365,938     |  50206_34 |            -1.066 | 0.7857 | 0.5581 | 1.472e-07 |
| GL018717 | 314,346     |  48659_40 |            -1.471 | 0.8571 | 0.5795 | 1.555e-06 |
| GL018723 | 261,534     | 235050_26 |            -2.152 | 0.9773 | 0.8333 | 5.513e-06 |
| GL018723 | 266,003     | 206301_14 |            -2.047 | 0.9792 | 0.8585 | 4.882e-05 |
| GL018739 | 75,851      |  41475_63 |           -0.7295 | 0.7273 | 0.5625 | 2.309e-05 |
| GL018761 | 75,552      |  36037_47 |            -2.377 | 0.9894 | 0.8962 | 0.0001096 |
| GL018878 | 365,646     | 218710_78 |           -0.1719 | 0.6667 | 0.6275 | 3.554e-05 |
| GL019096 | 125,145     |  11492_22 |            -1.717 | 0.9592 | 0.8085 | 4.587e-05 |
| GL019406 | 5,672       |    5908_8 |            -1.657 |   0.96 | 0.8208 | 0.0001212 |

Table: Summary of changes in the major (P) allele between the two timepoints. Changes in the log Odds ratio of observing the major allele are given, along with estimated population frequencies. Results from testing by allele count are given as raw p-values. All SNPs were considered as differentially associated with the two populations under FLK analysis to an FDR of 5\%.

![Estimated allele frequencies for SNPs found to be significant under analysis by allele only. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/flkModelOnly-1.png)



### SNPs Associated With Populations Under Analysis By Genotype Only



| Chr      | BP          |     snpID | Change in log(OR) | P_1996 | P_2012 | Genotype_p |
|:---------|:------------|----------:|------------------:|-------:|-------:|-----------:|
| 14       | 92,793,012  | 101831_18 |            -2.253 | 0.3778 |   0.06 |  2.741e-05 |
| 18       | 14,288,883  |  72030_54 |            -1.408 | 0.7347 | 0.4038 |  0.0001334 |
| 2        | 37,519,845  |  184763_8 |                -2 | 0.4348 | 0.0943 |  0.0001431 |
| 2        | 47,053,622  | 185586_34 |            -2.611 | 0.3571 | 0.0392 |  8.657e-05 |
| 5        | 12,946,260  | 158509_87 |             3.324 | 0.6429 | 0.9804 |  1.365e-05 |
| 6        | 27,103,573  | 157156_58 |            -1.588 | 0.6981 | 0.3208 |  6.578e-05 |
| 6        | 27,103,576  | 157157_45 |            -1.653 | 0.7115 | 0.3208 |   5.68e-05 |
| 7        | 2,702,916   | 151249_77 |            -2.091 | 0.3269 | 0.0566 |  9.629e-05 |
| 7        | 2,702,919   | 151251_13 |            -2.091 | 0.3269 | 0.0566 |  9.629e-05 |
| 7        | 38,250,131  | 151791_54 |            -2.582 | 0.3462 | 0.0385 |  9.501e-05 |
| 7        | 131,862,327 | 147965_18 |              -Inf | 0.3043 |      0 |  6.312e-06 |
| 9        | 48,543,229  | 134595_50 |             1.541 | 0.3061 | 0.6731 |   0.000107 |
| GL018754 | 1,365,061   |  37345_16 |            -1.453 | 0.7111 | 0.3654 |  5.665e-05 |
| GL018754 | 1,395,679   |  37349_74 |             -1.34 | 0.6939 | 0.3725 |  0.0001063 |
| GL018754 | 1,395,711   |   37350_8 |            -1.369 |    0.7 | 0.3725 |      1e-04 |
| GL018758 | 1,220,835   |  36566_45 |            -1.055 | 0.5962 | 0.3396 |  0.0001331 |
| GL018758 | 1,220,869   |  36567_12 |            -1.055 | 0.5962 | 0.3396 |  0.0001331 |

Table: Summary of changes in heterozygosity between the two timepoints. Changes in the log Odds ratio of observing heterozygotes are given, along with estimated population-level heterozygote frequencies. Results from testing by genotype count are given as raw p-values. All SNPs were considered as differentially associated with the two populations under the genotype count analysis to an FDR of 10\%, However, no significant changes in allele frequencies were detected using an FDR of 5\% for FLK analysis.

![Estimated genotype frequencies for SNPs found to be significant under analysis by genotype only. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/genotypeModelOnly-1.png)


# Output

Both sets of results were output as `genotypeResults.tsv` and `alleleResults.tsv`


```r
alleleResults %>% 
  left_join(allData) %>% 
  distinct(snpID, .keep_all = TRUE) %>% 
  dplyr::select(snpID, Chr, BP, p, adjP, FDR) %>%
  write_tsv(file.path("..", "results", "alleleResults.tsv"))
genotypeResults %>% 
  left_join(allData) %>% 
  distinct(snpID, .keep_all = TRUE) %>% 
  dplyr::select(snpID, Chr, BP, p, adjP, FDR) %>%
  write_tsv( file.path("..", "results", "genotypeResults.tsv"))
```


**R version 3.4.3 (2017-11-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_AU.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_AU.UTF-8_, _LC_COLLATE=en_AU.UTF-8_, _LC_MONETARY=en_AU.UTF-8_, _LC_MESSAGES=en_AU.UTF-8_, _LC_PAPER=en_AU.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_AU.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_stats4_, _grid_, _parallel_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_ape(v.5.0)_, _rtracklayer(v.1.38.2)_, _GenomicRanges(v.1.30.1)_, _GenomeInfoDb(v.1.14.0)_, _IRanges(v.2.12.0)_, _S4Vectors(v.0.16.0)_, _BiocGenerics(v.0.24.0)_, _bindrcpp(v.0.2)_, _rgdal(v.1.2-16)_, _ggmap(v.2.6.1)_, _sp(v.1.2-7)_, _qqman(v.0.1.4)_, _UpSetR(v.1.3.3)_, _magrittr(v.1.5)_, _readxl(v.1.0.0)_, _reshape2(v.1.4.3)_, _scales(v.0.5.0)_, _pander(v.0.6.1)_, _forcats(v.0.2.0)_, _stringr(v.1.2.0)_, _dplyr(v.0.7.4)_, _purrr(v.0.2.4)_, _readr(v.1.1.1)_, _tidyr(v.0.7.2)_, _tibble(v.1.4.1)_, _ggplot2(v.2.2.1)_ and _tidyverse(v.1.2.1)_

**loaded via a namespace (and not attached):** 
_nlme(v.3.1-131)_, _bitops(v.1.0-6)_, _matrixStats(v.0.52.2)_, _lubridate(v.1.7.1)_, _httr(v.1.3.1)_, _rprojroot(v.1.3-2)_, _tools(v.3.4.3)_, _backports(v.1.1.2)_, _R6(v.2.2.2)_, _lazyeval(v.0.2.1)_, _colorspace(v.1.3-2)_, _gridExtra(v.2.3)_, _mnormt(v.1.5-5)_, _compiler(v.3.4.3)_, _cli(v.1.0.0)_, _rvest(v.0.3.2)_, _Biobase(v.2.38.0)_, _xml2(v.1.1.1)_, _DelayedArray(v.0.4.1)_, _labeling(v.0.3)_, _psych(v.1.7.8)_, _digest(v.0.6.14)_, _Rsamtools(v.1.30.0)_, _foreign(v.0.8-69)_, _rmarkdown(v.1.8)_, _XVector(v.0.18.0)_, _jpeg(v.0.1-8)_, _pkgconfig(v.2.0.1)_, _htmltools(v.0.3.6)_, _highr(v.0.6)_, _maps(v.3.2.0)_, _rlang(v.0.1.6)_, _rstudioapi(v.0.7)_, _bindr(v.0.1)_, _jsonlite(v.1.5)_, _BiocParallel(v.1.12.0)_, _RCurl(v.1.95-4.10)_, _GenomeInfoDbData(v.1.0.0)_, _geosphere(v.1.5-7)_, _Matrix(v.1.2-12)_, _Rcpp(v.0.12.15)_, _munsell(v.0.4.3)_, _proto(v.1.0.0)_, _stringi(v.1.1.6)_, _yaml(v.2.1.16)_, _SummarizedExperiment(v.1.8.1)_, _zlibbioc(v.1.24.0)_, _plyr(v.1.8.4)_, _crayon(v.1.3.4)_, _lattice(v.0.20-35)_, _Biostrings(v.2.46.0)_, _haven(v.1.1.1)_, _mapproj(v.1.2-5)_, _hms(v.0.4.0)_, _knitr(v.1.18)_, _pillar(v.1.1.0)_, _rjson(v.0.2.15)_, _XML(v.3.98-1.9)_, _glue(v.1.2.0)_, _evaluate(v.0.10.1)_, _calibrate(v.1.7.2)_, _modelr(v.0.1.1)_, _png(v.0.1-7)_, _RgoogleMaps(v.1.4.1)_, _cellranger(v.1.1.0)_, _gtable(v.0.2.0)_, _assertthat(v.0.2.0)_, _broom(v.0.4.3)_ and _GenomicAlignments(v.1.14.1)_

