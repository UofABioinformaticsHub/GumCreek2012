# SNP Analysis
Steve Pederson  
`r format(Sys.Date(), "%b %d, %Y")`  




```r
library(tidyverse)
library(parallel)
library(pander)
library(scales)
library(reshape2)
library(readxl)
library(magrittr)
library(VennDiagram)
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

In addition, the metadata from sample collection was loaded.


```r
sampleMetadata <- file.path("..", "data", "Samples_Nov2012_GumCreek.xlsx") %>%
  read_excel(skip = 1) %>%
  mutate(ID = gsub("[Oo][Rr][Aa] ", "ora", ID))
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
outgroupPCA <- dataForPCA %>%
  prcomp(center = TRUE)
```

The two Gum Creek populations (1996 & 2012) clearly showed differences to the outgroup, however a strong "tail" was noted for some of the 2012 population along PC2.

The samples forming this tail (i.e. PC2 > 5) were identified and the collection region for these samples was investigated, using the `GPS waypoint` value and the map of waypoints.
The vast majority were found to come from the central collection region, and this collection region was added to the PCA plot.


```r
central <- seq(44, 61) %>%
  setdiff(c(45, 52, 53)) %>%
  c(74:76)
```


```r
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
  mutate(Central = `GPS waypoint` %in% central,
         Population = if_else(Population == "Oraparinna (2012)",
                              if_else(Central, "Oraparinna Central (2012)", "Oraparinna Outer (2012)" ), Population)) 
```

![PCA for all samples including the outgroup and indicating the sample collection region for the 2012 samples.](03_snpAnalysis_files/figure-html/finalPCA-1.png)


# Analysis

## Removal of SNPs Associated with Collection Region

The structure observed in the 2012 population in the PCA could possibly be explained by recent migration into this region.
As the samples collected in the outer regions appeared very similar to the 1996 population in the above plots, this would possibly indicate this a very recent event as the genetic influence of this has not spread through the wider area.
Although this may be due to other factors such as sampling bias, this structure was addressed by identifying SNPs which showed an association with the Central/Outer collection regions in this population.
In this way, any candidate SNPs obtained below will be less impacted by this structure, and may be more reflective of the intended variable under study, namely selective pressure, not migration.


```r
oraRegions <- sampleMetadata %>%
  mutate(Central = `GPS waypoint` %in% central) %>%
  dplyr::select(ID, Central) %>%
  as.data.frame() %>%
  column_to_rownames("ID") 
```


### Testing for Structure in 2012

This model tests:  
H<sub>0</sub>: No association between genotypes and collection region  
H<sub>A</sub>: An association exists between genotypes and collection region


```r
ora <- grep("ora", rownames(minorAlleleCounts), value = TRUE)
regionResults <- minorAlleleCounts %>%
  extract(ora,) %>%
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
Under H<sub>0</sub>, the number expected using &#945; = 0.05 would be 1016, and as this number was nearly double that expected, this was taken as evidence of this being a genuine point of concerning this data.

Type II errors were of principle concern in this instance, and as such every SNP with p < 0.05 in the above test was excluded from downstream analysis.


```r
regionSNPs <- names(which(regionResults < 0.05))
```

Under this additional filtering step, **the original set of 20336 SNPs will be reduced to 18393** for testing by genotype and allele frequency.


### Verification Of Removal

In order to verify that the removal of the above SNPs removed the undesired population structure from the 2012 population, the above PCA was repeated, excluding the SNPs marked for removal.
The previous structure noted in the data was no longer evident, and as such, these SNPs were marked for removal during analysis by genotype and allele frequency. 


```r
reducedPCA <- dataForPCA %>%
  extract(, !colnames(.) %in% regionSNPs) %>%
  prcomp()
```

![Principal Component Analysis after removal of SNPs showing evidence of structure within the 2012 population.](03_snpAnalysis_files/figure-html/reducedPCA-1.png)


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


|      Chr |          BP | snpID     |         p |     FDR |
|---------:|------------:|:----------|----------:|--------:|
|        7 | 131,862,327 | 147965_18 | 6.312e-06 | 0.08885 |
|        5 |  12,946,260 | 158509_87 | 1.365e-05 | 0.08885 |
|       13 | 125,904,892 | 104906_36 | 1.811e-05 | 0.08885 |
|       14 |  34,667,589 | 98522_63  | 2.712e-05 | 0.08885 |
|       14 |  92,793,012 | 101831_18 | 2.741e-05 | 0.08885 |
|       18 |  25,311,139 | 72765_47  | 4.157e-05 | 0.08885 |
|       18 |  25,311,176 | 72766_17  | 4.157e-05 | 0.08885 |
|       18 |  25,311,153 | 72765_61  | 5.227e-05 | 0.08885 |
| GL018754 |   1,365,061 | 37345_16  | 5.665e-05 | 0.08885 |
|        6 |  27,103,576 | 157157_45 |  5.68e-05 | 0.08885 |
| GL018881 |      24,230 | 21896_11  | 6.265e-05 | 0.08885 |
|        6 |  27,103,573 | 157156_58 | 6.578e-05 | 0.08885 |
|       17 |  69,838,525 | 80772_39  | 8.448e-05 | 0.08885 |
|        2 |  47,053,622 | 185586_34 | 8.657e-05 | 0.08885 |
|        7 |  38,250,131 | 151791_54 | 9.501e-05 | 0.08885 |
|        7 |   2,702,916 | 151249_77 | 9.629e-05 | 0.08885 |
|        7 |   2,702,919 | 151251_13 | 9.629e-05 | 0.08885 |
| GL018881 |      88,327 | 233206_29 | 9.849e-05 | 0.08885 |
| GL018754 |   1,395,711 | 37350_8   |     1e-04 | 0.08885 |

Table: SNPs with raw p-values < 1e-04 when analysing by genotype. All SNPs were considered significant using an FDR < 0.1

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



|      Chr |          BP | snpID     |         p |    adjP |     FDR |
|---------:|------------:|:----------|----------:|--------:|--------:|
|        4 |  84,940,235 | 167108_14 | 9.143e-07 | 0.01682 | 0.01075 |
|        4 |  84,940,214 | 167107_60 | 1.169e-06 |  0.0215 | 0.01075 |
|       14 |  92,793,012 | 101831_18 | 4.112e-06 | 0.07564 | 0.02059 |
| GL018713 |     365,938 | 50206_34  | 4.478e-06 | 0.08237 | 0.02059 |
|       14 |  34,667,589 | 98522_63  | 8.562e-06 |  0.1575 |  0.0235 |
| GL018881 |      24,230 | 21896_11  | 1.159e-05 |  0.2132 |  0.0235 |
|        7 | 131,862,327 | 147965_18 | 1.235e-05 |  0.2272 |  0.0235 |
| GL018881 |      88,327 | 233206_29 | 1.266e-05 |  0.2328 |  0.0235 |
|        7 |   2,702,916 | 151249_77 | 1.431e-05 |  0.2632 |  0.0235 |
|        7 |   2,702,919 | 151251_13 | 1.431e-05 |  0.2632 |  0.0235 |
|       17 |  69,838,525 | 80772_39  | 1.468e-05 |  0.2699 |  0.0235 |
| GL018704 |   4,853,505 | 53831_42  | 1.533e-05 |   0.282 |  0.0235 |
| GL018881 |      24,194 | 21896_47  | 2.106e-05 |  0.3874 |   0.029 |
|       18 |  25,311,153 | 72765_61  | 2.434e-05 |  0.4477 |   0.029 |
|       18 |  25,311,139 | 72765_47  | 2.523e-05 |   0.464 |   0.029 |
|       18 |  25,311,176 | 72766_17  | 2.523e-05 |   0.464 |   0.029 |

Table: SNPs considered as significant when analysing by genotype using an FDR cutoff of 0.05


## Comparison of Results


```r
fdr <- c(genotype = 0.1, allele = 0.1)
sigSNPs <- c(filter(genotypeResults,FDR < fdr["genotype"])$snpID,
             filter(alleleResults, FDR < fdr["allele"])$snpID) %>%
  unique %>%
  as.data.frame() %>%
  set_names("snpID") %>%
  mutate(genotype = snpID %in% filter(genotypeResults,FDR < fdr["genotype"])$snpID,
         allele = snpID %in% filter(alleleResults,FDR < fdr["allele"])$snpID)
```

![Overlap between the lists of SNPs considered as associated with the different populations under either of the two analytic approaches, using an FDR of 0.1.](03_snpAnalysis_files/figure-html/vennSNPs-1.png)

### SNPs Associated With Populations Under Both Approaches

The list of SNPs detected as associated with the population structure under both approaches is given below.


| snpID     | Chr      |          BP | Change in log(OR) | P_1996 | P_2012 |  Allele_p | Genotype_p |
|:----------|:---------|------------:|------------------:|-------:|-------:|----------:|-----------:|
| 127156_20 | 10       |  14,373,457 |            -1.291 |    0.5 | 0.2157 | 5.695e-05 |  0.0001159 |
| 98522_63  | 14       |  34,667,589 |            -1.897 | 0.9375 | 0.6923 | 8.562e-06 |  2.712e-05 |
| 101831_18 | 14       |  92,793,012 |             2.407 | 0.7444 |   0.97 | 4.112e-06 |  2.741e-05 |
| 80772_39  | 17       |  69,838,525 |             1.584 | 0.5952 | 0.8776 | 1.468e-05 |  8.448e-05 |
| 72765_47  | 18       |  25,311,139 |            -1.417 | 0.8333 | 0.5481 | 2.523e-05 |  4.157e-05 |
| 72765_61  | 18       |  25,311,153 |            -1.413 | 0.8333 |  0.549 | 2.434e-05 |  5.227e-05 |
| 72766_17  | 18       |  25,311,176 |            -1.417 | 0.8333 | 0.5481 | 2.523e-05 |  4.157e-05 |
| 68946_78  | 19       |  32,825,478 |            -1.547 | 0.9062 | 0.6731 | 5.891e-05 |  0.0001032 |
| 151249_77 | 7        |   2,702,916 |             2.277 | 0.7788 | 0.9717 | 1.431e-05 |  9.629e-05 |
| 151251_13 | 7        |   2,702,919 |             2.277 | 0.7788 | 0.9717 | 1.431e-05 |  9.629e-05 |
| 147965_18 | 7        | 131,862,327 |               Inf | 0.8478 |      1 | 1.235e-05 |  6.312e-06 |
| 36566_45  | GL018758 |   1,220,835 |            -1.162 | 0.5096 | 0.2453 | 0.0001083 |  0.0001331 |
| 36567_12  | GL018758 |   1,220,869 |            -1.162 | 0.5096 | 0.2453 | 0.0001083 |  0.0001331 |
| 21896_47  | GL018881 |      24,194 |            -1.296 | 0.5426 | 0.2451 | 2.106e-05 |  0.0001127 |
| 21896_11  | GL018881 |      24,230 |            -1.334 | 0.5521 | 0.2451 | 1.159e-05 |  6.265e-05 |
| 233206_29 | GL018881 |      88,327 |             1.466 | 0.5357 | 0.8333 | 1.266e-05 |  9.849e-05 |

Table: Summary of changes in the major (P) allele between the two timepoints. Changes in the log Odds ratio of observing the major allele are given, along with estimated population frequencies. Results from testing by genotype or allele are given as raw p-values. All SNPs were considered as differentially associated with the two populations under both analyses to an FDR of 10\%.


![Estimated genotype frequencies for SNPs found to be significant under both models. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/genotypesBothModels-1.png)



### SNPs Associated With Populations Under Analysis By Allele Only



| Chr      | BP          |     snpID | Change in log(OR) | P_1996 | P_2012 |  Allele_p |
|:---------|:------------|----------:|------------------:|-------:|-------:|----------:|
| 1        | 88,290,774  | 200447_80 |            -1.354 |   0.86 | 0.6132 | 6.994e-05 |
| 10       | 14,373,501  | 127157_45 |            -1.266 |    0.5 |   0.22 | 9.337e-05 |
| 13       | 128,331,083 | 105215_39 |            -1.199 | 0.6304 | 0.3396 | 5.926e-05 |
| 16       | 69,850,541  |  87407_11 |            -1.183 |   0.62 | 0.3333 | 7.197e-05 |
| 16       | 69,850,609  |  208806_6 |             -1.14 |   0.61 | 0.3333 | 0.0001291 |
| 3        | 53,556,574  | 175501_56 |            -1.324 | 0.8796 | 0.6604 | 0.0001568 |
| 4        | 84,940,214  | 167107_60 |             1.987 | 0.6354 | 0.9271 | 1.169e-06 |
| 4        | 84,940,235  | 167108_14 |             1.997 | 0.6383 | 0.9286 | 9.143e-07 |
| GL018704 | 4,853,505   |  53831_42 |            -1.346 | 0.7045 |  0.383 | 1.533e-05 |
| GL018713 | 365,938     |  50206_34 |            -1.533 | 0.7857 | 0.4419 | 4.478e-06 |
| GL018717 | 314,346     |  48659_40 |            -1.471 | 0.8571 | 0.5795 | 8.114e-05 |
| GL018739 | 75,851      |  41475_63 |            -1.232 | 0.7273 | 0.4375 | 9.437e-05 |
| GL018878 | 365,646     | 218710_78 |            -1.214 | 0.6667 | 0.3725 |  7.13e-05 |
| GL018881 | 22,826      |  21889_62 |            -1.143 | 0.5755 | 0.3019 | 9.673e-05 |
| GL018985 | 50,536      |   15651_5 |             1.432 | 0.6735 | 0.8962 | 0.0001242 |

Table: Summary of changes in the major (P) allele between the two timepoints. Changes in the log Odds ratio of observing the major allele are given, along with estimated population frequencies. Results from testing by allele count are given as raw p-values. All SNPs were considered as differentially associated with the two populations under the allele count analysis to an FDR of 10\%.

![Estimated allele frequencies for SNPs found to be significant under analysis by allele only. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/alleleModelOnly-1.png)



### SNPs Associated With Populations Under Analysis By Genotype Only



| Chr      | BP          |     snpID | Change in log(OR) | P_1996 | P_2012 | Genotype_p |
|:---------|:------------|----------:|------------------:|-------:|-------:|-----------:|
| 13       | 125,904,892 | 104906_36 |             2.278 |  0.093 |    0.5 |  1.811e-05 |
| 18       | 14,288,883  |  72030_54 |            -1.408 | 0.7347 | 0.4038 |  0.0001334 |
| 2        | 37,519,845  |  184763_8 |                -2 | 0.4348 | 0.0943 |  0.0001431 |
| 2        | 47,053,622  | 185586_34 |            -2.611 | 0.3571 | 0.0392 |  8.657e-05 |
| 5        | 12,946,260  | 158509_87 |             3.324 | 0.6429 | 0.9804 |  1.365e-05 |
| 6        | 27,103,573  | 157156_58 |            -1.588 | 0.6981 | 0.3208 |  6.578e-05 |
| 6        | 27,103,576  | 157157_45 |            -1.653 | 0.7115 | 0.3208 |   5.68e-05 |
| 7        | 38,250,131  | 151791_54 |            -2.582 | 0.3462 | 0.0385 |  9.501e-05 |
| 9        | 48,543,229  | 134595_50 |             1.541 | 0.3061 | 0.6731 |   0.000107 |
| GL018754 | 1,365,061   |  37345_16 |            -1.453 | 0.7111 | 0.3654 |  5.665e-05 |
| GL018754 | 1,395,679   |  37349_74 |             -1.34 | 0.6939 | 0.3725 |  0.0001063 |
| GL018754 | 1,395,711   |   37350_8 |            -1.369 |    0.7 | 0.3725 |      1e-04 |
| GL018802 | 439,482     | 204810_79 |             2.111 | 0.0833 | 0.4286 |  0.0001362 |

Table: Summary of changes in heterozygosity between the two timepoints. Changes in the log Odds ratio of observing heterozygotes are given, along with estimated population-level heterozygote frequencies. Results from testing by genotype count are given as raw p-values. All SNPs were considered as differentially associated with the two populations under the genotype count analysis to an FDR of 10\%, However, no significant changes in allele frequencies were detected using an FDR of 10\% for the analysis of allele counts. The possibility of the allele-level results being Type II errors cannot be discounted.

![Estimated genotype frequencies for SNPs found to be significant under analysis by genotype only. In each case the `A` allele represents the major allele in the 1996 population.](03_snpAnalysis_files/figure-html/genotypeModelOnly-1.png)


# Output

Both sets of results were output as `genotypeResults.tsv` and `alleleResults.tsv`


```r
alleleResults %>% 
  left_join(allData) %>% 
  distinct(snpID, .keep_all = TRUE) %>% 
  dplyr::select(snpID, Chr, BP, p, adjP, FDR) %>%
  write_tsv(file.path("..", "results", "alleleResults.tsv"))
```

```r
genotypeResults %>% 
  left_join(allData) %>% 
  distinct(snpID, .keep_all = TRUE) %>% 
  dplyr::select(snpID, Chr, BP, p, adjP, FDR) %>%
  write_tsv( file.path("..", "results", "genotypeResults.tsv"))
```



```r
pander(sessionInfo()) 
```

**R version 3.4.1 (2017-06-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_AU.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_AU.UTF-8_, _LC_COLLATE=en_AU.UTF-8_, _LC_MONETARY=en_AU.UTF-8_, _LC_MESSAGES=en_AU.UTF-8_, _LC_PAPER=en_AU.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_AU.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_grid_, _parallel_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_bindrcpp(v.0.2)_, _VennDiagram(v.1.6.17)_, _futile.logger(v.1.4.3)_, _magrittr(v.1.5)_, _readxl(v.1.0.0)_, _reshape2(v.1.4.2)_, _scales(v.0.5.0)_, _pander(v.0.6.1)_, _dplyr(v.0.7.2)_, _purrr(v.0.2.3)_, _readr(v.1.1.1)_, _tidyr(v.0.7.1)_, _tibble(v.1.3.4)_, _ggplot2(v.2.2.1)_ and _tidyverse(v.1.1.1)_

**loaded via a namespace (and not attached):** 
_haven(v.1.1.0)_, _lattice(v.0.20-35)_, _colorspace(v.1.3-2)_, _htmltools(v.0.3.6)_, _yaml(v.2.1.14)_, _rlang(v.0.1.2)_, _foreign(v.0.8-69)_, _glue(v.1.1.1)_, _modelr(v.0.1.1)_, _lambda.r(v.1.1.9)_, _bindr(v.0.1)_, _plyr(v.1.8.4)_, _stringr(v.1.2.0)_, _munsell(v.0.4.3)_, _gtable(v.0.2.0)_, _cellranger(v.1.1.0)_, _rvest(v.0.3.2)_, _psych(v.1.7.5)_, _evaluate(v.0.10.1)_, _labeling(v.0.3)_, _knitr(v.1.17)_, _forcats(v.0.2.0)_, _highr(v.0.6)_, _broom(v.0.4.2)_, _Rcpp(v.0.12.12)_, _backports(v.1.1.0)_, _jsonlite(v.1.5)_, _mnormt(v.1.5-5)_, _hms(v.0.3)_, _digest(v.0.6.12)_, _stringi(v.1.1.5)_, _rprojroot(v.1.2)_, _tools(v.3.4.1)_, _lazyeval(v.0.2.0)_, _futile.options(v.1.0.0)_, _pkgconfig(v.2.0.1)_, _xml2(v.1.1.1)_, _lubridate(v.1.6.0)_, _assertthat(v.0.2.0)_, _rmarkdown(v.1.6)_, _httr(v.1.3.1)_, _R6(v.2.2.2)_, _nlme(v.3.1-131)_ and _compiler(v.3.4.1)_

