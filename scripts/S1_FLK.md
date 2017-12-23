---
title: "Supplementary: FLK Analysis"
author: "Steve Pederson"
date: "29 September 2017"
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
library(magrittr)
library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(reshape2)
library(pander)
library(scales)
library(stringr)
```


```r
nCores <- min(12, detectCores() - 1)
```

# Introduction

This supplementary analysis uses the FLK model as outlined in [Bonhomme et al, Detecting Selection in Population Trees: The Lewontin and Krakauer Test Extended](https://dx.doi.org/10.1534%2Fgenetics.110.117275).
This test is robust to genetic drift, and is able to detect signatures of selection.
However as this is a linear evolution dataset within one population as opposed to a series of approximately parallel selection events occurring within a series of populations, the assumptions of the FLK model may be less than satisfied.
Hence this is presented as supplementary information.

A third population was included as an outlier group, representing a population separated by a large geographic distance. 35 samples were taken from [Schwensow et al]( http://onlinelibrary.wiley.com/doi/10.1111/mec.14228/full).

# Setup

## Genome & Genes

Whilst aligned to the RefSeq/NCBI genome build, for compatibility with Ensembl gene models, a mapping object was generated from the [RefSeq Assembly Report](ftp://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Oryctolagus_cuniculus/latest_assembly_versions/GCF_000003625.3_OryCun2.0/GCF_000003625.3_OryCun2.0_assembly_report.txt)


```r
ncbiReport <- file.path("..", "data", "GCF_000003625.3_OryCun2.0_assembly_report.txt") %>%
  read_delim(delim = "\t", col_names = FALSE, na = "na", comment = "#") %>%
  set_colnames(c("Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name")) %>%
  filter(!is.na(`RefSeq-Accn`))
ensUsesGB <- grepl("scaffold",ncbiReport$`Sequence-Role`)
ncbiReport$Ensembl <- ncbiReport$`Sequence-Name`
ncbiReport$Ensembl[ensUsesGB] <- gsub("\\.[0-9]$", "", ncbiReport$`GenBank-Accn`[ensUsesGB])
ncbi2Ens <- with(ncbiReport, data.frame(Chr = Ensembl, 
                                        Length = `Sequence-Length`,
                                        row.names = `RefSeq-Accn`))
```

Gene information was obtained from Ensembl Build 84 and loaded 


```r
ensGff <- file.path("..", "data", "Oryctolagus_cuniculus.OryCun2.0.84.gff3.gz")
ensGenes <- import.gff3(ensGff, feature.type = "gene", sequenceRegionsAsSeqinfo = TRUE) 
```


## Data Loading

The same stacks output as in [2_snpFiltering](2_snpFiltering.md) was loaded, retaining data from the third population (Turretfield).


```r
allData <- file.path("..", "data", "batch_3.sumstats.tsv.gz") %>%
  gzfile() %>%
  read_delim(delim = "\t", col_names = TRUE, skip = 3) %>%
  dplyr::select(-contains("Batch"), -contains("Pi"), -contains("Fis")) %>%
  mutate(RefSeq = gsub("^gi\\|.+\\|ref\\|(.+)\\|$", "\\1", Chr),
         Chr = as.character(ncbi2Ens[RefSeq, "Chr"]),
         snpID = paste(`Locus ID`, Col, sep="_"),
         Private = as.logical(Private)) %>%
  filter(Chr != "X",
         Col > 4) %>%
  dplyr::select(Chr, BP, `Locus ID`, Col, snpID, 
                `Pop ID`, contains("Nuc"),  N, P, 
                contains("Obs"), contains("Exp"), Private)
```


```r
allSnpsGR <- allData %>%
  distinct(snpID, .keep_all = TRUE) %>%
  makeGRangesFromDataFrame(ignore.strand = TRUE,
                           keep.extra.columns = TRUE,
                           seqinfo = seqinfo(ensGenes), 
                           seqnames.field = "Chr", 
                           start.field = "BP", 
                           end.field = "BP") %>%
  sort()
mcols(allSnpsGR) <- DataFrame(snpID = allSnpsGR$snpID)
```


## FLK

The script `FLK.R` was obtained from the [QGSP]( https://qgsp.jouy.inra.fr/archives/FLK/FLK.R) and loaded into `R`


```r
source("FLK.R")
```


## Neutral Loci


```r
hitsIn100kb <- findOverlaps(
  allSnpsGR %>%
    resize(width = 200001,fix = "center") %>%
    trim(),
    ensGenes
    )
```


```r
neutLoci <- allSnpsGR %>%
  extract(setdiff(seq_along(.), queryHits(hitsIn100kb))) %>%
  mcols() %>%
  extract2("snpID")
```


```r
neutData <- allData %>%
  filter(snpID %in% neutLoci) %>%
  distinct(Chr, BP, `Pop ID`, .keep_all = TRUE)
```


As a broad definition of neutral loci, the set of SNPs > 100kb from any known gene were selected.
After removal of duplicate loci, this gave a set of 26,097 SNPs.

Neutral SNPs were again checked to ensure that the `P` allele was the same as used in the 1996 population.


```r
neutData <- neutData %>%
  split(f = .$snpID) %>%
  mclapply(function(x){
    if (length(unique(x$`P Nuc`)) != 3){
     pop2Change <- which(x$`P Nuc` != x$`P Nuc`[1])
     x$`P Nuc`[pop2Change] <- x$`P Nuc`[1]
     x$`Q Nuc`[pop2Change] <- x$`Q Nuc`[1]
     x$P[pop2Change] <- 1 - x$P[pop2Change] 
    }
    x
  }, mc.cores = nCores) %>%
  bind_rows() %>%
  arrange(Chr, BP, `Pop ID`)
```


```r
N <- allData %>% 
  group_by(`Pop ID`) %>% 
  summarise(maxN = max(N)) %>%
  mutate(minN = round(0.95*maxN))
```

The list of Neutral SNPs were then restricted to those identified in >95% of individuals within all populations.


```r
neutMatrix <- neutData %>%
  left_join(N) %>%
  group_by(snpID) %>%
  filter(all(N > minN)) %>%
  mutate(`Pop ID` = c("Gum Creek (1996)", "Oraparinna (2012)", "Turretfield (2010)")[`Pop ID`]) %>%
  dcast(snpID~`Pop ID`, value.var = "P") %>%
  filter(`Gum Creek (1996)` != 1) %>%
  column_to_rownames("snpID") %>%
  as.matrix() %>%
  t()
```

This gave a final list of 1171 neutral SNPs for estimation of Reynolds Distance

## Reynolds Distance 


```r
reynDist <- reynolds(neutMatrix)
```


------------------------------------------------------------------------------------
&nbsp;                     Gum Creek (1996)   Oraparinna (2012)   Turretfield (2010)
------------------------ ------------------ ------------------- --------------------
**Gum Creek (1996)**                      0             0.01062              0.08543

**Oraparinna (2012)**               0.01062                   0              0.09209

**Turretfield (2010)**              0.08543             0.09209                    0
------------------------------------------------------------------------------------

Table: Reynolds Distance as calculated using the neutral loci as defined above.

## Co-ancestry Matrix (F_ij)


```r
F_ij <- Fij(neutMatrix, "Turretfield (2010)", reynDist)
```

![Neighbout-joining Tree for all 3 populations.](S1_FLK_files/figure-html/njTree-1.png)

## FLK Results


```r
testData <-  file.path("..", "data", "filteredSNPs.tsv.gz") %>%
  gzfile() %>%
  read_delim(delim = "\t") 
```

The set of 2.0336\times 10^{4} SNPs previously obtained for testing after all filtering steps were then tested using the FLK model.
P-values obtained under FLK were adjusted using Bonferroni's method to obtain a set of high-confidence SNPs, then using Benjamini-Hochberg's FDR to provide a larger set for pathway testing.


```r
flkResults <- testData %>%
  dplyr::select(snpID, `Pop ID`, P) %>%
  mutate(`Pop ID` = c("Gum Creek (1996)", "Oraparinna (2012)")[`Pop ID`]) %>%
  acast(`Pop ID` ~ snpID ) %>%
  FLK(Fij = F_ij) %>%
  rownames_to_column("snpID") %>%
  as_tibble() %>%
  dplyr::select(snpID, Ht, contains("F.LK")) %>%
  mutate(FDR = p.adjust(F.LK.p.val, method = "fdr"),
         P_bonf = p.adjust(F.LK.p.val, method = "bonferroni")) %>%
  arrange(F.LK.p.val)
```

A total of 6 SNPs retained significance after the Bonferroni adjustment with a total of 62 being considered in the larger set of SNPs with an FDR of 0.05.


```r
testData %>%
  dplyr::select(snpID, Chr, BP, `Pop ID`, `P Nuc`, `Q Nuc`, P) %>%
  mutate(`Pop ID` = c("Gum Creek (1996)", "Oraparinna (2012)")[`Pop ID`],
         SNP = paste(`P Nuc`, `Q Nuc`, sep = "/")) %>%
  dcast(snpID + Chr + BP + SNP  ~ `Pop ID`, value.var = "P") %>%
  right_join(flkResults) %>%
  as_data_frame() %>%
  arrange(F.LK.p.val) %>%
  write_tsv( file.path("..", "results", "flkResults.tsv"))
```


