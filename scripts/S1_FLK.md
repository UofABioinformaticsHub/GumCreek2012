# Supplementary: FLK Analysis
Steve Pederson  
29 September 2017  




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
However as this is a linear evolution dataset within one population as opposed to a series of approximately parallel selection event occurring within a series of populations, the assumptions of the FLK model may be less than satisfied.
Hence this is presented as supplementary information.

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


As a broad definition of neutral loci, the set of SNPs not within 100kb of a known gene were selected.
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
  mutate(`Pop ID` = c("Gum Creek", "Oraparinna", "Turretfield")[`Pop ID`]) %>%
  dcast(snpID~`Pop ID`, value.var = "P") %>%
  filter(`Gum Creek` != 1) %>%
  column_to_rownames("snpID") %>%
  as.matrix() %>%
  t()
```

This gave a final list of 1171 neutral SNPs for estimation of Reynolds Distance

## Reynolds Distance and F_ij


```r
reynDist <- reynolds(neutMatrix)
```


--------------------------------------------------------
     &nbsp;        Gum Creek   Oraparinna   Turretfield 
----------------- ----------- ------------ -------------
  **Gum Creek**        0        0.01062       0.08543   

 **Oraparinna**     0.01062        0          0.09209   

 **Turretfield**    0.08543     0.09209          0      
--------------------------------------------------------

Table: Reynolds Distance as calculated using the neutral loci as defined above.


```r
F_ij <- Fij(neutMatrix, "Turretfield", reynDist)
```

![Neighbout-joining Tree for all 3 populations.](S1_FLK_files/figure-html/njTree-1.png)

## Co-ancestry Matrix

## FLK Results