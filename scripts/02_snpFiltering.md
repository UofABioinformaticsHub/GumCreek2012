---
title: "SNP Filtering"
author: "Steve Pederson"
date: "29 January, 2018"
output: 
  html_document: 
    keep_md: yes
    toc: yes
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
nCores <- min(12, detectCores() - 1)
```


# Introduction

This file takes the raw output from `Stacks` and filters the complete set of SNPs to arrive at the final list of SNPs used for downstream analysis.
The two populations are referred to below as Population `1` (*Gum Creek*, 1996) and Population `2` (*Oraparinna*, 2012).

## Outline of selection process for SNPs to be analysed

The filtering steps involved are:

1. Removal of SNPs located on sex chromosomes
2. Removal of SNPs within the *PstI* restriction site errors
3. Removal of *fixed* alleles
4. Removal of SNPs not observed in 1996 as these would be likely due to migration not selection
5. Removal of SNPs identified in <75% of both populations
6. Removal of SNPs within the same GBS tag, with identical allele frequencies
7. Removal of SNPs called by `Stacks` to be in separate tags, but which share identical genomic co-ordinates
8. Selection of SNPs within 40kb of a gene

`Stacks` identifies `P` (major) and `Q` (minor) alleles within each population, and an additional processing step was included to ensure the `P` allele in both populations was the SNP determined to be the `P` allele in the 1996 population.

# Data Setup

#### Genome Setup

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

#### Loading SNP Data


```r
sumData <- file.path("..", "data", "batch_3.sumstats.tsv.gz") %>%
  gzfile() %>%
  read_delim(delim = "\t", col_names = TRUE, skip = 3) %>%
  dplyr::select(-contains("Batch"), -contains("Pi"), -contains("Fis")) %>%
  mutate(RefSeq = gsub("^gi\\|.+\\|ref\\|(.+)\\|$", "\\1", Chr),
         Chr = as.character(ncbi2Ens[RefSeq, "Chr"]),
         snpID = paste(`Locus ID`, Col, sep="_"),
         Private = as.logical(Private)) %>%
  filter(`Pop ID` != 3,
         Chr != "X",
         Col > 4) %>%
  dplyr::select(Chr, BP, `Locus ID`, Col, snpID, 
                `Pop ID`, contains("Nuc"),  N, P, 
                contains("Obs"), contains("Exp"), Private)
```

The `stacks` summary output was loaded immediately removing:

- SNPs on Chromosome X
- SNPs within the restriction site
- Data from the Turretfield outgroup, which was not required for this stage of the analysis

This object contained information for 155,455 SNPs across 65,676 GBS tags.

#### Conversion of SNPs to GRanges

The genomic locations of all SNPs were then prepared as a `GRanges` object.


```r
snp2GR <- sumData %>%
  distinct(snpID, .keep_all = TRUE) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE,
                           seqinfo = seqinfo(ensGenes),
                           seqnames.field = "Chr",
                           start.field = "BP",
                           end.field = "BP",
                           ignore.strand = TRUE) %>%
  extract(, c("Locus ID", "snpID")) %>%
  set_names(.$snpID) %>%
  sort
```


# SNP Filtering And Selection

## Uninformative SNPs

- Fixed SNPs were defined as `P > 0.95` in both populations or `P = 0.5` in both populations


```r
fixedSNPs <- sumData %>%
  group_by(snpID) %>%
  summarise(minP = min(P),
            maxP = max(P),
            both50 = all(minP == 0.5, maxP == 0.5)) %>%
  filter(minP > 0.95 | both50) %>%
  extract2("snpID")
```

- SNPs not detected in the 1996 population were identified as these were likely due to migration or some other factor, not selection


```r
uniq2Ora <- filter(sumData, `Pop ID` == 1, P == 1)$snpID
```

- SNPs were also excluded if not identified in > 75% of both populations


```r
minSamp <- sumData %>% 
  group_by(`Pop ID`) %>% 
  summarise(N = max(N)) %>%
  mutate(minN = ceiling(0.75*N)) %>%
  dplyr::select(`Pop ID`, N, minN)
```


| Pop ID |  N | Min Required |
|:-------|---:|-------------:|
| 1      | 55 |           42 |
| 2      | 53 |           40 |

Table: Minimum sample sizes required in each population.


```r
lowSamp <- sumData %>%
  left_join(dplyr::select(minSamp, -N)) %>%
  mutate(keep = N >= minN) %>%
  group_by(snpID) %>%
  summarise(keep = sum(keep) == 2) %>%
  filter(!keep) %>%
  extract2("snpID")
```

These three filtering steps identified SNPs for exclusion as seen in the following table:

| Reason For Removal  | SNPs Marked For Removal       |
|:------------------- | -------------------------------:|
| Fixed Alleles       | 30,336    |
| Not Found in 1996   | 9,153     |
| Identified in < 75% | 77,280      |
| **Total**           | **94,679** |

Table: SNPs identified for removal

## Duplicate SNPs


```r
bestInTag <- sumData %>% 
  filter(!snpID %in% c(fixedSNPs, uniq2Ora, lowSamp)) %>%
  dplyr::select(`Locus ID`, Col, snpID, `Pop ID`, N) %>%
  dcast(`Locus ID` + Col + snpID ~ `Pop ID`) %>% 
  distinct(`Locus ID`, `1`, `2`, .keep_all = TRUE) %>%
  extract2("snpID")
```

- Where multiple SNPs were present in a tag with identical allele frequencies and sample sizes, one was selected at random
- This marked 18,509 SNPs as containing duplicated information


```r
bestDuplicate <- sumData %>% 
  filter(snpID %in% bestInTag) %>%
  group_by(Chr, BP, snpID) %>%
  summarise(N = sum(N)) %>%
  ungroup() %>%
  mutate(`Locus ID` = gsub("([0-9]+)_[0-9]","\\1", snpID)) %>%
  arrange(Chr, BP, N, `Locus ID`) %>%
  distinct(Chr, BP, .keep_all = TRUE) %>%
  extract2("snpID")
```

- Some SNPs were identified at the same genomic location in two loci and these were assumed to be due to differing cut sites, or tags aligned to complementary strands
- If samples sizes were identical, otherwise the SNP with the largest total samples was chosen
- This gave 31,386 candidate SNPs for downstream analysis.


## Selection of SNPs Near Genes


```r
finalSNPs <- snp2GR %>%
  subset(snpID %in% bestDuplicate) %>%
  resize(80001, fix = "center") %>%
  trim() %>%
  subsetByOverlaps(ensGenes) %>%
  names()
```


```r
nGenes <- snp2GR %>%
  extract(finalSNPs) %>%
  resize(80001, fix = "center") %>%
  trim() %>%
  findOverlaps(ensGenes) %>%
  subjectHits() %>%
  unique() %>%
  length()
```

Only SNPs within 40kb of a known gene were then retained, giving the **final total of 20,336 SNPs** for downstream analysis, which potentially provided markers for 10,078 genes.

## Genomic Coverage After SNP Filtering


```r
winSize <- 1e05
windows <- tileGenome(seqinfo(snp2GR)[as.character(1:21),], tilewidth = 1e05)
windowCounts <- countOverlaps(windows, subset(snp2GR, snpID %in% finalSNPs)) 
```

Using 100kb tiled windows across the autosomes, this final list of SNPs gave 16827 windows without a SNP, and 4534 windows with at least one SNP.
This means that 21.2% of windows were sampled to a high enough quality using the GBS approach.
Of the windows with a SNP, the median number of SNPs in each window was 2, with the most densely sampled window containing 31 SNPs.
This is in comparison to the dataset before filtering of low-quality/duplicate SNPs, in which 55% of the 100kb windows contained at least one SNP.


# Correction of P & Q Nucleotides

P and Q alleles are defined separately for each population, and need to be redefined such that these alleles are common across populations.
Define each SNP such that the P (i.e. most common) allele in Population 1 is considered as the reference allele.


```r
sumData %<>%
  filter(snpID %in% finalSNPs) %>%
  split(f = .$snpID) %>%
  mclapply(function(x){
    if (length(unique(x$`P Nuc`)) == 2){ # If there are 2 P alleles
      x[2, c("P Nuc", "Q Nuc")] <- x[2, c("Q Nuc", "P Nuc")]
      x$P[2] <- 1 - x$P[2]
    }
    x
  }, mc.cores = nCores) %>%
  bind_rows
```



```r
gzOut <- file.path("..", "data", "filteredSNPs.tsv.gz") %>%
  gzfile("w")
sumData %>%
  dplyr::select(-Private) %>%
  write_tsv(gzOut)
close(gzOut)
```

This final object was exported as the object `filteredSNPs.tsv.gz`


# Editing of genpop data

The `genpop` file as output by `stacks` was loaded into `R` and converted to genotypes manually.


- First extract all the setup information


```r
gpFile <- file.path("..", "data", "batch_3.genepop.gz") %>% gzfile()
gpData <- readLines(gpFile)
close(gpFile)
popStart <- grep("^pop", gpData) + 1
popEnd <- c(popStart[-1] - 2, length(gpData))
nPops <- length(popStart)
sampleRows <- seq_along(popStart) %>%
  sapply(function(x){seq(popStart[x], popEnd[x])}) %>%
  unlist()
sampleIDs <- gsub("(gc|ora|TF|tf|Y|pt)([0-9ABC]+),.+", "\\1\\2", gpData[sampleRows])
snpIDs <- gpData[2] %>% str_split(",") %>% extract2(1)
```

- Now get the alleles using the existing allele encoding
- Then restrict to the SNPs being analysed


```r
alleles <- gsub("(gc|ora|TF|tf|Y|pt)[0-9ABC]+,\\t", "", gpData[sampleRows]) %>%
  str_split_fixed("\t", n = length(snpIDs)) %>%
  replace(. == "0000", NA) %>%
  set_colnames(snpIDs) %>%
  set_rownames(sampleIDs) %>%
  extract(,finalSNPs)
```

- Convert to minor allele counts


```r
minorAlleleCounts <- alleles %>%
  apply(MARGIN = 2, function(x){
    allNums <- sort(unique(x[!is.na(x)])) # The existing allele combinations
    allIDs <- unique(c(substr(allNums, 1, 2),
                       substr(allNums, 3, 4))) # The existing allele numbers
    # Create all possible genotypes
    gTypes <- c(paste0(allIDs[1], allIDs[1]),
                  paste0(allIDs[1], allIDs[2]),
                  paste0(allIDs[2], allIDs[2]))
    # Find the maximum allele based on homozygotes across all populations
    homs <- vapply(gTypes[c(1, 3)],function(y){sum(grepl(y, x), na.rm = TRUE)}, integer(1))
    P <- which.max(homs)
    # Reverse the genotypes & use as factor levels
    if (P == 2) gTypes <- rev(gTypes)
    as.integer(factor(x, levels = gTypes)) - 1
  }) %>%
  set_rownames(sampleIDs)
```


```r
gzOut <- file.path("..", "data", "minorAlleleCounts.tsv.gz") %>%
  gzfile("w")
minorAlleleCounts %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  write_tsv(gzOut)
close(gzOut)
```

This file was then exported as `minorAlleleCounts.tsv.gz`.
An additional genepop file was created for use with Bayescan and for estimation of effective population sizes.


```r
customGenepopFile <- file.path("..", "data", "filteredSNPs.genepop.gz") 
# The header
write_lines(x = "Created by custom script; Genepop version 4.1.3; Nov 05, 2017", 
           path = gzfile(customGenepopFile))
# Allele names
colnames(alleles) %>% 
  paste(collapse = ",") %>%
  write_lines(gzfile(customGenepopFile), append = TRUE)
# The first population
write_lines("pop", gzfile(customGenepopFile), append = TRUE)
# Allele data for each gc sample
grep("gc", rownames(alleles), value = TRUE) %>%
  lapply(function(x){
    alls <- replace(alleles[x,], is.na(alleles[x,]), "0000") %>%
      paste(collapse = "\t")
    paste(x, alls, sep = ",") %>%
      write_lines(gzfile(customGenepopFile), append = TRUE)
  })
# The next population
write_lines("pop", gzfile(customGenepopFile), append = TRUE)
# Allele data for each ora sample
grep("ora", rownames(alleles), value = TRUE) %>%
  lapply(function(x){
    alls <- replace(alleles[x,], is.na(alleles[x,]), "0000") %>%
      paste(collapse = "\t")
    paste(x, alls, sep = ",") %>%
      write_lines(gzfile(customGenepopFile), append = TRUE)
  })
```


**R version 3.4.3 (2017-11-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_AU.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_AU.UTF-8_, _LC_COLLATE=en_AU.UTF-8_, _LC_MONETARY=en_AU.UTF-8_, _LC_MESSAGES=en_AU.UTF-8_, _LC_PAPER=en_AU.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_AU.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_parallel_, _stats4_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_bindrcpp(v.0.2)_, _scales(v.0.5.0)_, _pander(v.0.6.1)_, _reshape2(v.1.4.3)_, _rtracklayer(v.1.38.3)_, _GenomicRanges(v.1.30.1)_, _GenomeInfoDb(v.1.14.0)_, _IRanges(v.2.12.0)_, _S4Vectors(v.0.16.0)_, _BiocGenerics(v.0.24.0)_, _magrittr(v.1.5)_, _forcats(v.0.2.0)_, _stringr(v.1.2.0)_, _dplyr(v.0.7.4)_, _purrr(v.0.2.4)_, _readr(v.1.1.1)_, _tidyr(v.0.7.2)_, _tibble(v.1.4.2)_, _ggplot2(v.2.2.1)_ and _tidyverse(v.1.2.1)_

**loaded via a namespace (and not attached):** 
_Rcpp(v.0.12.15)_, _lubridate(v.1.7.1)_, _lattice(v.0.20-35)_, _Rsamtools(v.1.30.0)_, _Biostrings(v.2.46.0)_, _assertthat(v.0.2.0)_, _rprojroot(v.1.3-2)_, _digest(v.0.6.14)_, _psych(v.1.7.8)_, _R6(v.2.2.2)_, _cellranger(v.1.1.0)_, _plyr(v.1.8.4)_, _backports(v.1.1.2)_, _evaluate(v.0.10.1)_, _httr(v.1.3.1)_, _pillar(v.1.1.0)_, _zlibbioc(v.1.24.0)_, _rlang(v.0.1.6)_, _lazyeval(v.0.2.1)_, _readxl(v.1.0.0)_, _rstudioapi(v.0.7)_, _Matrix(v.1.2-12)_, _rmarkdown(v.1.8)_, _BiocParallel(v.1.12.0)_, _foreign(v.0.8-69)_, _RCurl(v.1.95-4.10)_, _munsell(v.0.4.3)_, _DelayedArray(v.0.4.1)_, _broom(v.0.4.3)_, _compiler(v.3.4.3)_, _modelr(v.0.1.1)_, _pkgconfig(v.2.0.1)_, _mnormt(v.1.5-5)_, _htmltools(v.0.3.6)_, _SummarizedExperiment(v.1.8.1)_, _GenomeInfoDbData(v.1.0.0)_, _matrixStats(v.0.53.0)_, _XML(v.3.98-1.9)_, _crayon(v.1.3.4)_, _GenomicAlignments(v.1.14.1)_, _bitops(v.1.0-6)_, _grid(v.3.4.3)_, _nlme(v.3.1-131)_, _jsonlite(v.1.5)_, _gtable(v.0.2.0)_, _cli(v.1.0.0)_, _stringi(v.1.1.6)_, _XVector(v.0.18.0)_, _xml2(v.1.2.0)_, _tools(v.3.4.3)_, _Biobase(v.2.38.0)_, _glue(v.1.2.0)_, _hms(v.0.4.1)_, _yaml(v.2.1.16)_, _colorspace(v.1.3-2)_, _rvest(v.0.3.2)_, _knitr(v.1.18)_, _bindr(v.0.1)_ and _haven(v.1.1.1)_


