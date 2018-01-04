---
title: "Analysis Of Genes Associated With Significant SNPs"
author: "Steve Pederson"
date: "04 January, 2018"
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
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
library(GO.db)
library(parallel)
library(magrittr)
library(pander)
library(scales)
library(reshape2)
library(stringr)
```



```r
nCores <- min(12, detectCores() - 1)
theme_set(theme_bw())
evalLongChunks <- FALSE
```

# Introduction

This analysis takes the SNPs found to be associated with the two populations and determines which genes are nearby.
An analysis of the GO terms associated with these SNPS is also undertaken.

## Outline of analysis

This file covers three main areas

1. Preparation of the complete mappings from Ensembl gene IDs to GO terms
2. Finding SNPs within and near genes
3. Identification of enriched GO terms using Fisher's Exact Test on the sets of genes within 2kb and 40kb of the significant SNPs

# Data Setup

## Genes

`GenomicRanges` were formed defining all genes and exons according to the gene models in Build 84 of Ensembl.


```r
gffFile <- file.path("..", "data", "Oryctolagus_cuniculus.OryCun2.0.84.gff3.gz") 
ensGenes <- import.gff(gffFile, feature.type = "gene", sequenceRegionsAsSeqinfo = TRUE) %>%
  sortSeqlevels() 
ensExons <-  import.gff(gffFile, feature.type = "exon", sequenceRegionsAsSeqinfo = TRUE) %>%
  sortSeqlevels()
```

## GO Terms

### Mapping Genes to GO Terms

- The set of terms defining the root nodes was initially formed


```r
rootGO <- list(
  BP = "GO:0008150",
  CC = "GO:0005575",
  MF = "GO:0003674"
  )
```

- The complete set of GO IDs mapping to these genes was then obtained from biomaRt.


```r
mart <- useEnsembl("ensembl", dataset ="ocuniculus_gene_ensembl", version = "84")
ens2Go <- getBM(attributes = c("ensembl_gene_id","go_id", "namespace_1003"), 
                filters = "ensembl_gene_id", 
                values = ensGenes$gene_id, 
                mart = mart) %>%
  filter(go_id != "") %>% 
  dplyr::rename(ontology = namespace_1003) %>%
  mutate(ontology = c(biological_process = "BP", 
                      cellular_component = "CC", 
                      molecular_function = "MF")[ontology])
```

- GO Terms in the biomaRt results which failed to map to the database were removed.


```r
ens2Go %<>% split(f = .$ontology)
ens2Go$BP %<>% filter(go_id %in% get(rootGO$BP, GOBPOFFSPRING))
ens2Go$CC %<>% filter(go_id %in% get(rootGO$CC, GOCCOFFSPRING))
ens2Go$MF %<>% filter(go_id %in% get(rootGO$MF, GOMFOFFSPRING))
ens2Go %<>% bind_rows()
```
- The file was then saved to disk and became the reference database for subsequent analyses.


```r
gzOut <- file.path("..", "data", "ens2GO.biomaRt.84.tsv.gz") %>% gzfile("w")
write_delim(ens2Go, gzOut, delim ="\t")
close(gzOut)
```


```r
ens2Go <- file.path("..", "data","ens2GO.biomaRt.84.tsv.gz") %>%
  gzfile() %>%
  read_tsv()
```

### Building the Complete Set of Mappings

- Paths were formed mapping all GO Terms obtained from `biomaRt` back to the root nodes.


```r
allGOAncestor <- list(
  BP = filter(ens2Go, ontology == "BP") %>%
    extract2("go_id") %>%
    unique() %>%
    sapply(get, GOBPANCESTOR, simplify = FALSE),
  CC = filter(ens2Go, ontology == "CC") %>%
    extract2("go_id") %>%
    unique() %>%
    sapply(get, GOCCANCESTOR, simplify = FALSE),
  MF = filter(ens2Go, ontology == "MF") %>%
    extract2("go_id") %>%
    unique() %>%
    sapply(get, GOMFANCESTOR, simplify = FALSE)
)
```

- Each gene was then mapped to the complete set of GO terms and saved in an external table.


```r
ALLGO2ENS <- ens2Go %>%
  split(f = .$ensembl_gene_id) %>%
  mclapply(function(x){
    g <- unique(x$ensembl_gene_id)
    BP <- allGOAncestor$BP[filter(x, ontology == "BP")$go_id] %>%
      unlist() %>%
      c(filter(x, ontology == "BP")$go_id) %>%      
      unique()
    CC <- allGOAncestor$CC[filter(x, ontology == "CC")$go_id] %>%
      unlist() %>%
      c(filter(x, ontology == "CC")$go_id) %>%      
      unique()
    MF <- allGOAncestor$MF[filter(x, ontology == "MF")$go_id] %>%
      unlist() %>%
      c(filter(x, ontology == "MF")$go_id) %>%
      unique() 
    data_frame(ensembl_gene_id = g, go_id = c(BP, CC, MF))
  },mc.cores = nCores) %>%
  bind_rows() %>%
  filter(go_id != "all") %>%
  mutate(ontology = Ontology(go_id))
```


```r
gzOut <- file.path("..", "data", "ALLGO2ENS.tsv.gz") %>% gzfile("w")
write_delim(ALLGO2ENS, gzOut, delim = "\t")
close(gzOut)
```


```r
ALLGO2ENS <- file.path("..", "data", "ALLGO2ENS.tsv.gz") %>%
  gzfile() %>%
  read_tsv()
singleGeneGO <- ALLGO2ENS %>% 
  group_by(go_id) %>% 
  tally %>% 
  filter(n == 1) %>%
  extract2("go_id")
```

- This object contained all 1,264,162 possible gene to GO mappings for this build of Ensembl.
- After building the database once, this object then became the reference object for all downstream analysis.
- A subset of 4237 GO terms which mapped to only a single gene were also noted for exclusion from testing, as no significant results would be possible from these terms.

### Defining GO Term Levels

GO Terms were then defined in terms of their minimum distance to the Root Nodes, down to those requiring up to 3 steps.


```r
firstLevelGO <- list(
  BP = c(rootGO$BP, get(rootGO$BP, GOBPCHILDREN)),
  CC = c(rootGO$CC, get(rootGO$CC, GOCCCHILDREN)),
  MF = c(rootGO$MF, get(rootGO$MF, GOMFCHILDREN))
) %>%
  lapply(unique)
secondLevelGO <- list(
  BP = c(firstLevelGO$BP, lapply(firstLevelGO$BP, get, GOBPCHILDREN)),
  CC = c(firstLevelGO$CC, lapply(firstLevelGO$CC, get, GOCCCHILDREN)),
  MF = c(firstLevelGO$MF, lapply(firstLevelGO$MF, get, GOMFCHILDREN))
) %>%
  lapply(unlist) %>%
  lapply(unique) %>%
  lapply(function(x){x[!is.na(x)]})
thirdLevelGO <- list(
  BP = c(secondLevelGO$BP, lapply(secondLevelGO$BP, get, GOBPCHILDREN)),
  CC = c(secondLevelGO$CC, lapply(secondLevelGO$CC, get, GOCCCHILDREN)),
  MF = c(secondLevelGO$MF, lapply(secondLevelGO$MF, get, GOMFCHILDREN))
) %>%
  mclapply(unlist, mc.cores = nCores) %>%
  mclapply(unique, mc.cores = nCores) %>%
  mclapply(function(x){x[!is.na(x)]}, mc.cores = nCores)
```

## SNPs


```r
resultFiles <- file.path("..", "results", c("flkResults.tsv", 
                                            "alleleResults.tsv",
                                            "genotypeResults.tsv")) %>%
  set_names(c("flk", "alleles", "genotypes"))
results <- resultFiles %>%
  lapply(read_tsv) 
```


```r
resultsGR <- full_join(
  results$flk %>%
    dplyr::select(snpID, Chr, BP, SNP, `Gum Creek (1996)`, `Oraparinna (2012)`,
                  FLK_p = F.LK.p.val, FLK_FDR = FDR),
  results$alleles %>%
    dplyr::select(snpID, Chr, BP, alleles_p = p, alleles_FDR = FDR)
) %>%
  full_join(
  results$genotypes %>%
    dplyr::select(snpID, Chr, BP, genotypes_p = p, genotypes_FDR = FDR)
  ) %>%
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE, 
    ignore.strand = TRUE,
    seqnames.field = "Chr", 
    start.field = "BP", 
    end.field = "BP", 
    seqinfo = seqinfo(ensGenes)) %>%
  sort()
```


```r
resultsGR$Sig <- with(resultsGR, FLK_FDR < 0.05 | genotypes_FDR < 0.1)
```


# Relationship of SNPs to Genes

## SNPs within genes


```r
snpInGenes <- resultsGR %>% 
  subset(Sig) %>%
  findOverlaps(ensGenes)
snpInExons <- resultsGR %>% 
  subset(Sig) %>%
  findOverlaps(ensExons)
```

A total of 28 significant SNPs were found to be within the start and end points of genes.
4 of these were within the coding regions of these genes.


| Gene ID            | Name     |      Chr |          BP | snpID      |
|:-------------------|:---------|---------:|------------:|:-----------|
| ENSOCUG00000003383 | MGAT4A   |        2 |  93,561,889 | 214772_69  |
| ENSOCUG00000017054 | TEX33    |        4 |  84,940,214 | 167107_60  |
| ENSOCUG00000017054 | TEX33    |        4 |  84,940,235 | 167108_14  |
| ENSOCUG00000011060 | MYO1B    |        7 | 131,862,327 | 147965_18  |
| ENSOCUG00000024842 | CTIF     |        9 |  89,862,459 | 211772_50  |
| ENSOCUG00000024842 | CTIF     |        9 |  89,862,481 | 211772_28  |
| ENSOCUG00000024842 | CTIF     |        9 |  89,862,497 | 211772_12  |
| ENSOCUG00000015851 | DYM      |        9 |  90,037,571 | 137949_40  |
| ENSOCUG00000015851 | DYM      |        9 |  90,037,630 | 137951_16  |
| ENSOCUG00000000804 | SCRN1    |       10 |  14,373,457 | 127156_20  |
| ENSOCUG00000000804 | SCRN1    |       10 |  14,373,501 | 127157_45  |
| ENSOCUG00000012252 | DPP6     |       13 |  15,938,245 | 107422_35  |
| ENSOCUG00000007988 | MUC4     |       14 |  92,793,012 | 101831_18* |
| ENSOCUG00000016518 | KIAA1217 |       16 |     735,716 | 87819_88   |
| ENSOCUG00000007461 | NAV1     |       16 |  69,850,541 | 87407_11   |
| ENSOCUG00000015494 | NIN      |       17 |  69,838,525 | 80772_39   |
| ENSOCUG00000011533 | RHOBTB1  |       18 |  25,311,139 | 72765_47*  |
| ENSOCUG00000007069 | FGFR2    |       18 |  68,517,919 | 208110_71  |
| ENSOCUG00000007069 | FGFR2    |       18 |  68,517,942 | 208111_21  |
| ENSOCUG00000013175 | FBLN5    |       20 |  14,022,385 | 64526_24*  |
| ENSOCUG00000001722 | KSR2     |       21 |  14,643,101 | 62998_77   |
| ENSOCUG00000021564 |          | GL018704 |   4,853,505 | 53831_42*  |
| ENSOCUG00000009435 | XKR5     | GL018713 |     365,938 | 50206_34   |
| ENSOCUG00000000596 | PLA2G16  | GL018717 |     314,346 | 48659_40   |
| ENSOCUG00000017106 | KAZN     | GL018739 |      75,851 | 41475_63   |
| ENSOCUG00000010092 | CPNE9    | GL018802 |     439,482 | 204810_79  |
| ENSOCUG00000027400 | PLAUR    | GL018881 |      88,327 | 233206_29  |
| ENSOCUG00000023502 |          | GL019002 |      95,028 | 14915_82   |

Table: Genes with significant SNPs located between their start and end positions. SNPs within exons are indicated with an additional asterisk.

- It was also noted that proteins produced from `MYO1B` and `USP46`are known interacting proteins in human biology

## Genes within 20kb


```r
allGenesIn20kb <- subsetByOverlaps(
  ensGenes, 
  resultsGR %>%
    resize(width = 40001,fix = "center") %>%
    trim()
    )
sigGenesIn20kb <- subsetByOverlaps(
  ensGenes, 
  resultsGR %>%
    subset(Sig)%>%
    resize(width = 40001,fix = "center") %>%
    trim()
    )
```

As a an initial conservative approach, the set of genes within 20kb of a SNP was investigated.
This produced a list of 7365 genes in total, of which 43 overlapped SNPs of interest.



```r
go2Test20kb <- ALLGO2ENS %>%
  filter(ensembl_gene_id %in% sigGenesIn20kb$gene_id,
         !go_id %in% unlist(thirdLevelGO)) %>%
  extract2("go_id") %>%
  unique %>%
  setdiff(singleGeneGO)
```

A set of 522 GO terms mapping to more than one gene, and with at least 4 steps back to the ontology roots, were assigned to the genes within 20kb of the 46 candidate SNPs.
These were then tested for enrichment using the set of genes within 20kb of the remaining 17892 SNPs

Due to the closely related nature of some SNPs in this dataset the possibility of two SNPs being within 20kb of the same genes was extremely high, and this approach would ensure that each gene only appears once despite the possibility of mapping to multiple SNPs within the dataset.

Fisher's Exact Test was then used to test for enrichment by counting the genes mapped to each GO term within the set of significant, and non-significant SNPs, and comparing to the genes not mapped to each GO term in both sets of SNPs.
All p values were then adjusted using Benjamini-Hochberg's FDR.



```r
nSigGenes <- length(sigGenesIn20kb)
nNotSigGenes <- length(allGenesIn20kb) - nSigGenes
goResults20kb <- go2Test20kb %>%
  lapply(function(x){
    # Form a matrix of the counts
    mat <- filter(ALLGO2ENS, go_id == x) %>%
      mutate(Sig = ensembl_gene_id %in% sigGenesIn20kb$gene_id) %>%
      acast(Sig~., fun.aggregate = length, value.var = "Sig")
    # Ensure it is a matrix, then set row/colnames
    if (length(mat) == 1) mat <- c(0, mat)
    mat <- cbind(mat, withoutGO = c(nNotSigGenes, nSigGenes) - mat)
    colnames(mat)[1] <- "withGO"
    rownames(mat) <- c("notNearSNP", "nearSNP")
    # Fisher's exact test
    ft <- fisher.test(mat)
    data_frame(go_id = x,
               expected = nSigGenes * mat["notNearSNP", "withGO"] / nNotSigGenes,
               observed = mat["nearSNP", "withGO"],
               p = ft$p.value)
  }) %>%
  bind_rows() %>%
  mutate(FDR = p.adjust(p, "fdr")) %>%
  arrange(p)
```

A total of 6 GO terms were considered as enriched using the criteria of an FDR-adjusted p-value < 0.05 and with observed numbers greater than that predicted by the ratio in the non-significant SNP genes.


| GO ID      | Term                                   | Ont | Observed | Expected |         p |      FDR |
|:-----------|:---------------------------------------|:----|---------:|---------:|----------:|---------:|
| GO:0010359 | regulation of anion channel activity   | BP  |        2 | 0.005873 | 9.953e-05 | 0.003996 |
| GO:0004792 | thiosulfate sulfurtransferase activity | MF  |        2 | 0.005873 | 9.953e-05 | 0.003996 |
| GO:0071294 | cellular response to zinc ion          | BP  |        2 |  0.01762 | 0.0003293 |  0.01011 |
| GO:0010043 | response to zinc ion                   | BP  |        2 |  0.04111 |  0.001168 |  0.02439 |
| GO:0016783 | sulfurtransferase activity             | MF  |        2 |  0.04698 |  0.001455 |   0.0292 |
| GO:0048471 | perinuclear region of cytoplasm        | CC  |        7 |    1.791 |   0.00198 |  0.03691 |

Table: Gene Ontologies considered as enriched amongst the set of 43 genes within 20kb of the significant SNPs. The number of genes matching each term is given in the 'Observed' column.

## Genes within 40kb


```r
allGenesIn40kb <- subsetByOverlaps(
  ensGenes, 
  resultsGR %>%
    resize(width = 80001,fix = "center") %>%
    trim()
    )
sigGenesIn40kb <- subsetByOverlaps(
  ensGenes, 
  resultsGR %>%
    subset(Sig)%>%
    resize(width = 80001,fix = "center") %>%
    trim()
    )
```

The same approach was then repeated for genes within 40kb of each SNP.
This produced a list of 9636 genes in total, of which 68 overlapped SNPs of interest.



```r
go2Test40kb <- ALLGO2ENS %>%
  filter(ensembl_gene_id %in% sigGenesIn40kb$gene_id,
         !go_id %in% unlist(thirdLevelGO)) %>%
  extract2("go_id") %>%
  unique %>%
  setdiff(singleGeneGO)
```

A set of 754 GO terms mapping to more than one gene, and with at least 4 steps back to the ontology roots were assigned to the genes within 40kb of the 46 candidate SNPs.
These were then tested for enrichment using the set of genes within 40kb of the remaining 17892 SNPs



```r
nSigGenes <- length(sigGenesIn40kb)
nNotSigGenes <- length(allGenesIn40kb) - nSigGenes
goResults40kb <- go2Test40kb %>%
  lapply(function(x){
    # Form a matrix of the counts
    mat <- filter(ALLGO2ENS, go_id == x) %>%
      mutate(Sig = ensembl_gene_id %in% sigGenesIn40kb$gene_id) %>%
      acast(Sig~., fun.aggregate = length, value.var = "Sig")
    # Ensure it is a matrix, then set row/colnames
    if (length(mat) == 1) mat <- c(0, mat)
    mat <- cbind(mat, withoutGO = c(nNotSigGenes, nSigGenes) - mat)
    colnames(mat)[1] <- "withGO"
    rownames(mat) <- c("notNearSNP", "nearSNP")
    # Fisher's exact test
    ft <- fisher.test(mat)
    data_frame(go_id = x,
               expected = nSigGenes * mat["notNearSNP", "withGO"] / nNotSigGenes,
               observed = mat["nearSNP", "withGO"],
               p = ft$p.value)
  }) %>%
  bind_rows() %>%
  mutate(FDR = p.adjust(p, "fdr")) %>%
  arrange(p)
```

A total of 7 GO terms were considered as enriched using the criteria of an FDR-adjusted p-value < 0.05 and with observed numbers greater than that predicted by the ratio in the non-significant SNP genes.


| GO ID      | Term                                   | Ont | Observed | Expected |         p |      FDR |
|:-----------|:---------------------------------------|:----|---------:|---------:|----------:|---------:|
| GO:0071294 | cellular response to zinc ion          | BP  |        3 |  0.01421 | 3.328e-06 | 0.002509 |
| GO:0010043 | response to zinc ion                   | BP  |        3 |  0.04264 | 2.739e-05 |  0.01033 |
| GO:0010359 | regulation of anion channel activity   | BP  |        2 | 0.007107 | 0.0001465 |   0.0221 |
| GO:0004792 | thiosulfate sulfurtransferase activity | MF  |        2 | 0.007107 | 0.0001465 |   0.0221 |
| GO:0015793 | glycerol transport                     | BP  |        2 | 0.007107 | 0.0001465 |   0.0221 |
| GO:0015791 | polyol transport                       | BP  |        2 |  0.01421 | 0.0002918 |  0.03666 |
| GO:0070295 | renal water absorption                 | BP  |        2 |  0.02132 |  0.000484 |  0.04843 |

Table: Gene Ontologies considered as enriched amongst the set of 68 genes within 40kb of the significant SNPs. The number of genes matching each term is given in the 'Observed' column.

Similar sets of terms were detected at both 20kb and 40kb.
The appearance of terms connected to the Zinc ion may be of note as the link between zinc and clearance of HCV has recently been established, via the IFN-&#947; pathway

### Export of Results

The set of genes within 40kb of each SNP of interest were then exported.


```r
hitsIn40kb <- findOverlaps(
  resultsGR %>%
    subset(Sig)%>%
    resize(width = 80001,fix = "center") %>%
    trim(),
    ensGenes
    )
```


```r
tsvOut <- file.path("..", "results", "GenesIn40kb.tsv")
list(
  resultsGR %>%
    subset(Sig) %>%
    extract(queryHits(hitsIn40kb)) %>%
    as.data.frame() %>%
    dplyr::select(snpID, Chr = seqnames, BP = start, genotypes_p,  FLK_p) ,
  ensGenes %>%
    extract(subjectHits(hitsIn40kb)) %>%
    as.data.frame() %>%
    dplyr::select(GeneStart = start, GeneEnd = end, strand, NearGene = gene_id, GeneName = Name) 
) %>%
  bind_cols %>%
  as_data_frame() %>%
    rowwise()%>%
  mutate(LocusID = gsub("([0-9]*)_[0-9]*", "\\1", snpID),
         dist2Gene = if_else(BP > GeneStart && BP < GeneEnd, 0L, min(abs(c(BP - GeneStart, BP - GeneEnd)))),
         GeneWidth = abs(GeneStart - GeneEnd)) %>%
  dplyr::select(LocusID, snpID, Chr, BP, NearGene, GeneName, GeneStrand = strand, GeneStart, GeneEnd, GeneWidth, dist2Gene, genotypes_p, FLK_p) %>%
  mutate(minP = min(genotypes_p, FLK_p)) %>%
  arrange(minP) %>%
  dplyr::select(-minP) %>%
  as.data.frame() %>%
  write_tsv(tsvOut)
```


```r
pander(sessionInfo()) 
```

**R version 3.4.3 (2017-11-30)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_AU.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_AU.UTF-8_, _LC_COLLATE=en_AU.UTF-8_, _LC_MONETARY=en_AU.UTF-8_, _LC_MESSAGES=en_AU.UTF-8_, _LC_PAPER=en_AU.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_AU.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_parallel_, _stats4_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_bindrcpp(v.0.2)_, _reshape2(v.1.4.2)_, _scales(v.0.5.0)_, _pander(v.0.6.1)_, _magrittr(v.1.5)_, _GO.db(v.3.5.0)_, _AnnotationDbi(v.1.40.0)_, _Biobase(v.2.38.0)_, _biomaRt(v.2.34.0)_, _rtracklayer(v.1.38.2)_, _GenomicRanges(v.1.30.0)_, _GenomeInfoDb(v.1.14.0)_, _IRanges(v.2.12.0)_, _S4Vectors(v.0.16.0)_, _BiocGenerics(v.0.24.0)_, _forcats(v.0.2.0)_, _stringr(v.1.2.0)_, _dplyr(v.0.7.4)_, _purrr(v.0.2.4)_, _readr(v.1.1.1)_, _tidyr(v.0.7.2)_, _tibble(v.1.3.4)_, _ggplot2(v.2.2.1)_ and _tidyverse(v.1.2.1)_

**loaded via a namespace (and not attached):** 
_httr(v.1.3.1)_, _bit64(v.0.9-7)_, _jsonlite(v.1.5)_, _modelr(v.0.1.1)_, _assertthat(v.0.2.0)_, _blob(v.1.1.0)_, _GenomeInfoDbData(v.0.99.1)_, _cellranger(v.1.1.0)_, _Rsamtools(v.1.30.0)_, _progress(v.1.1.2)_, _yaml(v.2.1.15)_, _RSQLite(v.2.0)_, _backports(v.1.1.1)_, _lattice(v.0.20-35)_, _glue(v.1.2.0)_, _digest(v.0.6.12)_, _XVector(v.0.18.0)_, _rvest(v.0.3.2)_, _colorspace(v.1.3-2)_, _htmltools(v.0.3.6)_, _Matrix(v.1.2-12)_, _plyr(v.1.8.4)_, _psych(v.1.7.8)_, _XML(v.3.98-1.9)_, _pkgconfig(v.2.0.1)_, _broom(v.0.4.3)_, _haven(v.1.1.0)_, _zlibbioc(v.1.24.0)_, _BiocParallel(v.1.12.0)_, _SummarizedExperiment(v.1.8.0)_, _lazyeval(v.0.2.1)_, _cli(v.1.0.0)_, _mnormt(v.1.5-5)_, _crayon(v.1.3.4)_, _readxl(v.1.0.0)_, _memoise(v.1.1.0)_, _evaluate(v.0.10.1)_, _nlme(v.3.1-131)_, _xml2(v.1.1.1)_, _foreign(v.0.8-69)_, _prettyunits(v.1.0.2)_, _tools(v.3.4.3)_, _hms(v.0.4.0)_, _matrixStats(v.0.52.2)_, _munsell(v.0.4.3)_, _DelayedArray(v.0.4.1)_, _Biostrings(v.2.46.0)_, _compiler(v.3.4.3)_, _rlang(v.0.1.4)_, _grid(v.3.4.3)_, _RCurl(v.1.95-4.8)_, _rstudioapi(v.0.7)_, _bitops(v.1.0-6)_, _rmarkdown(v.1.8)_, _gtable(v.0.2.0)_, _DBI(v.0.7)_, _R6(v.2.2.2)_, _GenomicAlignments(v.1.14.1)_, _lubridate(v.1.7.1)_, _knitr(v.1.17)_, _bit(v.1.1-12)_, _bindr(v.0.1)_, _rprojroot(v.1.2)_, _stringi(v.1.1.6)_ and _Rcpp(v.0.12.14)_
