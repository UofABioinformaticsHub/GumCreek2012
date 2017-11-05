# Analysis Of Genes Associated With Significant SNPs
Steve Pederson  
27th September 2017  




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

- This object contained all 1,289,778 possible gene to GO mappings for this build of Ensembl.
- After building the database once, this object then became the reference object for all downstream analysis.
- A subset of 4243 GO terms which mapped to only a single gene were also noted for exclusion from testing, as no significant results would be possible from these terms.

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
resultFiles <- file.path("..", "results", c("alleleResults.tsv", 
                                            "genotypeResults.tsv")) %>%
  set_names(c("alleles", "genotypes"))
results <- resultFiles %>%
  lapply(read_tsv) 
```


```r
resultsGR <- do.call(
  "full_join", 
  names(results) %>%
    lapply(function(x){
      names(results[[x]]) <- gsub("^(p|FDR)$", paste(x, "\\1", sep = "_"), names(results[[x]]))
      dplyr::select(results[[x]], -adjP)
    })) %>%
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
resultsGR$Sig <- resultsGR$alleles_FDR < 0.1 | resultsGR$genotypes_FDR < 0.1
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

A total of 27 significant SNPs were found to be within the start and end points of genes.
5 of these were within the coding regions of these genes.


| Gene ID            | Name    |      Chr |          BP | snpID      |
|:-------------------|:--------|---------:|------------:|:-----------|
| ENSOCUG00000017054 | TEX33   |        4 |  84,940,214 | 167107_60  |
| ENSOCUG00000017054 | TEX33   |        4 |  84,940,235 | 167108_14  |
| ENSOCUG00000005399 | EPO     |        6 |  27,103,573 | 157156_58* |
| ENSOCUG00000005399 | EPO     |        6 |  27,103,576 | 157157_45* |
| ENSOCUG00000001407 | MAGI2   |        7 |  38,250,131 | 151791_54  |
| ENSOCUG00000011060 | MYO1B   |        7 | 131,862,327 | 147965_18  |
| ENSOCUG00000005605 | LDLRAD4 |        9 |  48,543,229 | 134595_50  |
| ENSOCUG00000000804 | SCRN1   |       10 |  14,373,457 | 127156_20  |
| ENSOCUG00000000804 | SCRN1   |       10 |  14,373,501 | 127157_45  |
| ENSOCUG00000007988 | MUC4    |       14 |  92,793,012 | 101831_18* |
| ENSOCUG00000007461 | NAV1    |       16 |  69,850,541 | 87407_11   |
| ENSOCUG00000007461 | NAV1    |       16 |  69,850,609 | 208806_6   |
| ENSOCUG00000015494 | NIN     |       17 |  69,838,525 | 80772_39   |
| ENSOCUG00000011533 | RHOBTB1 |       18 |  25,311,139 | 72765_47*  |
| ENSOCUG00000008478 | STXBP4  |       19 |  32,825,478 | 68946_78   |
| ENSOCUG00000009435 | XKR5    | GL018713 |     365,938 | 50206_34   |
| ENSOCUG00000000596 | PLA2G16 | GL018717 |     314,346 | 48659_40   |
| ENSOCUG00000017106 | KAZN    | GL018739 |      75,851 | 41475_63   |
| ENSOCUG00000017031 | USP46   | GL018754 |   1,365,061 | 37345_16   |
| ENSOCUG00000010092 | CPNE9   | GL018802 |     439,482 | 204810_79  |
| ENSOCUG00000005387 | XRCC1   | GL018881 |      22,826 | 21889_62   |
| ENSOCUG00000005387 | XRCC1   | GL018881 |      24,194 | 21896_47   |
| ENSOCUG00000005387 | XRCC1   | GL018881 |      24,230 | 21896_11   |
| ENSOCUG00000027400 | PLAUR   | GL018881 |      88,327 | 233206_29  |

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
This produced a list of 7407 genes in total, of which 46 overlapped SNPs of interest.



```r
go2Test20kb <- ALLGO2ENS %>%
  filter(ensembl_gene_id %in% sigGenesIn20kb$gene_id,
         !go_id %in% unlist(thirdLevelGO)) %>%
  extract2("go_id") %>%
  unique %>%
  setdiff(singleGeneGO)
```

A set of 639 GO terms mapping to more than one gene, and with at least 4 steps back to the ontology roots, were assigned to the genes within 20kb of the 44 candidate SNPs.
These were then tested for enrichment using the set of genes within 20kb of the remaining 18349 SNPs

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

A total of 4 GO terms were considered as enriched using the criteria of an FDR-adjusted p-value < 0.05 and with observed numbers greater than that predicted by the ratio in the non-significant SNP genes.


| GO ID      | Term                                   | Ont | Observed | Expected |         p |      FDR |
|:-----------|:---------------------------------------|:----|---------:|---------:|----------:|---------:|
| GO:0004792 | thiosulfate sulfurtransferase activity | MF  |        2 | 0.006249 | 0.0001128 | 0.006371 |
| GO:0071294 | cellular response to zinc ion          | BP  |        2 |  0.01875 | 0.0003729 |  0.01402 |
| GO:0016783 | sulfurtransferase activity             | MF  |        2 |  0.04374 |  0.001321 |  0.03518 |
| GO:0010043 | response to zinc ion                   | BP  |        2 |  0.04374 |  0.001321 |  0.03518 |

Table: Gene Ontologies considered as enriched amongst the set of 46 genes within 20kb of the significant SNPs. The number of genes matching each term is given in the 'Observed' column.

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
This produced a list of 9688 genes in total, of which 71 overlapped SNPs of interest.



```r
go2Test40kb <- ALLGO2ENS %>%
  filter(ensembl_gene_id %in% sigGenesIn40kb$gene_id,
         !go_id %in% unlist(thirdLevelGO)) %>%
  extract2("go_id") %>%
  unique %>%
  setdiff(singleGeneGO)
```

A set of 867 GO terms mapping to more than one gene, and with at least 4 steps back to the ontology roots were assigned to the genes within 40kb of the 44 candidate SNPs.
These were then tested for enrichment using the set of genes within 40kb of the remaining 18349 SNPs



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

A total of 3 GO terms were considered as enriched using the criteria of an FDR-adjusted p-value < 0.05 and with observed numbers greater than that predicted by the ratio in the non-significant SNP genes.


| GO ID      | Term                                   | Ont | Observed | Expected |         p |      FDR |
|:-----------|:---------------------------------------|:----|---------:|---------:|----------:|---------:|
| GO:0071294 | cellular response to zinc ion          | BP  |        3 |  0.01477 | 3.733e-06 | 0.003236 |
| GO:0010043 | response to zinc ion                   | BP  |        3 |   0.0443 |  3.07e-05 | 0.008873 |
| GO:0004792 | thiosulfate sulfurtransferase activity | MF  |        2 | 0.007383 | 0.0001581 |  0.03427 |

Table: Gene Ontologies considered as enriched amongst the set of 71 genes within 40kb of the significant SNPs. The number of genes matching each term is given in the 'Observed' column.

Similar sets of terms were detected at both 20kb and 40kb.
The appearance of terms connected to the Zinc ion may be of note as the link between zinc and clearance of HCV has recently been established, via the IFN-&#947; pathway

### Export of Results

The set of genes within 40kb of each SNP of nterest were then exported.


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
    dplyr::select(snpID, Chr = seqnames, BP = start, genotypes_p,  alleles_p) ,
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
  dplyr::select(LocusID, snpID, Chr, BP, NearGene, GeneName, GeneStrand = strand, GeneStart, GeneEnd, GeneWidth, dist2Gene, genotypes_p, alleles_p) %>%
  mutate(minP = min(genotypes_p, alleles_p)) %>%
  arrange(minP) %>%
  dplyr::select(-minP) %>%
  as.data.frame() %>%
  write_tsv(tsvOut)
```


```r
pander(sessionInfo()) 
```

**R version 3.4.2 (2017-09-28)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_AU.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_AU.UTF-8_, _LC_COLLATE=en_AU.UTF-8_, _LC_MONETARY=en_AU.UTF-8_, _LC_MESSAGES=en_AU.UTF-8_, _LC_PAPER=en_AU.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_AU.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_parallel_, _stats4_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_bindrcpp(v.0.2)_, _stringr(v.1.2.0)_, _reshape2(v.1.4.2)_, _scales(v.0.5.0)_, _pander(v.0.6.1)_, _magrittr(v.1.5)_, _GO.db(v.3.4.1)_, _AnnotationDbi(v.1.38.2)_, _Biobase(v.2.36.2)_, _biomaRt(v.2.32.1)_, _rtracklayer(v.1.36.6)_, _GenomicRanges(v.1.28.6)_, _GenomeInfoDb(v.1.12.3)_, _IRanges(v.2.10.5)_, _S4Vectors(v.0.14.7)_, _BiocGenerics(v.0.22.1)_, _dplyr(v.0.7.4)_, _purrr(v.0.2.4)_, _readr(v.1.1.1)_, _tidyr(v.0.7.2)_, _tibble(v.1.3.4)_, _ggplot2(v.2.2.1)_ and _tidyverse(v.1.1.1)_

**loaded via a namespace (and not attached):** 
_httr(v.1.3.1)_, _bit64(v.0.9-7)_, _jsonlite(v.1.5)_, _modelr(v.0.1.1)_, _assertthat(v.0.2.0)_, _blob(v.1.1.0)_, _GenomeInfoDbData(v.0.99.0)_, _cellranger(v.1.1.0)_, _Rsamtools(v.1.28.0)_, _yaml(v.2.1.14)_, _RSQLite(v.2.0)_, _backports(v.1.1.1)_, _lattice(v.0.20-35)_, _glue(v.1.1.1)_, _digest(v.0.6.12)_, _XVector(v.0.16.0)_, _rvest(v.0.3.2)_, _colorspace(v.1.3-2)_, _htmltools(v.0.3.6)_, _Matrix(v.1.2-11)_, _plyr(v.1.8.4)_, _psych(v.1.7.8)_, _XML(v.3.98-1.9)_, _pkgconfig(v.2.0.1)_, _broom(v.0.4.2)_, _haven(v.1.1.0)_, _zlibbioc(v.1.22.0)_, _BiocParallel(v.1.10.1)_, _SummarizedExperiment(v.1.6.5)_, _lazyeval(v.0.2.0)_, _mnormt(v.1.5-5)_, _readxl(v.1.0.0)_, _memoise(v.1.1.0)_, _evaluate(v.0.10.1)_, _nlme(v.3.1-131)_, _forcats(v.0.2.0)_, _xml2(v.1.1.1)_, _foreign(v.0.8-69)_, _tools(v.3.4.2)_, _hms(v.0.3)_, _matrixStats(v.0.52.2)_, _munsell(v.0.4.3)_, _DelayedArray(v.0.2.7)_, _Biostrings(v.2.44.2)_, _compiler(v.3.4.2)_, _rlang(v.0.1.2)_, _grid(v.3.4.2)_, _RCurl(v.1.95-4.8)_, _bitops(v.1.0-6)_, _rmarkdown(v.1.6)_, _gtable(v.0.2.0)_, _DBI(v.0.7)_, _R6(v.2.2.2)_, _GenomicAlignments(v.1.12.2)_, _lubridate(v.1.6.0)_, _knitr(v.1.17)_, _bit(v.1.1-12)_, _bindr(v.0.1)_, _rprojroot(v.1.2)_, _stringi(v.1.1.5)_ and _Rcpp(v.0.12.13)_
