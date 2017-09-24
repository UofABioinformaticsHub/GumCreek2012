# Analysis Of Genes Associated With Significant SNPs
Steve Pederson  
`r format(Sys.Date(), "%b %d, %Y")`  




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
```



```r
nCores <- min(12, detectCores() - 1)
theme_set(theme_bw())
evalLongChunks <- FALSE
```

# Introduction

This analysis takes the SNPs found to be associated with the two populations and determines which genes are nearby.
An analysis of the GO terms associated with these SNPS is also undertaken.

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
```

- This object contained all 1,289,778 possible gene to GO mappings for this build of Ensembl.
- After building the database once, this object then became the reference object for all downstream analysis.

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

Table: Genes with significant SNPs located between their start and end positions. SNPs within exons are indicated with anadditional asterisk.

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

As a an intial conservative approach, the set of genes within 20kb of a SNP was investigated.
This produced a list of 7407 genes in total, of which 46 overlapped SNPs of interest.

The set of GO terms assigned to these genes in `biomaRt` were extended back to the ontology roots.

