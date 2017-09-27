# Simulation of Genetic Drift
Steve Pederson  
28 September 2017  




```r
library(pander)
library(magrittr)
library(dplyr)
library(scales)
```


```r
logit <- binomial()$linkfun
inv.logit <- binomial()$linkinv
```



```r
nGen <- 16
nSim <- 5000
litterSize <- 30
migRate <- 0.15
alpha <- 0.001
```

# Introduction

This analysis takes the filtered SNPs from under analysis and simulates genetic drift under no selective pressure, in order to compare detected changes in allele frequencies to the range of values predicted by the drift model.

## Outline of Drift Simulation

This analysis uses a custom set of scripts collected as an R package and available from https://github.com/steveped/driftSim.
This is able to be installed using `devtools::install_github("steveped/driftSim").

The drift model starts with an allele frequency and follows the process:

### Initialisation


```r
Ne <- c(`1996` = 268.9, `2012` = 127.6)
pSurv <- 0.1
f0 <- seq(0.5, 0.95, by = 0.05)
nPops <- c(4, 8)
sigma <- c(0.5, 0.8)
```


1. The effective population size was defined as Ne = 268.9 representing the initial population in 1996. 
2. A starting major allele frequency was selected as one of f0 = _0.5_, _0.55_, _0.6_, _0.65_, _0.7_, _0.75_, _0.8_, _0.85_, _0.9_ and _0.95_. Each rabbit was simulated as heterozygous or homozygous for either the major or minor allele using the initial starting frquency
3. An initial survival rate was defined as p = 0.1, with this functioning as an initial bottleneck, and rabbits were assigned as survivors or fatalities with this probability.
4. These populationwas consideredas the *central population*, and this initialisation process was then repeated for either 4 or 8 neighbouring populations of the same size. However, variability was added to the initial allele frequencies on the logit scale using values of &#963; = _0.5_ and _0.8_. The same bottleneck procedure as applied to each of these populations.



```r
nSamp <- 1e06
x_0 <- runif(nSamp, 0.5, 0.95)
corSigma <- c(0.5, 0.8, 1)
y_0 <- corSigma %>%
  lapply(function(sigma){
    inv.logit(rnorm(nSamp, logit(x_0), sigma))
  })
```

In order to express these values for &#963; in terms of correlations, 1,000,000 random values were simulated for a starting frequency ($f_0$) anywhere between 0.5 and 0.95.
1,000,000 _neighbouring population starting values_ were randomly sampled around this using the values of $\sigma =$ _0.5_, _0.8_ and _1_, with the random sampling taking place on the logit scale.


| Sigma | Correlation |
|:-----:|:-----------:|
|  0.5  |    0.81     |
|  0.8  |    0.65     |
|   1   |    0.57     |

Table: Approximate correlations between starting allele frequencies in two populations, for three chosen values of $\sigma$



