https://github.com/annweideman/canopy2/assets/52984250/6bd2366d-ceff-4d39-a0a1-aefb8611b8e8

# Important Note (added on 12/3/23): 
If you receive the following error while using the `get_gene_expression()` function, which internally uses `biomaRt::useEnsembl()`:
```
Error in `collect()`:
! Failed to collect lazy table.
Caused by error in `db_collect()`:
! Arguments in `...` must be used.
Problematic argument:
 ..1 = Inf
 Did you misspell an argument name?
Run `rlang::last_trace()` to see where the error occurred.
```

This is due to a bug relating to BiocFileCache compatibility with the new version dbplyr. To correct, run the following:

```
library(BiocManager)
install("BiocFileCache")
```

[https://support.bioconductor.org/p/9154901/](https://stat.ethz.ch/pipermail/bioc-devel/2023-October/020003.html)

# Method
Canopy2: tumor phylogeny inference using bulk DNA and single-cell RNA sequencing

# Author
Ann Marie Weideman

# Maintainer
Ann Marie Weideman, anndo1(at)umbc.edu

# Description
An R package that uses Bayesian and computational methods to infer the tumor phylogenetic tree from bulk DNA whole exome sequencing (WES) and single-cell RNA sequencing (scRNA-seq).

The distinguishing features of Canopy2 when compared to available methods, including predecessor Canopy, are: 1) joint inference using single nucleotide variants derived from bulk DNA WES and scRNA-seq, 2) separation of zeros categorized as non-cancerous (cells without mutations), stochastic (mutations not expressed due to bursting), and technical (expressed mutations not picked up by sequencing), and 3) a three-tier output used to comprehensively examine the tumor evolutionary history and intratumor heterogeneity. The output allows one to infer the cell-of-origin (cancer initiating cell), temporal order of the point mutations, the mutational profiles of the single-cells, and the composition of the bulk samples that comprise these single-cells. 

**Inputs:**
  * `Rs`:   $M \times N$ (mutation by single-cell) matrix of alternative read counts from single-cell RNA sequencing data. We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts for this pipeline can be accessed in our Zenodo repository (https://zenodo.org/record/7931384).
  * `Xs`:  $M \times N$ (mutation by single-cell) matrix of total read counts (benign + mutated) from single-cell RNA sequencing data.
  *  `Rb`:  $M \times T$ (mutation by bulk sample) matrix Alternative read counts from bulk DNA whole exome sequencing (WES). We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts for this pipeline can be accessed in our Zenodo repository (https://zenodo.org/record/7931384).
  *  `Xb`: Total read counts (benign + mutated) from bulk DNA whole exome sequencing (WES).
  *  `alpha` and `beta`:  $1 \times M$ (number of mutations) vectors of mutation specific gene activation ($\alpha$) and gene deactivation ($\beta$) rates. These can be estimated using either the BPSC methodology in function `get_burstiness_bpsc()` (recommended) or the SCALE methodology with function `get_burstiness_scale()`.
  *  `kappa` and `tau`:   $1 \times 1$ scalars used to compute the sequencing error by $\frac{\kappa}{\kappa + \tau}$. The average error rate of next-generation sequencing is reported to be 0.1\% per nucleotide, so these parameters have been set to default values of $\kappa=1$ and $\tau = 999$ such that $\frac{\kappa}{\kappa + \tau} = 0.001$ (or 0.1\%). 
  *  `Klist`:   A range of possible numbers of subclones to use for DIC calculations. For each $k$ in Klist, the tree with highest mean posterior density across all post-burn-in iterations and all chains is selected, and then the deviance information criterion (DIC) (of type Gelman or Spiegelhalter) is used to select the optimal configuration. We suggest examining the graph of BIC versus number of subclones to determine the optimal number of subclones. 
 
**Outputs:** Three-tier configuration consisting of
* `Z`: $M \times K$ clonal tree configuration matrix, a binary matrix which assigns mutations to clones
* `Ps`: $K \times N$ cell-to-clone assignment matrix, a binary matrix which assigns single-cells to clones
* `Pb`: $K \times T$ sample-to-clone assignment matrix, which indicates what percentage of each bulk sample consists of each clone

# Installation
Install the current release from CRAN (not recommended). Not published on CRAN as of 6/24/23.

```
install.packages("canopy2")
```

Install the devel version from GitHub (HIGHLY recommended, as this will allow you to install any corrected bugs post-publication to CRAN)

```
# Dependencies
list.of.packages<-c("ape", "bayesplot", "devtools", "BiocManager", "coda", "DirichletReg", "ggh4x", "ggplot2", "ggplotify", "ggpubr", "ggthemes",
"ggtree", "gridExtra", "gtools", "Rdpack", "SCALE", "viridis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

BiocManager::install("biomaRt")
devtools::install_github("nghiavtr/BPSC")
devtools::install_github("annweideman/canopy2")

```
# Vignette
[Canopy2 Vignette]https://htmlpreview.github.io/?https://github.com/annweideman/canopy2/blob/main/vignettes/canopy2.html




