# canopy2
Canopy2: tumor phylogeny inference using bulk DNA and single-cell RNA sequencing

# Author
Ann Marie Weideman

# Maintainer
Ann Marie Weideman, anndo1(at)umbc.edu

# Description
An R package that uses Bayesian and computational methods to infer the tumor phylogenetic tree from bulk DNA whole exome sequencing (WES) and single-cell RNA sequencing (scRNA-seq).

The distinguishing features of Canopy2 when compared to available methods, including predecessor Canopy, are: 1) joint inference using single nucleotide variants derived from bulk DNA WES and scRNA-seq, 2) separation of zeros categorized as non-cancerous (cells without mutations), stochastic (mutations not expressed due to bursting), and technical (expressed mutations not picked up by sequencing), and 3) a three-tier output used to comprehensively examine the tumor evolutionary history and intratumor heterogeneity. The output allows one to infer the cell-of-origin (cancer initiating cell), temporal order of the point mutations, the mutational profiles of the single-cells, and the composition of the bulk samples that comprise these single-cells. 

**Inputs:**
  * `Rs`:   $M \times N$ (mutation by single-cell) matrix of alternative read counts from single-cell RNA sequencing data. We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts for this pipeline can be accessed in our Zenodo repository [include Zenodo link here].
  * `Xs`:  $M \times N$ (mutation by single-cell) matrix of total read counts (benign + mutated) from single-cell RNA sequencing data.
  *  `Rb`:  $M \times T$ (mutation by bulk sample) matrix Alternative read counts from bulk DNA whole exome sequencing (WES). We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts for this pipeline can be accessed in our Zenodo repository [include Zenodo link here].
  *  `Xb`: Total read counts (benign + mutated) from bulk DNA whole exome sequencing (WES).
  *  `alpha` and `beta`:  $1 \times M$ (number of mutations) vectors of mutation specific gene activation ($\alpha$) and gene deactivation ($\beta$) rates. These can be estimated using either the BPSC methodology in function `get_burstiness_bpsc()` (recommended) or the SCALE methodology with function `get_burstiness_scale()`.
  *  `kappa` and `tau`:   $1 \times 1$ scalars used to compute the sequencing error by $\frac{\kappa}{\kappa + \tau}$. The average error rate of next-generation sequencing is reported to be 0.1\% per nucleotide, so these parameters have been set to default values of $\kappa=1$ and $\tau = 999$ such that $\frac{\kappa}{\kappa + \tau} = 0.001$ (or 0.1\%). 
  *  `Klist`:   A range of possible numbers of subclones to use for DIC calculations. For each $k$ in Klist, the tree with highest mean posterior density across all post-burn-in iterations and all chains is selected, and then the deviance information criterion (DIC) (of type Gelman or Spiegelhalter) is used to select the optimal configuration. We suggest examining the graph of BIC versus number of subclones to determine the optimal number of subclones. 
 
**Outputs:** Three-tier configuration consisting of
* `Z`: $M \times K$ clonal tree configuration matrix, a binary matrix which assigns mutations to clones
* `Ps`: $K \times N$ cell-to-clone assignment matrix, a binary matrix which assigns single-cells to clones
* `Pb`: $K \times T$ sample-to-clone assignment matrix, which indicates what percentage of each bulk sample consists of each clone

#Installation
Install the current release from CRAN (not recommended). NOT PUBLISHED ON CRAN YET SO WON'T WORK (5/9/23)

```
install.packages("canopy2")
```

Install the devel version from GitHub (HIGHLY recommended, as this will allow you to install any corrected bugs post-publication to CRAN)

```
install.packages(c("ape", "bayesplot", "biomaRt", "BPSC", "coda", "DirichletReg", "ggh4x", "ggplot2", "ggplotify", "ggpubr", "ggthemes",
"ggtree", "gridExtra", "gtools", "Rdpack", "SCALE", "viridis"))
devtools::install_github("annweideman/canopy2/package")
```




