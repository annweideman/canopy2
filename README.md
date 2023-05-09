# canopy2
Canopy2: tumor phylogeny inference using bulk DNA and single-cell RNA sequencing

# Author
Ann Marie Weideman

# Maintainer
Ann Marie Weideman, anndo1(at)umbc.edu

# Description
An R package that uses Bayesian and computational methods to infer the tumor phylogenetic tree from bulk DNA whole exome sequencing (WES) and single-cell RNA sequencing (scRNA-seq).

The distinguishing features of Canopy2 when compared to available methods, including predecessor Canopy, are: 1) joint inference using single nucleotide variants derived from bulk DNA WES and scRNA-seq, 2) separation of zeros categorized as non-cancerous (cells without mutations), stochastic (mutations not expressed due to bursting), and technical (expressed mutations not picked up by sequencing), and 3) a three-tier output used to comprehensively examine the tumor evolutionary history and intratumor heterogeneity. The output allows one to infer the cell-of-origin (cancer initiating cell), temporal order of the point mutations, the mutational profiles of the single-cells, and the composition of the bulk samples that comprise these single-cells. 

Inputs: 
\begin{itemize}
  \item $\boldsymbol{R^s}$: Alternative read counts from single-cell RNA sequencing data. We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts from this pipeline can be accessed in our Zenodo repository [include Zenodo link here].
  \item $\boldsymbol{X^s}$: Total read counts (benign + mutated) from single-cell RNA sequencing data.
  \item $\boldsymbol{R^b}$: Alternative read counts from bulk DNA whole exome sequencing (WES). We provide details for a pre-processing pipeline for somatic variant calling in the paper; the scripts from this pipeline can be accessed in our Zenodo repository [include Zenodo link here].
  \item $\boldsymbol{X^b}$: Total read counts (benign + mutated) from bulk DNA whole exome sequencing (WES).
  \item $\boldsymbol{\alpha}$ and $\boldsymbol{\beta}$: A $1 \times M$ (number of mutations) vector of mutation specific gene activation ($\alpha$) and gene deactivation ($\beta$) rates. These can be estimated using either the BPSC methodology in function `get_burstiness_bpsc()' (recommended) or the SCALE metholodgy with function `get_burstiness_scale()'.
  \item $\kappa$ and $tau$: $1 \times 1$ scalars used to compute the sequencing error by $\frac{\kappa}{\kappa + \tau}$. The average error rate of next-generation sequencing is reported to be 0.1\% per nucleotide, so these parameters have been set to default values of $\kappa=1$ and $\tau = 999$ such that $\frac{\kappa}{\kappa + \tau} = 0.001$ (or 0.1\%). 
 \item Klist: A range of possible numbers of subclones to use for DIC calculations. For each $k$ in Klist, the tree with highest mean posterior density across all post-burn-in iterations and all chains is selected, and then the deviance information criterion (DIC) (of type Gelman or Spiegelhalter) is used to select the optimal configuration. We suggest examining the graph of BIC versus number of subclones to determine the optimal number of subclones. 
 
Output: 
