#' Pre-processed breast cancer (BC) data from patient BC03
#'
#' Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
#' (bulk RNA-seq), and bulk DNA whole exome sequencing data (bulk DNA WES) from
#' one patient (BC03) of 11 individuals (BC01-BC11) with primary breast cancer
#' \insertCite{Chung2017}{canopy2}. The single-nucleotide variant (SNV) pipeline
#' utilized to pre-process the data is described in the main text and supplement.
#' In the main study, the authors sequenced 549 primary breast cancer cells and
#' matched bulk tumors and/or pooled cells from 11 patients, including two lymph node
#' metastases (BC03LN, BC07LN). Their final dataset after filtering included 515
#' single cells and 14 bulk samples. <br />
#' Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
# * Blood: SRR3023076
# * Bulk DNA WES (tumor): SRR3023078, SRR3023080
# * Bulk RNA (tumor): SRR2973275
# * Bulk RNA (LN tumor): SRR2973276
# * Single-cell RNA (primary): SRR2973351-SRR2973383 and SRR5023442-SRR5023445
# * Single-cell RNA (LN): SRR2793384 - SRR2973436 and SRR5023446, SRR5023447
#'
#' \code{BC03_preproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{BC03_preproc@alt.qc} &nbsp; A 145 x 80 matrix of pre-processed alternative
#' read counts from the bulk and single-cell data for BC03. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
#' above labels from the SNV callset of patient BC03.
#'
#' \code{BC03_preproc@annovar} &nbsp; A 145 x 39 dataframe containing the
#' functional annotation of the somatic variants for BC03. For details regarding
#' the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{canopy2}.
#'
#' \code{BC03_preproc@featurecounts.qc} &nbsp; A 63677 x 75 matrix of
#' pre-processed mapped reads from single-cell gene expression data. The row
#' names denote the gene names. The column names from the third column onward
#' denote the scRNA-seq samples from the gene expression data for patient BC03.
#'
#' \code{BC03_preproc@mut} &nbsp; A 145 x 5 dataframe containing the chromosomal
#' positions of point mutations for BC03. The columns consist of: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{BC03_preproc@ref.qc} &nbsp; A 145 x 80 matrix of pre-processed
#' reference read counts from the bulk and single-cell data for BC03. The
#' row names denote the chromosomal position of the point mutations. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient BC03.
#'
#' \code{BC03_preproc@summary.qc} &nbsp; A 14 x 75 matrix of the
#' \code{featureCounts} output summary that reports assignment of alignments
#' to genomic features. The row names denote the assignment of the alignments.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient BC03.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (pre-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC03  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 55 (45)          	|
#'  |       | Primary tumor               | 1              	| 1 (pooled bulk RNA) | 37 (30)          	|
#'  |       | (total)                   	| 3              	| 2              	| 92 (75)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC03_preproc"

#' Post-processed breast cancer (BC) data from patient BC03
#'
#' Some further data munging of BC03_preproc to get BC03_postproc, which contains
#' the finalized read counts for Canopy2 functions and vignette.  <br />
#' Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
# * Blood: SRR3023076
# * Bulk DNA WES (tumor): SRR3023078, SRR3023080
# * Bulk RNA (tumor): SRR2973275
# * Bulk RNA (LN tumor): SRR2973276
# * Single-cell RNA (primary): SRR2973351-SRR2973383 and SRR5023442-SRR5023445
# * Single-cell RNA (LN): SRR2793384 - SRR2973436 and SRR5023446, SRR5023447
#'
#' \code{BC03_postproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{BC03_postproc@annovar.qc} &nbsp; A 33 x 39 dataframe containing the
#' post-processed functional annotations of the somatic variants. For details
#' regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{canopy2}. Only those
#' annotations corresponding to mutations with at least one non-zero count
#' across all cells were retained.
#'
#' \code{BC03_postproc@featurecounts.qc} &nbsp; A 63677 x 77 matrix of
#' post-processed mapped reads from the single-cell gene expression data.
#' The row names denote the gene names. The column names from the third column
#' onward denote the scRNA-seq samples from the gene expression data for patient
#' BC03. Only those cells that have at least one non-zero count across genes
#' were retained.
#'
#' \code{BC03_postproc@param.est} &nbsp; A list of 5 elements: \code{alpha},
#' \code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
#' three are of length equal to the number of mutations where \code{alpha} denotes the
#' mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
#' gene deactivation rates, and \code{scale} denotes the mutation-specific
#' transcription rates. The ids in \code{id.g} denote the positions corresponding
#' to the mutations that were estimable (returning non-NA or non-negative values),
#' and \code{pct.estimable} denotes the percentage that were estimable.
#'
#' \code{BC03_postproc@Rb} &nbsp; A 31 x 2 matrix of post-processed alternative
#' read counts from the bulk DNA WES data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the bulk samples
#' according to the above labels from the SNV callset of patient BC03.
#' Post-processing involved removing mutations (rows) with zero alternative
#' read counts across all cells from the single-cell data.
#'
#' \code{BC03_postproc@Rs} &nbsp; A 31 x 75 matrix of post-processed alternative
#' read counts from the scRNA-seq data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the scRNA-seq samples
#' according to the above labels from the SNV callset of patient BC03.
#' Post-processing involved i) removing mutations (rows) with zero alternative
#' read counts across all cells and ii) removing cells (columns) with zero total
#' (benign + mutated) read counts across all mutations.
#'
#' \code{BC03_postproc@Xb} &nbsp; A 31 x 2 matrix of post-processed total read
#' counts (benign + mutated) from the bulk DNA WES data. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk samples according to the above labels from the SNV callset of patient
#' BC03. Post-processing involved removing mutations (rows) with zero
#' alternative read counts across all cells from the single-cell data.
#'
#' \code{BC03_postproc@Xs} &nbsp; A 31 x 75 matrix of post-processed total read
#' counts (benign + mutated) from the scRNA-seq data. The row names denote the
#' chromosomal position of the point mutations. The column names denote the
#' scRNA-seq samples according to the above labels from the SNV callset of
#' patient BC03. Post-processing involved i) removing mutations (rows) with zero
#' alternative read counts across all cells and ii) removing cells (columns)
#' with zero total (benign + mutated) read counts across all mutations.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (post-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC03  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 55 (45)          	|
#'  |       | Primary tumor               | 1              	| 1 (pooled bulk RNA) | 37 (30)          	|
#'  |       | (total)                   	| 3              	| 2              	| 92 (75)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC03_postproc"
#'
#' Pre-processed breast cancer (BC) Data from Patient BC07
#'
#' Pre-processed single-cell RNA sequencing (scRNA-seq), bulk RNA sequencing
#' (bulk RNA-seq), and bulk DNA whole exome sequencing data (bulk DNA WES)
#' from one patient (BC07) of 11 individuals (BC01-BC11) with primary breast
#' cancer \insertCite{Chung2017}{canopy2}. The single-nucleotide variant (SNV)
#' pipeline utilized to pre-process the data is described in the main text and
#' supplement. In the main study, the authors sequenced 549 primary breast
#' cancer cells and matched bulk tumors and/or pooled cells from 11 patients,
#' including two lymph node metastases (BC03LN, BC07LN). Their final dataset
#' after filtering included 515 single cells and 14 bulk samples. <br />
#' GEO Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
# * Blood: SRR3023081
# * Bulk DNA WES (tumor): SRR3023082, SRR3023083
# * Bulk RNA (tumor): SRR2973277
# * Bulk RNA (LN tumor): SRR2973278
# * Single-cell RNA (primary): SRR2973437-SRR2973484, SRR5023558-SRR5023560
# * Single-cell RNA (LN): SRR2973485-SRR2973535, SRR5023561-SRR5023562
#'
#' \code{BC07_preproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{BC07_preproc@alt.qc} &nbsp; A 79 x 106 matrix of pre-processed alternative
#' read counts from the bulk and single-cell data for BC07. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood samples according to the
#' above labels from the SNV callset of patient BC07.
#'
#' \code{BC07_preproc@annovar} &nbsp; A 79 x 39 dataframe containing the
#' functional annotation of the somatic variants for BC07. For details regarding
#' the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{canopy2}.
#'
#' \code{BC07_preproc@featurecounts.qc} &nbsp; A 63677 x 101 matrix of
#' pre-processed mapped reads from single-cell gene expression data. The row
#' names denote the gene names. The column names from the third column onward
#' denote the scRNA-seq samples from the gene expression data for patient BC07.
#'
#' \code{BC07_preproc@mut} &nbsp; A 79 x 5 dataframe containing the chromosomal
#' positions of point mutations for BC07. The columns consist of: the chromosome
#' number \code{Chr}, the starting position of the point mutation \code{Start},
#' the ending position of the point mutation \code{End} (coincides with
#' \code{Start} because the mutation is a single nucleotide variant), the
#' reference allele \code{Ref}, and the alternate allele \code{Alt}.
#'
#' \code{BC07_preproc@ref.qc} &nbsp; A 79 x 106 matrix of pre-processed
#' reference read counts from the bulk and single-cell data for BC07. The
#' row names denote the chromosomal position of the point mutations. The
#' column names denote the bulk RNA-seq, scRNA-seq, bulk DNA WES, and blood
#' samples from the total read counts of patient BC07.
#'
#' \code{BC07_preproc@summary.qc} &nbsp; A 14 x 101 matrix of the
#' \code{featureCounts} output summary that reports assignment of alignments
#' to genomic features. The row names denote the assignment of the alignments.
#' The column names denote the scRNA-seq samples from the gene expression data
#' for patient BC07.
#'
#' @details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (pre-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC07  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 53 (52)          	|
#'  |       | Primary tumor               | 1              	| 1 (tumor bulk RNA) | 51 (49)          	|
#'  |       | (total)                   	| 3              	| 2              	| 104 (101)        	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC07_preproc"
#'
#' Post-processed breast cancer (BC) data from patient BC07
#'
#' Some further data munging of BC07_preproc to get BC07_postproc, which contains
#' the finalized read counts for Canopy2 functions and vignette. <br />
#' GEO Accession Number: GSE75688
#'
#' @format
#'
#' Labels:
# * Blood: SRR3023081
# * Bulk DNA WES (tumor): SRR3023082, SRR3023083
# * Bulk RNA (tumor): SRR2973277
# * Bulk RNA (LN tumor): SRR2973278
# * Single-cell RNA (primary): SRR2973437-SRR2973484, SRR5023558-SRR5023560
# * Single-cell RNA (LN): SRR2973485-SRR2973535, SRR5023561-SRR5023562
#'
#' \code{BC07_postproc} &nbsp; An object of class \code{S4}. Access internal
#' objects from slots using @.
#'
#' \code{BC07_postproc@annovar.qc} &nbsp; A 31 x 39 dataframe containing the
#' post-processed functional annotations of the somatic variants. For details
#' regarding the column headers, we refer readers to the extensive
#' ANNOVAR [documentation](https://annovar.openbioinformatics.org/) and the
#' original research article \insertCite{Wang2010}{canopy2}. Only those
#' annotations corresponding to mutations with at least one non-zero count
#' across all cells were retained.
#'
#' \code{BC07_postproc@featurecounts.qc} &nbsp; A 63677 x 102 matrix of
#' post-processed mapped reads from the single-cell gene expression data.
#' The row names denote the gene names. The column names from the third column
#' onward denote the scRNA-seq samples from the gene expression data for patient
#' BC07. Only those cells that have at least one non-zero count across genes
#' were retained.
#'
#' \code{BC07_postproc@param.est} &nbsp; A list of 5 elements: \code{alpha},
#' \code{beta}, \code{scale}, \code{id.g}, and \code{pct.estimable}. The first
#' three are of length equal to the number of mutations where \code{alpha} denotes the
#' mutation-specific gene activation rates, \code{beta} denotes the mutation-specific
#' gene deactivation rates, and \code{scale} denotes the mutation-specific
#' transcription rates. The ids in \code{id.g} denote the positions corresponding
#' to the mutations that were estimable (returning non-NA and non-negative values),
#' and \code{pct.estimable} denotes the percentage that were estimable.
#'
#' \code{BC07_postproc@Rb} &nbsp; A 30 x 2 matrix of post-processed alternative
#' read counts from the bulk DNA WES data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the bulk samples
#' according to the above labels from the SNV callset of patient BC07.
#' Post-processing involved removing mutations (rows) with zero alternative
#' read counts across all cells from the single-cell data.
#'
#' \code{BC07_postproc@Rs} &nbsp; A 30 x 100 matrix of post-processed alternative
#' read counts from the scRNA-seq data. The row names denote the chromosomal
#' position of the point mutations. The column names denote the scRNA-seq samples
#' according to the above labels from the SNV callset of patient BC07.
#' Post-processing involved i) removing mutations (rows) with zero alternative
#' read counts across all cells and ii) removing cells (columns) with zero total
#' (benign + mutated) read counts across all mutations.
#'
#' \code{BC07_postproc@Xb} &nbsp; A 30 x 2 matrix of post-processed total read
#' counts (benign + mutated) from the bulk DNA WES data. The row names denote
#' the chromosomal position of the point mutations. The column names denote the
#' bulk samples according to the above labels from the SNV callset of patient
#' BC07. Post-processing involved removing mutations (rows) with zero
#' alternative read counts across all cells from the single-cell data.
#'
#' \code{BC07_postproc@Xs} &nbsp; A 30 x 100 matrix of post-processed total read
#' counts (benign + mutated) from the scRNA-seq data. The row names denote the
#' chromosomal position of the point mutations. The column names denote the
#' scRNA-seq samples according to the above labels from the SNV callset of
#' patient BC07. Post-processing involved i) removing mutations (rows) with zero
#' alternative read counts across all cells and ii) removing cells (columns)
#' with zero total (benign + mutated) read counts across all mutations.
#'
#'@details
#'
#'  | ***Patient*** | ***Subtype*** | ***# bulk DNA WES*** | ***# bulk RNA-seq*** | ***# scRNA-seq (post-processed)***|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  |-----	|-----	|-----	|-----	|-----	|
#'  | BC07  | Blood                     	| 1              	| 0              	| 0                	|
#'  |       | Metastatic lymph node       | 1              	| 1 (pooled bulk RNA) | 53 (52)          	|
#'  |       | Primary tumor               | 1              	| 1 (tumor bulk RNA) | 51 (48)          	|
#'  |       | (total)                   	| 3              	| 2              	| 104 (100)         	|
#'
#' @md
#'
#' @import Rdpack
#' @references{
#'    \insertAllCited{}
#' }
#'
"BC07_postproc"
