# README for the 4DNvestigator toolbox

Scott Ronquist, scotronq@umich.edu. 1/18/19

# Overview
The 4DNvestigator is a user-friendly
MATLAB toolbox that analyzes time-series genome-wide chromosome
conformation capture data (Hi-C) and gene expression (RNA-seq).

Paper: in preparation

Availability: https://github.com/scotronq/4DNvestigator

# Hi-C and RNA-seq data types

the 4DNvestigator accepted the following input file formats:

|**Data Type**|**File Type**|**Program**|
|----|----|----|
|Hi-C|.hic|Juicer|
|RNA-seq|.genes.results|RSEM|

# Functions
- Larntz-Perlman procedure: Method for testing the equality of correlation
matrices. This is applied to Hi-C correlation matrices to determine the
significance of matrix differences.
- von Neuman entropy: Measures the entropy ("uncertainty") of a
multivariate system. Uncertainty is related to stemness. Here, we use this
measure to determine stemness of Hi-C samples.
- Chromatin partitioning: Partitioning of the genome into two distinct
groups based on either the Fiedler vector or principal component 1. This
partitioning corresponds to euchromatin and heterochromatin, or A/B
compartments.
- Differential expression: Differential expression measures the
significance of RNA-seq expression differences between samples.
- 4DN Feature Analyzer: This measures the amount of change in both genome
structure and function for specified genomic regions by mapping all time
points to a consistent low dimensional embedding, and quantifying the variance
of each loci within this space over time. Method specifics can be found in:
["Genome Architecture Mediates Transcriptional Control of Human Myogenic Reprogramming"](https://www.cell.com/iscience/fulltext/S2589-0042(18)30114-7)
- A/B switch RNA-seq: This function determines which genomic regions change
their chromatin structure from compartment "A" to compartment "B", with a
corresponding change in gene expression



