# README for the GSFAT toolbox

Scott Ronquist, scotronq@umich.edu. 12/22/18

# Overview
The Genome Structure Function Analysis Toolbox (GSFAT) is a user-friendly
MATLAB toolbox that analyzes time-series genome-wide chromosome
conformation capture data (Hi-C) and gene expression (RNA-seq).

Paper: TBD
Availability: https://github.com/scotronq/gsfat

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
- 4DN Feature Analyzer: 
- A/B switch RNA-seq: 

# In progress



