# 5xFAD-sNucSeq
Code used for the analysis of single nucleus RNA-seq (sNuc-Seq) data from WT and 5xFAD Alzheimer's disease mouse model, and an algorithmic approach to track changes in cell states across time. 

# 5xFAD-sNucSeq
Code used for the analysis of single nucleus RNA-seq (sNuc-Seq) data from WT and 5xFAD Alzheimer's disease mouse model, and an algorithmic approach to track changes in cell states across time. 

# Dynamic modeling of transitions between clusters using time course data:

FILE Time_Course_DATA.zip Contains:

RAW data used to run the MATLAB code including: (1) 2D embedding of cells from the 7 month old WT and AD mouse model 5xFAD 
(2) PCA values for all data across time points including the eignevector wights across celss 
(3) PCA values for all data across time points including the loading values across genes 
(4) Cluster assignement for each cell from the 7 month old WT and AD mouse model 5xFAD 

To run the code use:

PCgenes.all-to-AD.m and The matlab object: PCgenes.all-to-AD.mat

# Cluster analysis of single nucleus RNA-seq data in R

Functions for data initialization and clustering, variable genes, 2D embedding and visualization, using the Seurat Package and R code:

FUNCTIONS_SINGLE_CELL_ANALYSIS.r

Functions for diffusion map embedding and visualization of the results, using the Destiny Package and R code:

Diffusion_Map.r
