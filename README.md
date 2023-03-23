# UNC93B1-variants-cause-TLR7-dependent-autoimmunity

This repository contains jupyter notebooks containing the codes used for the analysis and to replicate the results of the single-cell data analysis from the manuscript:


C. Wolf et al., “UNC93B1 variants cause TLR7-dependent autoimmunity “(Submitted to Science, 2023).

The purpose of this repository is to provide additional details beyond the method section of the manuscript for the interested reader.

The codes contain the codes for the following analysis:


•	INF_preprocessing: Input for this script is the Cell Ranger “raw_feature_bc_matrix.h5” files.


•	INF_umap_and_clustering: The output of the previous script is used to perform dimensionality reduction using UMAP and clustering to assign each cell to a canonical immune cell type. After dissection of leukocytes into broad populations (such as CD4+ T cells, CD8+ T cells, myeloid cells, etc.) each of these subsets were further dissected in an iterative manner using the same computational pipeline.


•	INF_ploting

•	INF_GSEA_and_plots
