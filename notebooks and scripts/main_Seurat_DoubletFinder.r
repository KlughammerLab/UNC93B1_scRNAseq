#stting working directory
setwd("/home/mmokhtari/INF_analysis")

#loading libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(DoubletFinder)

#loading the donor_ids
donor_id <- read.table("/home/mmokhtari/INF_analysis/donor_ids.tsv",header = TRUE)
row.names(donor_id) <-donor_id$cell
head(donor_id)
genetic_doublets <- donor_id[which(donor_id$donor_id== "doublet"),1]

#loading count matrices
data_raw <- Read10X_h5("/home/mmokhtari//INF_analysis/raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)

#creating the Seurat object
data <- CreateSeuratObject(counts = data_raw)

Idents(data)<- row.names(data@meta.data)

#QC and filtering
data$mitoPercent <-PercentageFeatureSet(data,pattern = "^MT")
data_filtered <- subset(data, subset = nCount_RNA > 800 &nFeature_RNA > 500 &mitoPercent < 10)
data_filtered@meta.data$donor_id <- donor_id[row.names(data_filtered@meta.data),2]
data_filtered<- subset(x = data_filtered, subset = donor_id == c("doublet","unassigned"), invert = TRUE)

Idents(data_filtered) <- data_filtered@meta.data$orig.ident


# pre-process standard workflow
data_filtered <- NormalizeData(object = data_filtered)
data_filtered <- FindVariableFeatures(object = data_filtered)
data_filtered <- ScaleData(object = data_filtered)
data_filtered <- RunPCA(object = data_filtered)
ElbowPlot(data_filtered)


data_filtered <- FindNeighbors(object = data_filtered, dims = 1:13)
data_filtered <- FindClusters(object = data_filtered)
data_filtered <- RunUMAP(object = data_filtered, dims = 1:20)


## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_data <- paramSweep_v3(data_filtered, PCs = 1:13, sct = FALSE)
sweep.stats_data <- summarizeSweep(sweep.res.list_data, GT = FALSE)
bcmvn_data <- find.pK(sweep.stats_data)


ggplot(bcmvn_data, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()


pK <- bcmvn_data %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- data_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.09*nrow(data_filtered@meta.data))  ## Assuming 9% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
data_filtered <- doubletFinder_v3(data_filtered, 
                                     PCs = 1:13, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(data_filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_1405")

write.table(data_filtered@meta.data,"doublets_1.tsv", sep="\t")