library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/Desktop/2025-05-20_StJudeKinnexWorkshop/scRNA_mtx/")
data <- Read10X(data.dir = "08_scKinnex.isoquant.transcript.mtx")
data_obj <- CreateSeuratObject(counts = data, project = "tutorial", min.cells = 3, min.features = 200)
data_obj
data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^MT-")
VlnPlot(data_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data_obj <- NormalizeData(data_obj)

data_obj <- FindVariableFeatures(data_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(data_obj)
data_obj <- ScaleData(data_obj, features = all.genes)

data_obj <- RunPCA(data_obj, features = VariableFeatures(object = data_obj))

# Examine and visualize PCA results a few different ways
print(data_obj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(data_obj, dims = 1:2, reduction = "pca")

DimPlot(data_obj, reduction = "pca")
DimHeatmap(data_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data_obj, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(data_obj)


data_obj <- FindNeighbors(data_obj, dims = 1:10)
data_obj <- FindClusters(data_obj, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(data_obj), 5)


data_obj <- RunUMAP(data_obj, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(data_obj, reduction = "umap", label = TRUE)
DimPlot(data_obj, reduction = "umap", label = TRUE, group.by = "orig.ident")

# find all markers of cluster 4
cluster4.markers <- FindMarkers(data_obj, ident.1 = 4)
head(cluster4.markers, n = 5)
VlnPlot(data_obj, features = c("ENST00000368738.4"))
FeaturePlot(data_obj, features = c("ENST00000368738.4"))
