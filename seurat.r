
###Identify directory with single cell matrix or desired output folder
working_dir = directory = 'C:/Users/desha/Documents/Grad/Research/Singlecell/'
 setwd(working_dir)
###Essential software
library(dplyr)
library(Seurat)
library(patchwork)

# Load the Single cell feature matrix
pbmc.data <- Read10X(data.dir = "C:/Users/danie/OneDrive/Documents/Grad school/Research/singlecell/Singlecell_matrix/filtered_feature_bc_matrixT2")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "T2", min.cells = 3, min.features = 200)
pbmc

####Identify the percent of reads mapped to mitochondrial gene
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "MT")


###visualize UMI counts and Gene feature counts per cell

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


###Create a subset of cells filtered for number of UMI counts to eliminate doublets and dead cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)

###Normalization methods
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

###Identify variable features after normalization
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10870)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(pbmc)

##Scale data for downstream clustering
pbmc <- ScaleData(pbmc, features = all.genes)

###PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


DimPlot(pbmc, reduction = "pca")


###Create heatmaps for differential expressed genes
DimHeatmap(pbmc, dims = 1, cells = 5000, balanced = TRUE)


DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

####Determine number of dimensions to be used for UMAP and TSNE

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)


ElbowPlot(pbmc)

###Clustering and neighbor ID based on number of dims identified above.
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.4)

head(Idents(pbmc), 5)

###UMAP
pbmc <- RunUMAP(pbmc, dims = 1:15)


DimPlot(pbmc, reduction = "umap")


saveRDS(pbmc, file = "Merged_T1T2")
readRDS(pbmc, file = "Testis2")



###Identify markers in each cluster

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


###Visuzlize expression of marker features
VlnPlot(pbmc, features = c("ENSGACG00000005864"))


VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(pbmc1, features = c("ENSGACG00000005864", "ENSGACG00000003679", "ENSGACG00000005335", "ENSGACG00000014492")
	
WhichCells(pbmc, idents = "0")	

0.raw.data <- as.matrix(GetAssayData(pbmc, slot = "counts")[, WhichCells(pbmc, ident = "0")])


subset(pbmc, idents = c("3", "5", "4"), invert = TRUE)
