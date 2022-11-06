# ScRNA-Seq
TITLE- Gene expression profile at single-cell level of hematopoietic stem and progenitor cells (HSPCs) from individuals with chronic myelomonocytic leukemia (CMML)

GEO ACC NO- GSE211033

PLATFORM -Illumina NovaSeq 6000

# calling library

library(Matrix)

library(Seurat)

library(dplyr)

library("patchwork")

# read a file
data_dir <- 'C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/rna_seq/'

list.files(data_dir)

expression_matrix <- Read10X(data.dir = data_dir)

#create seurat object

leukemia = CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200) 

leukemia

expression_matrix[1:50, 1:10]

**#remove cell which have high expression of mitochondrial genes because either it is in replicating phase or in apoptosis phase**

leukemia[["percent.mt"]] = PercentageFeatureSet(leukemia, pattern = "^MT-")

tail(leukemia@meta.data)

VlnPlot(leukemia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

![image](https://user-images.githubusercontent.com/66779651/200185280-5481970d-6de6-43cc-98ee-de89aab10248.png)

# FeatureScatter
plot = FeatureScatter(leukemia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot

![image](https://user-images.githubusercontent.com/66779651/200185294-06958980-1f4e-4399-ae32-845d12bba5bf.png)


leukemia= subset(leukemia, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

leukemia

# NORMALIZE DATA

leukemia = NormalizeData(leukemia)

leukemia = FindVariableFeatures(leukemia, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

top10 = head(VariableFeatures(leukemia), 10)

top10

![image](https://user-images.githubusercontent.com/66779651/200185338-860166a2-543e-46b4-979f-b9cdf5c68fb5.png)

**plot variable features with and without labels**

plot1 = VariableFeaturePlot(leukemia)

plot1

![image](https://user-images.githubusercontent.com/66779651/200185344-e76e663d-bbeb-436a-8020-bb739a72ef4c.png)

plot2 = LabelPoints(plot = plot2, points = top10, repel = TRUE)

plot2
![image](https://user-images.githubusercontent.com/66779651/200185397-0c183cae-15c3-4f9a-8e27-2973d6236363.png)


all.genes = rownames(leukemia)

leukemia = ScaleData(leukemia, features = all.genes)

leukemia@assays$RNA@scale.data[1:50, 1:5]

![image](https://user-images.githubusercontent.com/66779651/200185454-acd69135-da2c-424d-b0bf-cb9a1bb1ddd9.png)

# Principal Component Analysis

leukemia = RunPCA(leukemia, features = VariableFeatures(object = leukemia))

DimHeatmap(leukemia, dims = 1:25, cells = 500, balanced = TRUE)

![image](https://user-images.githubusercontent.com/66779651/200185482-cbcd9d8c-efdc-46a9-9012-7f75ec03e13e.png)

**ranking of principle components based on the percentage of variance**

ElbowPlot(leukemia)

![image](https://user-images.githubusercontent.com/66779651/200185495-d9641edc-8cd8-4d80-8aca-0baf4a30c12c.png)

# CLUSTERING 

leukemia = FindNeighbors(leukemia, dims = 1:20)

leukemia = FindClusters(leukemia, resolution = 0.5)

head(leukemia@meta.data)

leukemia = RunUMAP(leukemia, dims = 1:20)

DimPlot(leukemia, reduction = "umap")

![image](https://user-images.githubusercontent.com/66779651/200185577-98644d44-ca4d-443e-ad70-1172ab5ae00c.png)

DimPlot(leukemia, reduction = "umap", label = T)
![image](https://user-images.githubusercontent.com/66779651/200185585-1821b230-fe2e-43ba-b37f-6b6ba6f6dffe.png)

leukemia.markers = FindAllMarkers(leukemia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(leukemia.markers)

a = leukemia.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a
genes = a %>% pull(gene)
genes

# CELL TYPE IDENTIFICATION
library(celldex)

**get reference datasets from celldex package**

hpca.ref <- celldex::HumanPrimaryCellAtlasData()

library(SingleCellExperiment)

**convert our Seurat object to single cell experiment (SCE) for convenience**

sce <- as.SingleCellExperiment(DietSeurat(leukemia))

sce

library(SingleR)

hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)

hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)

**summary of general cell type annotations** 

table(hpca.main$pruned.labels)

![image](https://user-images.githubusercontent.com/66779651/200187323-8dd2ba71-105c-4bb4-9cae-4a8b439160b2.png)

table(hpca.fine$pruned.labels)

![image](https://user-images.githubusercontent.com/66779651/200187342-87daaf76-5a64-442d-a506-a2588ff3e6f3.png)

**add the annotations to the Seurat object metadata so we can use them**

leukemia@meta.data$hpca.main   <- hpca.main$pruned.labels

leukemia@meta.data$hpca.fine   <- hpca.fine$pruned.labels

leukemia <- SetIdent(leukemia, value = "hpca.fine")

DimPlot(leukemia, label = T , repel = T, label.size = 3) + NoLegend()

DimPlot(leukemia, reduction = "umap", label = T,repel = T,pt.size = 0.5) +NoLegend()

![image](https://user-images.githubusercontent.com/66779651/200187085-8a1e1ad9-97b6-4bef-8d08-9345a34fbd43.png)


