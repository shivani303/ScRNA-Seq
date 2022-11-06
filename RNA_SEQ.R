library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
data_dir <- 'C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/rna_seq/'
list.files(data_dir)
expression_matrix <- Read10X(data.dir = data_dir)
leukemia = CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200) 
leukemia
expression_matrix[1:50, 1:10]

leukemia[["percent.mt"]] = PercentageFeatureSet(leukemia, pattern = "^MT-")
tail(leukemia@meta.data)
VlnPlot(leukemia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot = FeatureScatter(leukemia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot

leukemia= subset(leukemia, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
leukemia

leukemia = NormalizeData(leukemia)
leukemia = FindVariableFeatures(leukemia, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(leukemia), 10)
top10

plot1 = VariableFeaturePlot(leukemia)
plot1
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes = rownames(leukemia)
leukemia = ScaleData(leukemia, features = all.genes)

leukemia@assays$RNA@scale.data[1:50, 1:5]
leukemia = RunPCA(leukemia, features = VariableFeatures(object = leukemia))
DimHeatmap(leukemia, dims = 1:25, cells = 500, balanced = TRUE)

ElbowPlot(leukemia)
leukemia = FindNeighbors(leukemia, dims = 1:20)
leukemia = FindClusters(leukemia, resolution = 0.5)
head(leukemia@meta.data)
leukemia = RunUMAP(leukemia, dims = 1:20)
DimPlot(leukemia, reduction = "umap")
DimPlot(leukemia, reduction = "umap", label = T)
leukemia.markers = FindAllMarkers(leukemia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(leukemia.markers)

a = leukemia.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a
genes = a %>% pull(gene)
genes
library(celldex)
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(DietSeurat(leukemia))
sce
library(SingleR)

hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)

table(hpca.main$pruned.labels)
table(hpca.fine$pruned.labels)

leukemia@meta.data$hpca.main   <- hpca.main$pruned.labels
leukemia@meta.data$hpca.fine   <- hpca.fine$pruned.labels

leukemia <- SetIdent(leukemia, value = "hpca.fine")
DimPlot(leukemia, label = T , repel = T, label.size = 3) + NoLegend()
DimPlot(leukemia, reduction = "umap", label = T,repel = T,pt.size = 0.5) +NoLegend()

FeaturePlot(leukemia, features = genes[1:5])

FeaturePlot(leukemia, features = genes[15:20])












