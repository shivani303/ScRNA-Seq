library(Matrix)
library(Seurat)
library(dplyr)
library("patchwork")
ss_file= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/GSE185498_raw_counts_aggregated.csv")

#check duplicates
length(unique(ss_file$Gene.names)) == nrow(ss_file)

#remove duplicate value 
ss_file=ss_file[!duplicated(ss_file$Gene.names),]

file <- ss_file[,-1]
rownames(file) <- ss_file[,1]
View(file)

cell1= CreateSeuratObject(counts = file, min.cells = 3, min.features = 200) 

cell1

cell1[1:50, 1:10]

cell1[["percent.mt"]] = PercentageFeatureSet(cell1, pattern = "^MT-")
head(cell1@meta.data)
VlnPlot(cell1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

cell = subset(cell1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell

#  removing unwanted cells from the dataset, the next step is to normalize the data

cell = NormalizeData(cell)
cell = FindVariableFeatures(cell, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(cell), 10)
top10

plot1 = VariableFeaturePlot(cell)
plot1
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes = rownames(cell)
cell = ScaleData(cell, features = all.genes)

## Centering and scaling data matrix

cell@assays$RNA@scale.data[1:50, 1:5]

cell = RunPCA(cell, features = VariableFeatures(object = cell))

DimHeatmap(cell, dims = 1:15, cells = 500, balanced = TRUE)







