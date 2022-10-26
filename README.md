# ScRNA-Seq
TITLE- Expression of CD94 and CD96 on CD8+ T cells early after allogeneic stem cell transplantation is predictive of subsequent disease relapse and survival

GEO ACC NO- GSE185498

SAMPLE- 576

PLATFORM -Illumina NovaSeq 6000

# calling library

library(Matrix)

library(Seurat)

library(dplyr)

library("patchwork")

# read a file

ss_file= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/GSE185498_raw_counts_aggregated.csv")

**check duplicates**

length(unique(ss_file$Gene.names)) == nrow(ss_file)

**remove duplicate value**

ss_file=ss_file[!duplicated(ss_file$Gene.names),]

file <- ss_file[,-1]

rownames(file) <- ss_file[,1]

View(file)

**Initialize the Seurat object with the raw (non-normalized data)**

cell1= CreateSeuratObject(counts = file, min.cells = 3, min.features = 200) 

cell1

![image](https://user-images.githubusercontent.com/66779651/197959532-b62fb938-1b8d-490d-af55-a2f834ff00ec.png)

cell1[1:50, 1:10]

**The [[ operator can add columns to object metadata.**

cell1[["percent.mt"]] = PercentageFeatureSet(cell1, pattern = "^MT-")

head(cell1@meta.data)

![image](https://user-images.githubusercontent.com/66779651/197960033-3efda1b6-38a6-468a-84bb-8cc587d3e990.png)

VlnPlot(cell1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

![image](https://user-images.githubusercontent.com/66779651/197960317-7bab7d56-c984-4a37-81d5-d6a03769427e.png)

cell = subset(cell1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell

#  removing unwanted cells from the dataset, the next step is to normalize the data

cell = NormalizeData(cell)

cell = FindVariableFeatures(cell, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(cell), 10)
top10

![image](https://user-images.githubusercontent.com/66779651/197960446-ef2c02c5-d21c-40d2-8753-0381f766c808.png)

plot1 = VariableFeaturePlot(cell)

plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot2

![image](https://user-images.githubusercontent.com/66779651/197960585-75fbd49f-6233-4fda-84ac-88bf0cae08f3.png)

all.genes = rownames(cell)
cell = ScaleData(cell, features = all.genes)

## Centering and scaling data matrix

cell@assays$RNA@scale.data[1:50, 1:5]

cell = RunPCA(cell, features = VariableFeatures(object = cell))

![image](https://user-images.githubusercontent.com/66779651/197960841-e0fad3db-84a6-4843-9ac5-e67f1fd3d1fd.png)


DimHeatmap(cell, dims = 1:15, cells = 500, balanced = TRUE)
![image](https://user-images.githubusercontent.com/66779651/197960903-5c90320f-980c-4934-812f-db3c0f489994.png)
