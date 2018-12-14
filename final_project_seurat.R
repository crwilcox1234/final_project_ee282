#Seurat tsne
#library(devtools)
#install.packages('Seurat')
#install.packages("ggplot2")
library(ggplot2)
library(Seurat) # version 2.3.4
#devtools::install_github("tidyverse/ggplot2")
# The easiest way to get ggplot2 is to install the whole tidyverse:
#install.packages("tidyverse") # This way worked, must have tidyverse for ggplot2, which is needed for Seurat to install properly
library(dplyr)
library(Matrix)
setwd("D:/DATA/UC Irvine/JJEMERSON_CLASS_2018")
data <- read.table("final_matrix.txt",header = TRUE, row.names ="gene")
rnames <- rownames(data)
#rnames <- data[,0]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,1:ncol(data)])
data_create <- CreateSeuratObject(raw.data = mat_data, min.cells = 3, min.genes =2000)
data_create <- NormalizeData(object = data_create, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
length(data_create)
#filtered mitochondrial reads out of matrix already, so dont need to look at # mito reads
#mito.genes <- grep(pattern = "^MT-", x = rownames(x = data_create@data), value = TRUE)
#percent.mito <- Matrix::colSums(data_create@raw.data[mito.genes, ])/Matrix::colSums(data_create@raw.data)
#data_create <- AddMetaData(object=data_create, metadata= percent.mito, col.name="percent.mito")
#VlnPlot(object = plot, features.plot = c("nGene","nUMI","percent.mito"), nCol = 3)
data_create <- FindVariableGenes(object = data_create, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
data_create <- ScaleData(object = data_create)
length(data_create@var.genes)

data_create <- RunPCA(object = data_create, pc.genes = data_create@var.genes, pcs.compute = 20, do.print = TRUE, pcs.print = 1:5, 
                      genes.print = 5)

PCAPlot(object = data_create, dim.1 = 1, dim.2 = 2)

#PCAPlot(object = Normalize, dim.1 = 1, dim.2 = 2) 

data_create <- JackStraw(object = data_create, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = data_create, PCs = 1:13)
PCElbowPlot(object =  data_create)
data_create <- FindClusters(object = data_create, reduction.type = "pca", dims.use = 1:10, 
                            resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = data_create)
data_create <- RunTSNE(object = data_create, dims.use = 1:10, do.fast = TRUE,perplexity = 29)
TSNEPlot(object = data_create, do.label = TRUE) # If you dont want to label the clusters with numbers, make do.label = FALSE
saveRDS(data_create, file = "TZB_redo_11-20-18allcells.rds")