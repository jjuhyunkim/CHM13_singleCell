# --------------------------------- Parameters ------------------------------------#
args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
mainDir <- args[2]

# filtering threshold
nCount_RNA_upperCutoff <- 60000
nFeature_RNA_underCutoff <- 100
percentMT_upperCutoff <- 25
percentRibo_underCutoff <- 10

# Doublet Finder Parameters 
cluter = ""
doublet_frac <- 0.1663 # 0.075 = 7.5%
nN_num = 0.25
pK_num = 0.09
souporcell_doublet_fraction=0

# --------------------------------- Load Liirary ---------------------------------- #
library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(ggplot2)
library(sctransform)
library(patchwork)
library(dplyr)
library(ggplot2)

# ---------------------------------- Simple Code ---------------------------------- # 
# make DoubletFinder dir
dir.create(mainDir,"/DoubletFinder")

# load data from the filtered h5 file
data.file <- paste0(mainDir,"/cellbender/",sample,".cellbender_filtered.h5")
data.data <- Read10X_h5(filename = data.file, use.names = TRUE)

# create Seurat object
sdata <- CreateSeuratObject(counts = data.data, project = sample)
sdata

## Simple QC ## 
# store mitochondrial percentage in object meta data
sdata <- PercentageFeatureSet(sdata, pattern = "^MT-", col.name = "percent.mt")
sdata <- PercentageFeatureSet(sdata, pattern = "^RP[SL]", col.name = "percent.ribo")

# Subsampleing
sdata <- subset(sdata, subset = nCount_RNA < nCount_RNA_upperCutoff & nFeature_RNA > nFeature_RNA_underCutoff & percent.mt < percentMT_upperCutoff & percent.ribo > percentRibo_underCutoff)
sdata

# run sctransform
sdata <- SCTransform(sdata, vars.to.regress = c("percent.mt"), verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
sdata <- RunPCA(sdata, verbose = FALSE)
sdata <- RunUMAP(sdata, dims = 1:30, verbose = FALSE)

sdata <- FindNeighbors(sdata, dims = 1:30, verbose = FALSE)
sdata <- FindClusters(sdata, verbose = FALSE, resolution = c(0.1,0.3,0.6,0.9))

#---------------------------------------------- Running DoubletFinder ## --------------------------------------------------
# Parameters 
cluter = ""
doublet_frac <- 0.1663 # 0.075 = 7.5%
nN_num = 0.25
pK_num = 0.09

Nonidentified_doublet_frac <- doublet_frac - souporcell_doublet_fraction
print(paste0("doublets fraction from DoubletFinder : ",Nonidentified_doublet_frac))

## pK Identification (no ground-truth) --------------------
sweep.res <- paramSweep_v3(sdata, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## pK Identification (ground-truth) --------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico eeneotyping results 
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate ---------------------
annotations <- sdata@meta.data[,cluster]
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(doublet_frac*nrow(sdata@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sdata <- doubletFinder_v3(sdata, PCs = 1:10, pN = nN_num, pK = pK_num, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
sdata <- doubletFinder_v3(sdata, PCs = 1:10, pN = nN_num, pK = pK_num, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

saveRDS(sdata, paste0(mainDir,"/DoubletFinder/",sample,".Cellbender.DoubletFinder.SCTransform.rds"))
SaveH5Seurat(sdata, filename = paste0(mainDir,"/DoubletFinder/",sample,".Cellbender.DoubletFinder.SCTransform.h5Seurat"))
Convert(paste0(mainDir,"/DoubletFinder/",sample,".Cellbender.DoubletFinder.SCTransform.h5Seurat"), dest = "h5ad")


