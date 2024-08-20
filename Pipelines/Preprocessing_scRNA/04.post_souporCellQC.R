library(dplyr)
library(ggplot2)
library(stringr)
library(Matrix)
library(DropletUtils)
library(Seurat)

argv <- commandArgs(TRUE)
workdir <- argv[1]
batch <- argv[2]
file_10x <- argv[3]
metadata <- argv[4]
loadedcell <- argv[5]
pooled <- argv[6]

sdata.data=Read10X(data.dir=file_10x)
sdata=CreateSeuratObject(counts=sdata.data, project=batch)
rm(sdata.data)

meta <- read.table(metadata, header = T, row.names=1)[,c(1,2)]
sdata <- AddMetaData(sdata, metadata=meta)

sdata <- subset(sdata,status!="unassigned")

sdata <- PercentageFeatureSet(sdata, pattern = "^MT-", col.name = "percent.mt")
sdata <- SetIdent(sdata,value = sdata@meta.data$assignment)
sdata@meta.data$filter=ifelse(sdata@meta.data$status=="doublet","Doublet","Singlet")

### repeat for each assignment
pdf(paste0(workdir,"/",batch,".souporcell.identity.prefilter.pdf"))
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), group.by = "orig.ident")
VlnPlot(sdata, features = "nFeature_RNA", y.max=3000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", y.max=20000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", y.max=5000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "percent.mt", pt.size=0,group.by = "orig.ident")
VlnPlot(sdata, features = "percent.mt", y.max=25, pt.size=0,group.by = "orig.ident")
dev.off()

saveRDS(sdata , paste0(batch,"original.RDS"))
