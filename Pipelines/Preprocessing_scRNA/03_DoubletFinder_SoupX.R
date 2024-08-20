library(dplyr)
library(ggplot2)
library(stringr)
library(Matrix)
library(DoubletFinder)
library(DropletUtils)
library(Seurat)

argv <- commandArgs(TRUE)
workdir <- argv[1]
batch <- argv[2]
file_10x <- argv[3]
metadata <- argv[4]
loadedcell <- as.numeric(argv[5])
pooled <- as.numeric(argv[6])

under_nFeature <- as.numeric(argv[7])
upper_nFeature <- as.numeric(argv[8])
under_nCount <- as.numeric(argv[9])
upper_nCount <- as.numeric(argv[10])
upper_percent.mt <- as.numeric(argv[11])
souporcell_ambientRNA <- as.numeric(argv[2])

sdata_original <- argv[13]

print(sdata_original)

sdata <- readRDS(sdata_original)
sdata

### QC filter
sdata <- subset(sdata, subset = percent.mt < upper_percent.mt & nCount_RNA < upper_nCount & nCount_RNA > under_nCount & nFeature_RNA < upper_nFeature & nFeature_RNA > under_nFeature)
sdata

### repeat for each assignment

pdf(paste0(workdir,"/",batch,".after.manualQC.souporcell.identity.prefilter.pdf"))
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), group.by = "orig.ident")
VlnPlot(sdata, features = "nFeature_RNA", y.max=3000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", y.max=20000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "nCount_RNA", y.max=5000, pt.size=0, group.by = "orig.ident")
VlnPlot(sdata, features = "percent.mt", pt.size=0,group.by = "orig.ident")
VlnPlot(sdata, features = "percent.mt", y.max=25, pt.size=0,group.by = "orig.ident")
dev.off()


# ## SCT normalize for doubletfinder
### if needed, use cell cycle data as additional vars.to.regress
sdata <- SCTransform(sdata, vars.to.regress = c("percent.mt"))

sdata <- RunPCA(sdata)

pca = sdata@reductions$pca
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
pc_num <- 0
pc_cul <- 0
for (i in 1:length(varExplained)) {
  pc_cul <- pc_cul + varExplained[i]
    pc_num <- pc_num + 1
  if (pc_cul > 0.9) {
    break
  }
}
rm(pca,pc_cul,i,eigValues,varExplained)

pc_num <- max(20,min(pc_num,40))

sdata <- FindNeighbors(sdata, dims = 1:pc_num)
sdata <- FindClusters(sdata,resolution = 0.1,algorithm="louvain")
sdata <- RunUMAP(sdata, dims = 1:pc_num)

### visualize
pdf(paste0(workdir,"/",batch,".UMAP.cluster.pdf"))
DimPlot(sdata, reduction = 'umap', group.by = "SCT_snn_res.0.1", label = TRUE) + NoLegend() + ggtitle("UMAP leiden res 0.1") +theme(plot.title=element_text(hjust=0.5))
DimPlot(sdata, reduction = 'umap', group.by = "assignment") + ggtitle("souporcell identity") +theme(plot.title=element_text(hjust=0.5))
DimPlot(sdata, reduction = 'umap', group.by = "status") + ggtitle("souporcell call") +theme(plot.title=element_text(hjust=0.5))
dev.off()

### doubletfinder
sweep.res.list <- paramSweep_v3(sdata, PCs = 1:pc_num, sct = TRUE)
gt.calls=sdata@meta.data[rownames(sweep.res.list[[1]]),which(colnames(sdata@meta.data)=="filter")] ### make ground truth data from souporcell doublet calls
sweep.stats <- summarizeSweep(sweep.res.list, GT=TRUE, GT.calls = gt.calls)
bcmvn <- find.pK(sweep.stats)

### bcmvn visualize
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

pdf(paste0(workdir,"/",batch,".BCmvn.pdf"))
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b", col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
dev.off()

### @@@@@@@ change to adequate number @@@@@@
nonidentifiedmultiplet=loadedcell*0.000004597701/pooled

homotypic.prop <- modelHomotypic(sdata@meta.data$seurat_clusters)
nExp_poi <- round(nonidentifiedmultiplet*length(sdata@meta.data$seurat_clusters))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

sdata <- doubletFinder_v3(sdata, PCs = 1:pc_num, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
sdata <- doubletFinder_v3(sdata, PCs = 1:pc_num, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_choose,"_",nExp_poi), sct = TRUE)

### doubletfinder visualize
pdf(paste0(workdir,"/",batch,".DF.classification.pdf"))
DimPlot(sdata, reduction = 'umap', group.by = paste0("DF.classifications_0.25_",pK_choose,"_",nExp_poi), label = FALSE) + ggtitle("DF homotypic poisson UNadjusted") +theme(plot.title=element_text(hjust=0.5))
DimPlot(sdata, reduction = 'umap', group.by = paste0("DF.classifications_0.25_",pK_choose,"_",nExp_poi.adj), label = FALSE) + ggtitle("DF homotypic poisson adjusted") +theme(plot.title=element_text(hjust=0.5))
dev.off()

### write doubletfinder result
write.table(data.frame("##"=rownames(sdata@meta.data),sdata@meta.data),paste0(workdir,"/",batch,".human.filter1.backup.meta.tsv"),quote=FALSE,sep ='\t',row.names =FALSE,na="NA")
saveRDS(sdata,paste0(workdir,"/",batch,".filter1.backup.RDS"))


### filter doubletfinder and souporcell doublets and save as 10x format for soupX
### jaeyong used nExp_poi to filter doublets not nEXp_poi_adj
DF.result=FetchData(sdata,paste0("DF.classifications_0.25_",pK_choose,"_",nExp_poi))
sdata <- sdata[,which(DF.result=="Singlet")]
sdata=subset(sdata,filter!="Doublet")

### visualize filtered data with no additional normalization
pdf(paste0(workdir,"/",batch,".UMAP.cluster.doubletfiltered.pdf"))
DimPlot(sdata, reduction = 'umap', group.by = "SCT_snn_res.0.1", label = TRUE) + NoLegend() + ggtitle("UMAP leiden res 0.1") +theme(plot.title=element_text(hjust=0.5))
DimPlot(sdata, reduction = 'umap', group.by = "assignment") + ggtitle("souporcell identity") +theme(plot.title=element_text(hjust=0.5))
DimPlot(sdata, reduction = 'umap', group.by = "status") + ggtitle("souporcell call") +theme(plot.title=element_text(hjust=0.5))
dev.off()

### make metadata for soupX input
write.table(data.frame("##"=rownames(sdata@meta.data),sdata@meta.data),paste0(workdir,"/",batch,".human.soupx.meta.tsv"),quote=FALSE,sep ='\t',row.names =FALSE,na="NA")
write10xCounts(sdata@assays$RNA@counts, path=paste0(workdir,"/",batch,".doublet_filtered_matrix"), type="sparse", overwrite=TRUE,version="3")

saveRDS(sdata,paste0(workdir,"/",batch,".filter2.backup.RDS"))

library("SoupX")
alldroplet=Read10X(paste0(file_10x,"/../raw_feature_bc_matrix/"))
cluster=read.table(paste0(workdir,"/",batch,".human.soupx.meta.tsv"),header=T,row.names = "X..")
celldroplet=alldroplet[,colnames(alldroplet) %in% rownames(cluster)]

soupdata = SoupChannel(alldroplet,celldroplet)
rm(alldroplet,celldroplet)

soupdata = setClusters(soupdata,cluster$seurat_clusters)

### not used: soupX contamination estimation
###soupdata = autoEstCont(soupdata, tfidfMin = 1,soupQuantile = 0.9)

### get souporcell ambient RNA percent
soupdata$metaData$rho = souporcell_ambientRNA

### remove ambient RNA
out = adjustCounts(soupdata,roundToInt=TRUE)

### write matrix file
sdata <- CreateSeuratObject(counts=out)
saveRDS(sdata, paste0(workdir,"/",batch,".DoubletFinder.SoupX.matrix.RDS"))
writeMM(out,file=paste0(workdir,"/",batch,".DoubletFinder.SoupX.matrix.mtx"))
