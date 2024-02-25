library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(future)
library(grid)
library(RColorBrewer)
library(pheatmap)
library(igraph)
library(stringr)
library(cowplot)
library(scales)
library(reshape2)
library(grDevices)
library(tradeSeq)
source("script/self_function/save_pheatmap_pdf.R")

knitr::opts_chunk$set(fig.align="center", cache=TRUE,error=FALSE, #stop on error
fig.width=5, fig.height=5, autodep=TRUE,
results="markup", echo=TRUE, eval=TRUE)
#knitr::opts_knit$set(stop_on_error = 2L) #really make it stop
#knitr::dep_auto()
options(getClass.msg=FALSE)
graphics:::par(pch = 16, las = 1)
set.seed(12345) ## for reproducibility
library(SingleCellExperiment)

#Load single cell astrocyte data
da<-readRDS("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.rds")

da2<-readRDS("Astrocyte.Gfap.subset.pc10.res04.rds")
cells<-colnames(da2)

da$gfap_hclust[is.na(da$gfap_hclust)]<-"GM_Slc7a10"
da$gfap_hclust[da$gfap_hclust==1]<-"WM_Gfap"
da$gfap_hclust[da$gfap_hclust==2]<-"GM_Gfap"

rd<-as.data.frame(da@reductions$umap@cell.embeddings)
cells<-rownames(rd)[rd$UMAP_1>-4 & rd$UMAP_2<4]
options(repr.plot.width=4,repr.plot.height=4)
DimPlot(da,reduction = "umap",#cols = DiscretePalette(n=10,palette = "glasbey"),
        cells = cells,
        group.by = c("gfap_hclust"))+
     scale_color_manual(values = c('#ff7b25','#e6beff','#911eb4'))

rd<-da[,cells]@reductions$umap@cell.embeddings
cl<-da[,cells]$gfap_hclust

sce<-SingleCellExperiment(assays=List(counts=as.matrix(da@assays$SCT@counts[VariableFeatures(da),cells])))


#1.gene filtering
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 1) >= 10
})
sce <- sce[geneFilter, ]

#2.assign normalized data
assays(sce)$norm <- as.matrix(da@assays$SCT@data[rownames(sce),cells])#FQnorm(assays(sce)$counts)
plot(rd, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce) <- SimpleList(UMAP = rd)
colData(sce)$hcluster <- cl
#pdf("Gfap.astrocyte.subset.hcluster.UMAP.pdf")
#plot(rd, col = c('#ff7b25',
#                 #'#e6beff',
#    '#911eb4')[as.numeric(as.factor(cl))], pch=16, asp = 1)
#dev.off()

#3.using slingshot
sce <- slingshot(sce, clusterLabels = 'hcluster', reducedDim = 'UMAP')
colors <- colorRampPalette(brewer.pal(n = 9,name = "YlGn"))(100)#viridis_pal(option = "D")(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

options(repr.plot.width=5,repr.plot.height=5)
#pdf("Gfap.astrocyte.subset.trajectory.UMAP.pdf")
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=3.5, col='black')
#dev.off()

#identifying temporally dynamic genes
# fit negative binomial GAM
sce <- fitGAM(sce)
# test for dynamic expression
ATres <- associationTest(sce)
write.csv(ATres[ATres$meanLogFC>1.5,],"slingshot.dynamic.gene.logfc1.5.csv")

topgenes <- rownames(ATres[order(ATres$pvalue), ])[50:1]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- as.matrix(assays(sce)$norm)[topgenes, pst.ord]
heatpseud <- sce$slingPseudotime_1[pst.ord]
state<- sce$hcluster[pst.ord]
time<-da[,cells]$time[pst.ord]
col_df<-data.frame(row.names = names(state),state=state,time=time,pseud=heatpseud)
time<-c('#00A1D5FF','#B24745FF','#d6cbd3','#374E55FF')
names(time)<-c("1dpi","3dpi","7dpi","Uninjured")
state<-c('#911eb4','#ff7b25')
names(state)<-c("WM_Gfap","GM_Gfap")
an_col<-list(time=time,state=state)

options(repr.plot.width=10,repr.plot.height=8)
p<-pheatmap(heatdata,clustering_method = "ward.D",
         #color = viridis_pal(option = "F")(100),
         annotation_colors = an_col,
         annotation_col = col_df,
         #annotation_row = col_ha,
         cluster_cols = FALSE,show_colnames = FALSE,
         scale = "row",
         color = colorRampPalette(colors = c("#d6d4e0","#b8a9c9","#622569","#ffcc5c","#ffcc5c"))(50)
        )
#save_pheatmap_pdf(p,"Gfap.astrocyte.subset.trajectory.top50.genes.pheatmap.pdf")


