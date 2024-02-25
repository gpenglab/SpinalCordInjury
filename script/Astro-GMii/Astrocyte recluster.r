library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(clustree)
library(Matrix)
library(future)
library(grid)
library(RColorBrewer)
library(pheatmap)
library(igraph)
library(stringr)
library(cowplot)
library(scales)
library(reshape2)
library(clusterProfiler)
options(connectionObserver = NULL)
library("readxl")

spatial_marker<-read.csv("2subtype.marker.sig.enrich.in.sham.2domains.df.csv")
da<-readRDS("splitpool_3astr.labeltransfer.10x_L3.astrocyte.labelAstrcoyte.object.with7dpi.rds")

da2<-readRDS("sci.rds")
DefaultAssay(da2)<-"RNA"
test<-subset(da2,subset = celltype=="Astrocyte")
da_s<-CreateSeuratObject(counts = test@assays$RNA@counts,meta.data = test@meta.data,min.cells = 10)

rps_genes<-rownames(da_s@assays$RNA@counts)[grep("^Rps",rownames(da_s@assays$RNA@counts))]

rpl_genes<-rownames(da_s@assays$RNA@counts)[grep("^Rpl",rownames(da_s@assays$RNA@counts))]

mt_genes<-rownames(da_s@assays$RNA@counts)[grep("^mt-",rownames(da_s@assays$RNA@counts))]
# remove ribosome and mitocondrial genes
da_s<-da_s[!rownames(da_s@assays$RNA@counts)%in%c(rps_genes,rpl_genes,mt_genes),rownames(da_s@meta.data)]
da_s
### QC
# Visualize QC metrics as a violin plot
options(repr.plot.width=10,repr.plot.height=5)
VlnPlot(da_s, features = c("nFeature_RNA", "nCount_RNA", "percent_mt","percent_rp","percent_hbb"), ncol = 5)
#ggsave("qc.vlnplot.png",width = 10,height = 5,dpi = 300)

da_s<-SCTransform(da_s,do.scale = TRUE,
                  vars.to.regress = c("orig.ident","S.Score","G2M.Score","percent_mt","percent_rp"
                                     ),
                  residual.features = rownames(da_s@assays$RNA@counts),
                    return.only.var.genes = FALSE)
saveRDS(da_s,"L2.Astrocyte.rmMt_Rp.SCT.regICMR.doscale.rds")

cells<-rownames(da_s@meta.data)[da_s$nCount_RNA>2500]
test<-da_s[rownames(da_s@assays$RNA@counts),cells]
test<-FindVariableFeatures(test)
test<-RunPCA(test#,features = astr_marker
              )
options(repr.plot.width=10,repr.plot.height=4)
DimPlot(test,
        reduction = "pca",dims = c(1,2),group.by = c("time","sample_id"))
options(repr.plot.width=5,repr.plot.height=5)
ElbowPlot(test,ndims = 50)
test<-FindNeighbors(test,dims=1:16,k.param = 40)
for(i in c(0.1,0.125,0.15,0.2,0.3,0.35,0.4,0.45,0.5)){
    test<-FindClusters(test,resolution = i)
}

options(repr.plot.width=6,repr.plot.height=8)
clustree(test,prefix = "SCT_snn_res.")

test<-RunUMAP(test,dims = 1:16)

options(repr.plot.width=3.6,repr.plot.height=4)
DimPlot(test,reduction = "umap",#cols = DiscretePalette(n=10,palette = "glasbey"),#cells = cells,
        group.by = c("time"))+
     scale_color_manual(values = c('#374E55FF','#00A1D5FF','#B24745FF','#d6cbd3'
))
ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.time.UMAP.png",width = 3.6,height = 4,dpi = 300)

options(repr.plot.width=3,repr.plot.height=4)
DimPlot(test,reduction = "umap",#cols = DiscretePalette(n=10,palette = "glasbey"),#cells = cells,
        group.by = c("SCT_snn_res.0.2"))+
     scale_color_manual(values = c('#911eb4','#e6beff'))
ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.res02.2subtype.UMAP.png",width = 3,height = 4,dpi = 300)

pdf("lee.Astrocyte.Pla2g7_Slc7a10_Gfap_Igfbp2.UMAP.220820.pdf",width = 14,height = 5)
options(repr.plot.width=14,repr.plot.height=5)
FeaturePlot(test,pt.size = 1,ncol = 4,#cells = cells,
            features = c("Gfap","Pla2g7",
                         #"Slc1a2","Aldoc","Bcan","Clu","S100b","Aqp4",
                         "Slc7a10","Igfbp2"#,"Mdk","Itpkb"
                        ),
           cols= c(#"#cfe0e8",
    
                                                     # '#b8a9c9', 
    '#d6d4e0',
    "#f0f0f0",
    '#fff2df',
    "#eeac99",#"#ED797B",
                     "#c83349")
           )
dev.off()
#ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.Gfap_Slc7a10_Vim.data.UMAP.png",width = 12,height = 5,dpi = 300)

options(repr.plot.width=12,repr.plot.height=5)
FeaturePlot(test,pt.size = 1,ncol = 3,#cells = cells,
            features = c("Gfap","Slc7a10","Bcan"
                        ),
           cols= c(#"#cfe0e8",
    
                                                     # '#b8a9c9', 
    '#d6d4e0',
    "#f0f0f0",
    '#fff2df',
    "#eeac99",#"#ED797B",
                     "#c83349")
           )
#ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.Gfap_Slc7a10_Vim.data.UMAP.png",width = 12,height = 5,dpi = 300)

pdf("lee.Astrocyte.Pla2g7_Slc7a10_Gfap_Igfbp2.Vlnplot.220820.pdf",width = 10,height = 4)
options(repr.plot.width=10,repr.plot.height=4)
VlnPlot(test,group.by = "subtype",ncol = 5,cols = c('#911eb4','#e6beff'),
        features=c("Gfap","Slc7a10","Pla2g7","Igfbp2"),
        pt.size = 0)#+
        #scale_fill_manual(values = c('#911eb4','#e6beff'))
dev.off()
#ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.res02.some_genes.UMAP.png",width = 10,height = 6,dpi = 300)

###migration score
cm<-read_excel("GO_term_positiveregulationofcellmigration_20220620_052409.xlsx")
cmgenes<-unique(cm$Symbol)
cmgenes<-intersect(cmgenes,rownames(test@assays$SCT@data))
test<-AddModuleScore(test,features = list(cmgenes),name = "positive_cell_migration")
options(repr.plot.width=4,repr.plot.height=4)
VlnPlot(test,group.by = "SCT_snn_res.0.2",cols = c('#911eb4','#e6beff'),
        features=c("dp_marker1"),
        pt.size = 0)+
    geom_boxplot(width=0.25,fill="white",outlier.colour = NA)
ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.k40.res02.2subtype.GM_Gfap_both_pos_deg_score.UMAP.png",width = 3.6,height = 4,dpi = 300)

test$subtype<-"Astr_Gfap"
test$subtype[test$SCT_snn_res.0.2==1]<-"Astr_Slc7a10"

#Idents(test)<-test$subtype
options(repr.plot.width=4,repr.plot.height=4)
VlnPlot(test,group.by = "subtype",
        features=c("positive_cell_migration1"),
        pt.size = 0)+
    geom_boxplot(width=0.35,fill="white",outlier.colour = NA)+
    scale_fill_manual(values = c('#911eb4','#e6beff'))
ggsave("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.positive_cell_migration_score.across.subtype.vlnplot.png",width = 4,height = 4,dpi = 300)

pdf("astro.2subtypes.2marker.vlnplot_220913.pdf",width = 8,height = 4)
options(repr.plot.width=8,repr.plot.height=4)
VlnPlot(test,group.by = "subtype",ncol = 4,
        features=c("Pla2g7","Slc7a10","Gfap","Dhrs1"),
        pt.size = 0)
dev.off()

### subset gfap cells
test_g<-subset(test,subset=SCT_snn_res.0.2==0)
# remove small cluster(pc16, res01)
test_g<-subset(test_g,subset= SCT_snn_res.0.1 ==0)
test_g<-FindVariableFeatures(test_g)
test_g<-RunPCA(test_g#,features = astr_marker
              )
test_g<-FindNeighbors(test_g,dims=1:10,k.param = 20
                     )
for(i in c(0.1,0.125,0.15,0.2,0.3,0.35,0.4,0.45,0.5)){
    test_g<-FindClusters(test_g,resolution = i)
}
test_g<-RunUMAP(test_g,dims = 1:10)

saveRDS(test_g,"Astrocyte.Gfap.subset.pc10.res04.rds")







