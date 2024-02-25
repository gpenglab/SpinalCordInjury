library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(ggalluvial)
library(Matrix)
library(future)
library(grid)
library(RColorBrewer)
library(circlize)
library(igraph)
library(stringr)
library(cowplot)
library(clusterProfiler)
options(connectionObserver = NULL)
library(biomaRt)

### load 10X reference
me<-read.csv("sc/obs_sci.csv")
da<-readRDS("sc/sci.rds")
DefaultAssay(da)<-"RNA"
cell<-intersect(me$X,colnames(da))
rownames(me)<-me$X
me<-me[,-1]
me<-me[colnames(da),]
table(colnames(da)==rownames(me))
#merge neuron subtypes
neuron<-rownames(da@meta.data[da@meta.data$celltype=="Neuron",])
me[neuron,"celltype"]<-"Neuron"
me[neuron,"L1_taxon"]<-"Neuron"
me[neuron,"L2_taxon"]<-"Neuron"
me[neuron,"L3_taxon"]<-"Neuron"
#remove NA cells
na_cell<-rownames(me[is.na(me$L3_taxon),])
me_rm<-me[!rownames(me)%in%na_cell,]
da_rm<-da[rownames(da),rownames(me_rm)]
da_rm@meta.data<-me_rm
da_ad<-da_rm
saveRDS(da_rm,"10X.L3_taxon.na.processed.rds")
da_rm<-CreateSeuratObject(counts = da_rm@assays$RNA@counts,meta.data = da_rm@meta.data)
#subset astrocyte
da1_astr<-da_rm[rownames(da_rm),rownames(da_rm@meta.data[da_rm@meta.data$L3_taxon=="Astrocyte",])]
table(da1_astr$L3_taxon)
saveRDS(da1_astr,"10X.L3.Astrcoyte.RNA.rds")
### load split-pool reference
ref <- readRDS("sc/ref.spine.with.unresolved.seurat.clustered.rds")
#remove Unknown cells and rename cell types
ref_rmUnknown<-ref[rownames(ref),
          rownames(ref@meta.data[!ref@meta.data$annotation %in% 
                                               c('32 Unresolved                                      ',
                                                '33 Unresolved                                      ',
                                                '2 Unassigned                                       ',
                                                '3 Unassigned                                       '#,
                                                #'15 Cerebrospinal Fluid-Contacting Neurons (CSF-cNs)',
                                                #'4 Astroctye - Unassigned                           ',
                                                #'14 Committed Oligodendroctye Precursor Cells       '
                                                ),])]
cell1_meta<-ref_rmUnknown@meta.data
cell1_meta[grepl("Excitatory",cell1_meta$annotation,fixed = T),"annotation"]<-"Excitatory_neuron"

cell1_meta[grepl("Inhibitory",cell1_meta$annotation,fixed = T),"annotation"]<-"Inhibitory_neuron"
cell1_meta[grepl("motor neurons",cell1_meta$annotation,fixed = T),"annotation"]<-"motor_neuron"
cell1_meta[grepl("5 Astrocyte - Gfap                                 ",cell1_meta$annotation,fixed = T),"annotation"]<-"Astrocyte_Gfap"
cell1_meta[grepl("13 OPC                                             ",cell1_meta$annotation,fixed = T),"annotation"]<-"OPC"
cell1_meta[grepl("7 Astro - Svep1                                    ",cell1_meta$annotation,fixed = T),"annotation"]<-"Astro_Svep1"
cell1_meta[grepl("11 Oligo Mature                                    ",cell1_meta$annotation,fixed = T),"annotation"]<-"Oligo_Mature"
cell1_meta[grepl("15 Cerebrospinal Fluid-Contacting Neurons (CSF-cNs)",cell1_meta$annotation,fixed = T),"annotation"]<-"CSF-cNs"
cell1_meta[grepl("14 Committed Oligodendroctye Precursor Cells       ",cell1_meta$annotation,fixed = T),"annotation"]<-"COPC"
cell1_meta[grepl("4 Astroctye - Unassigned                           ",cell1_meta$annotation,fixed = T),"annotation"]<-"Astroctye_Unassigned"
cell1_meta[grepl("12 Oligodendroctye Myelinating                     ",cell1_meta$annotation,fixed = T),"annotation"]<-"Oligodendroctye_Myelinating"
cell1_meta[grepl("8 Endothelial                                      ",cell1_meta$annotation,fixed = T),"annotation"]<-"Endothelial"
cell1_meta[grepl("6 Astrocyte - Slc7a10                              ",cell1_meta$annotation,fixed = T),"annotation"]<-"Astrocyte_Slc7a10"
cell1_meta[grepl("9 VLMC                                             ",cell1_meta$annotation,fixed = T),"annotation"]<-"VLMC"
cell1_meta[grepl("10 Microglia                                       ",cell1_meta$annotation,fixed = T),"annotation"]<-"Microglia"
cell1_meta[grepl("1 Epyndemal                                        ",cell1_meta$annotation,fixed = T),"annotation"]<-"Ependymal"
ref_rmUnknown@meta.data<-cell1_meta
cell<-ref_rmUnknown$annotation
cell<-as.factor(cell)
da2<-Reference(ref_rmUnknown@assays$RNA@counts,cell)
saveRDS(da2,"sc.reference.ALS.rmUnkonwn.rds")
meta2<-data.frame(da2@cell_types,row.names = colnames(da2@counts))
colnames(meta2)[1]<-"cell_types"
da_s<-CreateSeuratObject(counts = da2@counts,meta.data = meta2
                        )
da_s<-NormalizeData(da_s)
Idents(da_s)<-da_s$cell_types
### signatures of cell types
cell_marker<-FindAllMarkers(da_s,logfc.threshold = 1,only.pos = TRUE)
write.csv(cell_marker,"Rosenberg.celltype_marker.fc1.pos.csv")
#subsert astrocyte
da2_astr<-da_s[rownames(da_s),rownames(da_s@meta.data[da_s@meta.data$cell_types %in% c('Astrocyte_Gfap','Astrocyte_Slc7a10',"Astro_Svep1"),])]
saveRDS(da2_astr,"split-pool.astrocyte.rds")

###LabelTransfer Rosenberg astrocyte subtype to Milich
query<-da1_astr
query<-NormalizeData(query,verbose = F)
query<-FindVariableFeatures(query,selection.method = "vst")
query<-ScaleData(query,verbose=F)
ref<-da2_astr
ref<-NormalizeData(ref,verbose = F)
ref<-FindVariableFeatures(ref,selection.method = "vst")
ref<-ScaleData(ref,verbose = F)
ref<-RunPCA(ref,npcs = 30,verbose = F)
ref<-RunUMAP(ref,reduction = "pca",dims = 1:30)
anchors<-FindTransferAnchors(reference = ref,query = query,dims = 1:30)
predictions<-TransferData(anchorset=anchors,refdata=ref$cell_types,dims = 1:30)
query<-AddMetaData(query,metadata = predictions)
options(repr.plot.width=10,repr.plot.height=6)
VlnPlot(query,c("Gfap","Slc7a10","Svep1"),group.by = "predicted.id",pt.size = 0)
ggsave("splitpool_3astr.labeltransfer.10x_L3.astr.markers.vlnplot_220509.png",width = 10,height = 6)
cell1<-rownames(query@meta.data[query@meta.data$predicted.id=="Astro_Svep1",])
cell2<-rownames(query@meta.data[query@meta.data$predicted.id=="Astrocyte_Gfap",])
cell3<-rownames(query@meta.data[query@meta.data$predicted.id=="Astrocyte_Slc7a10",])

##reassign cell label in Milich
da3<-da_rm[rownames(da2_rm),#rownames(da_rm@meta.data[!da_rm@meta.data$time=="7dpi"])
        ]
da3$labeltranfer<-da3$L3_taxon
da3@meta.data[cell1,"labeltranfer"]<-"Astro_Svep1"
da3@meta.data[cell2,"labeltranfer"]<-"Astrocyte_Gfap"
da3@meta.data[cell3,"labeltranfer"]<-"Astrocyte_Slc7a10"
rmcell<-names(da3$labeltranfer[!is.na(da3$labeltranfer)])
da3_filt<-da3[rownames(da3),rmcell]
Idents(da3_filt)<-da3_filt$labeltranfer
options(repr.plot.width=15,repr.plot.height=15)
VlnPlot(da3_filt,c("Gfap","Slc7a10"),group.by = "labeltranfer",pt.size = 0,ncol = 1)
ggsave("splitpool_3astr.labeltransfer.10x_L3.astrocyte.marker.across.allcelltypes.vlnplot_220509.png",width=15,height=15)
saveRDS(da3_filt,"splitpool_3astr.labeltransfer.10x_L3.astrocyte.with7dpi.rds")
## subset astrocyte object
da3_astr<-da3_filt[rownames(da3_filt),rownames(da3_filt@meta.data[da3_filt@meta.data$labeltranfer %in% c('Astrocyte_Slc7a10','Astrocyte_Gfap','Astro_Svep1'),])]
saveRDS(da3_astr,"splitpool_3astr.labeltransfer.10x_L3.astrocyte.labelAstrcoyte.object.with7dpi.rds")





