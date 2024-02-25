library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(spacexr)
library(Matrix)
library(future)
library(grid)
library(RColorBrewer)
library(circlize)
library(stringr)
library(ComplexHeatmap)
library(tidyverse)
source("script/self_function/SPOTlightHelper-Copy1.R")
source("script/self_function/DotPlot_ct.R")
get_marker_data <- function(cell_type_names, cell_type_means, gene_list) {
  marker_means = cell_type_means[gene_list,]
  marker_norm = marker_means / rowSums(marker_means)
  marker_data = as.data.frame(cell_type_names[max.col(marker_means)])
  marker_data$max_epr <- apply(cell_type_means[gene_list,],1,max)
  colnames(marker_data) = c("cell_type",'max_epr')
  rownames(marker_data) = gene_list
  marker_data$log_fc <- 0
  epsilon <- 1e-9
  for(cell_type in unique(marker_data$cell_type)) {
    cur_genes <- gene_list[marker_data$cell_type == cell_type]
    other_mean = rowMeans(cell_type_means[cur_genes,cell_type_names != cell_type])
    marker_data$log_fc[marker_data$cell_type == cell_type] <- log(epsilon + cell_type_means[cur_genes,cell_type]) - log(epsilon + other_mean)
  }
  return(marker_data)
}

###Create RCTD spatial query
da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
counts<-as.data.frame(da@assays$SCT@counts)
coords<-da@meta.data[,c('x.axis',"y.axis")]
nUMI<-colSums(counts)
puck<-SpatialRNA(coords,counts,nUMI)
saveRDS(puck,"RCTD/WT.replace_v2.SCT.puck.rds")

###create RCTD reference
da3_filt<-readRDS("sc/splitpool_3astr.labeltransfer.10x_L3.astrocyte.with7dpi.rds")
#downsampling
set.seed(1234)
da3_filt_down<-subset(da3_filt,downsample=1000)
table(Idents(da3_filt_down))
options(repr.plot.width=15,repr.plot.height=15)
VlnPlot(da3_filt_down,c("Gfap","Slc7a10"),group.by = "labeltranfer",pt.size = 0,ncol = 1)
ggsave("2astr.labeltransfer.L3.astrocyte.marker.across.allcelltypes.with7dpi.downsample1000.vlnplot.png",width=15,height=15)
cell<-da3_filt_down$labeltranfer
cell<-as.factor(cell)
reference<-Reference(da3_filt_down@assays$RNA@counts,cell)
saveRDS(reference,"2reference.2astr.labeltransfer.L3.astrocyte_RCTD.reference.seed1234.down1000.with7dpi.rds")

###run RCTD
myRCTD <- create.RCTD(puck, reference,CELL_MIN_INSTANCE = 10, 
                      UMI_min = (min(puck@nUMI)-1),UMI_max = 200000000,fc_cutoff_reg = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS(myRCTD,'myRCTD_visium_multi.rds')

### get the cell type highest marker list and df
sc_exp<-myRCTD@cell_type_info$info[[1]]
sc_df<-get_marker_data(myRCTD@cell_type_info$info[[2]],myRCTD@cell_type_info$info[[1]],myRCTD@internal_vars$gene_list_bulk)
sc_df$gene<-rownames(sc_df)
sc_exp_norm<-sc_exp[sc_df$gene,]
sc_exp_norm<-sc_exp_norm / rowSums(sc_exp_norm)
sc_exp_norm$gene<-rownames(sc_exp_norm)
dev_marker_df<-cbind(sc_exp_norm,sc_df)
write.csv(dev_marker_df,"fcreg1.deconvolution.celltype.marker.normalized.data.frame.csv")

###get confident deconvolution dataframe
result<-myRCTD@results
sub_weights <- data.frame(result[[1]]$sub_weights)
sub_weights$ct <- row.names(sub_weights)
colnames(sub_weights)[1]<-1

sub_weights_conf <- data.frame(result[[1]]$sub_weights[result[[1]]$conf_list])
sub_weights_conf$ct <- row.names(sub_weights_conf)
colnames(sub_weights_conf)[1]<-1
for(i in seq(2, length(result))){
    tmp <- data.frame(result[[i]]$sub_weights[result[[i]]$conf_list])
    tmp$ct <- row.names(data.frame(tmp))
    colnames(tmp)[1]<-i
    sub_weights_conf <- merge(x = sub_weights_conf, y =tmp, by = "ct", all = T)
    
    tmp2 <- data.frame(result[[i]]$sub_weights)
    tmp2$ct <- row.names(data.frame(tmp2))
    colnames(tmp2)[1]<-i
    sub_weights <- merge(x = sub_weights, y =tmp2, by = "ct", all = T)
}
sub_weights[is.na(sub_weights)] <- 0
sub_weights_conf[is.na(sub_weights_conf)] <- 0
dec_mtx<-sub_weights
rownames(dec_mtx)<-dec_mtx$ct
dec_mtx<-dec_mtx[,-1]
dec_mtx.t<-t(dec_mtx)
dec_mtx_conf<-sub_weights_conf
rownames(dec_mtx_conf)<-dec_mtx_conf$ct
dec_mtx_conf<-dec_mtx_conf[,-1]
dec_mtx_conf.t<-t(dec_mtx_conf)
### normalize confident dec_mtx
empty_cells<-rownames(dec_mtx_conf.t)[rowSums(dec_mtx_conf.t)==0]
length(empty_cells)
dec_mtx_conf.t_nor<-normalize_weights(dec_mtx_conf.t)
dec_mtx_conf.t_nor[empty_cells,]<-0
table(rowSums(dec_mtx_conf.t_nor)==0)
rownames(dec_mtx.t)<-colnames(puck@counts)
rownames(dec_mtx_conf.t_nor)<-colnames(puck@counts)

rownames(dec_mtx.t)<-colnames(puck@counts)
rownames(dec_mtx_conf.t_nor)<-colnames(puck@counts)

colnames(dec_mtx.t)<-gsub("[[:punct:]]",".",colnames(dec_mtx.t))
colnames(dec_mtx.t)<-gsub(" ",".",colnames(dec_mtx.t))

colnames(dec_mtx_conf.t_nor)<-gsub("[[:punct:]]",".",colnames(dec_mtx_conf.t_nor))
colnames(dec_mtx_conf.t_nor)<-gsub(" ",".",colnames(dec_mtx_conf.t_nor))

write.csv(dec_mtx.t, "WT.fcreg1_subweights.decon.mtx.csv")
write.csv(dec_mtx_conf.t_nor, "WT.fcreg1_subweights_conf_normalized.decon.mtx.csv")

### merged into 16 cell types
df_merged<-as.data.frame(dec_mtx_conf.t_nor)
df_merged$Vascular<-rowSums(df_merged[,c("A.Endothelial","C.Endothelial",'Pericyte','Tip.Cell','U.Vascular','V.Endothelial','VSMC')])
df_merged$OPC<-rowSums(df_merged[,c('Div.OPC','OPC.A','OPC.B','Pre.Oligo')])
df_merged$Macrophage<-rowSums(df_merged[,c('Border.Associated.Mac','Chemotaxis.Inducing.Mac','Inflammatory.Mac')])
#H_Microglia<-apply(df_merged[,c('Homeostatic.Microglia')],1,sum)
df_merged$R.Microglia<-rowSums(df_merged[,c('Dividing.Microglia','Inflammatory.Microglia','Migrating.Microglia')])
df_merged$Myeloid<-rowSums(df_merged[,c('Dividing.Myeloid','Interferon.Myeloid')])
df_merged$Ependymal<-rowSums(df_merged[,c('Astroependymal','Ependymal.A','Ependymal.B')])
###220724 merge dendritic and Astr-Svep1 into others
df_merged$others<-rowSums(df_merged[,c('Dendritic','Astro.Svep1')])
df_merged<-df_merged[,-which(colnames(df_merged)%in%c("A.Endothelial","C.Endothelial",'Pericyte','Tip.Cell','U.Vascular','V.Endothelial','VSMC',
                                                        'Div.OPC','OPC.A','OPC.B','Pre.Oligo',
                                                        'Border.Associated.Mac','Chemotaxis.Inducing.Mac','Inflammatory.Mac',
                                                        #'Homeostatic.Microglia',
                                                        'Dividing.Microglia','Inflammatory.Microglia','Migrating.Microglia',
                                                        'Dividing.Myeloid','Interferon.Myeloid',
                                                        'Astroependymal','Ependymal.A','Ependymal.B',
                                                      'Dendritic','Astro.Svep1'
                                                     ))]
write.csv(df_merged,"WT.fcreg1_subweights.decon.conf_nor.merged15celltypes.csv")
meta<-read.csv("WT_replace_v2/res02_220310/WT.SCT.pc20.k50.res02.meta.data.csv")
rownames(meta)<-meta[,1]
meta_dec<-cbind(meta,df_merged)
write.csv(meta_dec,"15celltypes.dec_conf_nor.meta.csv")






