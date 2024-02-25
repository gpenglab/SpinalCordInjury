library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
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

da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
meta_dec<-read.csv("15celltypes.dec_conf_nor.meta.csv")
rownames(meta_dec)<-meta_dec[,1]
meta_dec<-meta_dec[,-1]

###quantify cell type each timepoint and position
celltypes<-c('Neuron','Oligodendrocyte','OPC','Homeostatic.Microglia','R.Microglia',
             'Macrophage','Myeloid','Monocyte','Neutrophil',
             'Fibroblast','Ependymal','Vascular',
                'Astrocyte.Gfap','Astrocyte.Slc7a10',"others"
    #'Astro.Svep1'
            )
# aggregate celltype proportion in each timepoint
temp<-meta_dec[,c("time",celltypes)]
temp$time<-factor(temp$time,levels = c("WT_sham","WT_3h","WT_24h","WT_72h"))
temp<-aggregate(temp[,celltypes],by = list(time=temp$time),FUN = sum)
temp<-meta_dec[,c("time",celltypes)]
temp$time<-factor(temp$time,levels = c("WT_sham","WT_3h","WT_24h","WT_72h"))
temp<-aggregate(temp[,celltypes],by = list(time=temp$time),FUN = sum)
rownames(temp)<-temp[,1]
temp<-temp[,-1]
write.csv(temp,"15celltypes.time.total_number.220825.csv")

###cell type spatial distribution
sf<-"15celltypes_spatial/"
if(!dir.exists(sf))
    dir.create(sf)
options(repr.plot.width=20,repr.plot.height=16)
for(i in celltypes){
    ggplot(meta_dec,aes_string("x.ad","y.ad",color=i))+
    geom_point(size=1)+
    scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(50))+
    xlab(paste0(""))+
    ylab(paste0(""))+
    theme(panel.background = element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
         )
    ggsave(paste0(sf,i,".spatialplot.png"),width=20,height=16,dpi=300)
}

pdf("Neuron_221011.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(meta_dec,aes_string("x.ad","y.ad",color="Neuron"))+
    geom_point(size=1)+
    scale_color_gradientn(colours = c('#f0efef','#c94c4c'))+
    xlab(paste0(""))+
    ylab(paste0(""))+
    theme(panel.background = element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
         )
dev.off()

pdf("oligodendrocyte_221011.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(meta_dec,aes_string("x.ad","y.ad",color='Oligodendrocyte'))+
    geom_point(size=1)+
    scale_color_gradientn(colours = c('#f0efef','#c94c4c'))+
    xlab(paste0(""))+
    ylab(paste0(""))+
    theme(panel.background = element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
         )
dev.off()

#pdf("Astrocyte.Slc7a10_221011.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(meta_dec,aes_string("x.ad","y.ad",color='Astrocyte.Slc7a10'))+
    geom_point(size=1)+
    scale_color_gradientn(colours = c('#f0efef','#c94c4c'))+
    xlab(paste0(""))+
    ylab(paste0(""))+
    theme(panel.background = element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
         )
#dev.off()

#pdf("Astrocyte.Gfap_221011.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(meta_dec,aes_string("x.ad","y.ad",color='Astrocyte.Gfap'))+
    geom_point(size=1)+
    scale_color_gradientn(colours = c('#f0efef','#c94c4c'),limits=c(0,0.5),na.value = "#c94c4c"
                         )+
    xlab(paste0(""))+
    ylab(paste0(""))+
    theme(
        panel.background = element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
         )
#dev.off()

da@meta.data<-meta_dec
saveRDS(da,"WT.merge.SCT_RCTD.c5.rds")

col<-c('#000075','#008080','#ffe119','#ffd8b1','#4363d8','#e6194b','#f58231','#808000','#fffac8','#46f0f0',
      '#800000','#aaffc3','#911eb4','#e6beff','#fabebe')
names(col)<-celltypes

dir.create("spatialPie_15celltypes")
slice <- Images(wt)
options(repr.plot.width = 15, repr.plot.height = 15)
for(i in 1:length(slice)){
    p1 <- SpatialScatterPie(wt,
                  cell_types_all = celltypes,#cell_types_all,
                  images=slice[i],       
                  cols = col,
                  pie.scale = 0.4)
    p1
    ggsave(paste("spatialPie_15celltypes/WT.deconvByRCTD_spatialPie.",slice[i],".png", sep = ""), p1, 
           width = 15, height = 15,dpi=300)
}

### remove others
celltype_ad<-celltypes[celltypes!="others"]
celltype_ad

###dotplot at each time and position
da$orig.ident<-factor(da$orig.ident,levels = c('WT_sham_H_R2_1mm',
                                               'WT_sham_H_R2',
                                               'WT_sham_T_210323',
                                               'WT_sham_T_R2_1mm',
                                               'WT_3h_H_R2_1mm',
                                               'WT_3h_H_R2',
                                               'WT_3h_T_R2',
                                               'WT_3h_T_210330_1mm',
                                               'WT_24h_H_R1_1mm',
                                               'WT_24h_H_R1',
                                               'WT_24h_T_201231',
                                               'WT_24h_T_R1_1mm',
                                               'WT_72h_H_R1_1mm',
                                               'WT_72h_H_210323',
                                               'WT_72h_T_R2',
                                               'WT_72h_T_R1_1mm'))

pdf("WT_merge_replace_v2.SCT.reg6.RCTD.fcreg1.14celltype_section.Dotplot_221011.pdf", width = 6, height = 5.5)
p3 <- DotPlot(da, features = celltype_ad, group.by = 'orig.ident', scale = F,cols = c('#66C2A5','#FC8D62'),dot.min = 1e-5) + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 12, margin = margin(10,0,0,0)))
options(repr.plot.width = 6, repr.plot.height = 5.5)
p3
dev.off()
#WM
da_wm<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="WM",])]
p3 <- DotPlot(da_wm, features = celltype_ad, group.by = 'orig.ident', scale = F,dot.min = 1e-5,cols = c('#66C2A5','#FC8D62')) + 
    theme_bw()+
    scale_colour_gradientn(colours=c('#66C2A5','#FC8D62'),breaks=c(0,0.4,0.8),limits=c(0,0.8)
                           )
options(repr.plot.width = 6.5, repr.plot.height = 4.5)
p3
ggsave("WT_merge_replace_v2.SCT.reg6.RCTD.fcreg1.14celltype_WM.ident.Dotplot.png", width = 6.5, height = 4.5,dpi = 300)
#MG
da_wm<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="MG",])]
p3 <- DotPlot(da_wm, features = celltype_ad, group.by = 'orig.ident', scale = F,dot.min = 1e-5,cols = c('#66C2A5','#FC8D62')) + 
    theme_bw()+
    scale_colour_gradientn(colours=c('#66C2A5','#FC8D62'),breaks=c(0,0.4,0.8),limits=c(0,0.8)
                           )
options(repr.plot.width = 6.5, repr.plot.height = 4.5)
p3
ggsave("WT_merge_replace_v2.SCT.reg6.RCTD.fcreg1.14celltype_MG.ident.Dotplot.png", width = 6.5, height = 4.5,dpi = 300)

#VH
da_wm<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="VH",])]
p3 <- DotPlot(da_wm, features = celltype_ad, group.by = 'orig.ident', scale = F,dot.min = 1e-5,cols = c('#66C2A5','#FC8D62')) + 
      theme_bw()+
      scale_colour_gradientn(colours=c('#66C2A5','#FC8D62'),breaks=c(0,0.4,0.8),limits=c(0,0.8)
                           )
options(repr.plot.width = 6.5, repr.plot.height = 4.5)
p3
ggsave("WT_merge_replace_v2.SCT.reg6.RCTD.fcreg1.14celltype_VH.ident.Dotplot.png", width = 6.5, height = 4.5,dpi = 300)

#DH
da_wm<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="DH",])]
p3 <- DotPlot(da_wm, features = celltype_ad, group.by = 'orig.ident', scale = F,dot.min = 1e-5,cols = c('#66C2A5','#FC8D62')) + 
      theme_bw()+
      scale_colour_gradientn(colours=c('#66C2A5','#FC8D62'),breaks=c(0,0.4,0.8),limits=c(0,0.8)
                           )
    #theme(axis.text.x = element_text(angle = 90, size = 12, margin = margin(10,0,0,0)),plot.title = element_text(hjust = 0.5,size = 15))+
    #ggtitle(label =  "Dorsal horn")
options(repr.plot.width = 6.5, repr.plot.height = 4.5)
p3
ggsave("WT_merge_replace_v2.SCT.reg6.RCTD.fcreg1.14celltype_DH.ident.Dotplot.png", width = 6.5, height = 4.5,dpi = 300)

###stacked area chart
df<-da@meta.data[,c("time",celltypes)]
df2<-df %>% gather(key = "cell_type",value = "total_proportion",-c(1))
df3<-df2 %>% group_by(time,cell_type) %>% summarise(n=sum(total_proportion)) %>% mutate(percentage=n/sum(n))
df3$time2<-NA
df3[grepl("WT_sham",df3$time,fixed = T),"time2"]<-1
df3[grepl("WT_3h",df3$time,fixed = T),"time2"]<-2
df3[grepl("WT_24h",df3$time,fixed = T),"time2"]<-3
df3[grepl("WT_72h",df3$time,fixed = T),"time2"]<-4
df3$cell_type<-factor(df3$cell_type,levels = celltypes)
options(repr.plot.width=10,repr.plot.height=6)
ggplot(df3,aes(x=time2,y=percentage,fill=cell_type))+
    geom_area(alpha=1,size=0.2,color="black")+
    scale_fill_manual(values =col
    )+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
         )
ggsave("WT.15celltypes.time.stacked.area.plot.png",width=10,height = 6,dpi=300)

###stack barplot of 72h_0.5mm
tmp<-meta_dec[meta_dec$time=='WT_72h',]
tmp2<-aggregate(tmp[,celltypes],list(tmp$pos),mean)
tmp3<-aggregate(tmp[,celltypes],list(tmp$sample),mean)
colnames(tmp2)[1]<-"distance"
colnames(tmp3)[1]<-"sample"
tmp2_l<-gather(tmp2,celltype,total_num,celltypes,factor_key = TRUE)
tmp2_l$distance<-factor(tmp2_l$distance,levels=c('WT_H_1mm','WT_H_0.5mm',
                                                 'WT_T_0.5mm','WT_T_1mm'))

tmp3_l<-gather(tmp3,celltype,total_num,celltypes,factor_key = TRUE)
tmp3_l$sample<-factor(tmp3_l$sample,
                        levels=c('WT_72h_H_R1_1mm_1','WT_72h_H_R1_1mm_2',
                                 'WT_72h_H_R1_1mm_3','WT_72h_H_R1_1mm_4',
                            'WT_72h_H_210323_1','WT_72h_H_210323_2',
                              'WT_72h_H_210323_3','WT_72h_H_210323_4',
                                 'WT_72h_T_R2_1','WT_72h_T_R2_2',
                                 'WT_72h_T_R2_3','WT_72h_T_R2_4',
                                 'WT_72h_T_R1_1mm_1','WT_72h_T_R1_1mm_2',
                                 'WT_72h_T_R1_1mm_3','WT_72h_T_R1_1mm_4'))
ggplot(tmp2_l) + 
  geom_bar(mapping = aes(x=distance, y=total_num,fill=celltype),
           position="fill", stat = "identity") + 
  scale_fill_manual(values=col) +
  labs(title = "celltype composition per distance at 72hpi 0.5mm", 
       x = "distance", 
       y = "celltype") + 
  theme_classic() + 

  #coord_flip() + 
  theme(axis.text = element_text(size = 8, hjust = 1))

ggplot(tmp3_l) + 
  geom_bar(mapping = aes(x=sample, y=total_num,fill=celltype),
           position="fill", stat = "identity") + 
  scale_fill_manual(values=col) +
  labs(title = "celltype composition per distance at 72hpi 0.5mm", 
       x = "sample", 
       y = "celltype") + 
  theme_classic() + 

  #coord_flip() + 
  theme(axis.text = element_text(size = 8, hjust = 1),
        axis.text.x = element_text(angle = 30))
