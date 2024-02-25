library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(clustree)
library(RCTD)
library(Matrix)
library(future)
library(grid)
library(RColorBrewer)
library(circlize)
library(igraph)
library(stringr)
library(cowplot)
library(biomaRt)
library(scales)
library(colorspace)
library(biomaRt)
library(clusterProfiler)
options(connectionObserver = NULL)
library(reshape2)
library(corrplot)
source("script/self_function/save_pheatmap_pdf.R")
source("script/self_function/cor_diagram.R")
###read and merge all sample data
sample_path <- "/data/h5_spatial_20210604/WT/"
f<-dir(sample_path)
f<-f[c(14,13,15,16,6,5,8,7,2,1,3,4,10,9,12,11)]
temp<-Seurat::Load10X_Spatial(paste0(sample_path,f[1]),slice = f[1])
da<-temp
Images(da)
da$orig.ident<-rep(f[1],dim(da@meta.data)[1])
da$barcode<-paste(da$orig.ident,rownames(da@meta.data),sep = "__")
da$x.axis<-da@images$WT_sham_H_R2_1mm@coordinates$row
da$y.axis<-da@images$WT_sham_H_R2_1mm@coordinates$col
da$imagerow<-da@images$WT_sham_H_R2_1mm@coordinates$imagerow
da$imagecol<-da@images$WT_sham_H_R2_1mm@coordinates$imagecol
da<-RenameCells(da,new.names = da$barcode)
head(da@meta.data)
for(i in 2:length(f)){
    temp<-Seurat::Load10X_Spatial(paste0(sample_path,f[i]),slice = f[i])
    temp$orig.ident<-rep(f[i],dim(temp@meta.data)[1])
    temp$barcode<-paste(temp$orig.ident,rownames(temp@meta.data),sep = "__")
    temp$x.axis<-temp@images[[f[i]]]@coordinates$row
    temp$y.axis<-temp@images[[f[i]]]@coordinates$col
    temp$imagerow<-temp@images[[f[i]]]@coordinates$imagerow
    temp$imagecol<-temp@images[[f[i]]]@coordinates$imagecol
    temp<-RenameCells(temp,new.names = temp$barcode)
    #head(da@meta.data)
    da<-merge(da,temp)
}
da
da$distance<-NA
da@meta.data[grepl("1mm",da@meta.data$orig.ident,fixed = T),]$distance<-"1mm"
da@meta.data[!grepl("1mm",da@meta.data$orig.ident,fixed = T),]$distance<-"0.5mm"
Images(da)
saveRDS(da,"WT.merge.replace_v2.raw.rds")

###split and define 4 repeat sections
meta<-da@meta.data
dim(meta)
for (i in unique(da$orig.ident)){  
  tiss <- meta[meta$orig.ident==i,]
  head(tiss)
  #tiss$V3<-as.numeric(tiss$V3)
  #tiss$V4<-as.numeric(tiss$V4)
  x.b<-max(tiss$x.axis)-min(tiss$x.axis)
  h<-hist(tiss$x.axis,breaks =x.b )
  n0<- which(h$counts==0)
  length(n0)
  
  if ( length(n0)==0) {
    print(paste0(i,"_x_not.divide"))
    #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_x_not.divide")
    #write.table(not,file =paste0(not,".txt") )
  } else {
    x.s<-mean(h$mids[n0])
    s24<-filter(tiss,tiss$x.axis < x.s )
    s13<-filter(tiss,tiss$x.axis > x.s)
    
    y.b13<-max(s13$y.axis)-min(s13$y.axis)
    h<-hist(s13$y.axis,breaks =y.b13)
    n0<- which(h$counts==0)
    
    if (length(n0)==0) {
      print(paste0(i,"_s13_y_not.divide"))
      #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_y_not.divide")
      #write.table(not,file =paste0(not,".txt") )
    } else {
      
      y.s13<-mean(h$mids[n0]) 
      
      s1<-filter(tiss, tiss$x.axis > x.s & tiss$y.axis > y.s13)
      s1$sample <-"1"
      
      s3<-filter(tiss, tiss$x.axis > x.s & tiss$y.axis < y.s13)
      s3$sample <-"3"      
      
      
    }
      
   y.b24<-max(s24$y.axis)-min(s24$y.axis)
   h<-hist(s24$y.axis,breaks =y.b24)
   n0<- which(h$counts==0)   
   if (length(n0)==0) {
      print(paste0(i,"_s24_y_not.divide"))
      #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_y_not.divide")
      #write.table(not,file =paste0(not,".txt") )
    } else {
      
      y.s24<-mean(h$mids[n0]) 
      
      s2<-filter(tiss, tiss$x.axis < x.s & tiss$y.axis > y.s24)
      s2$sample <-"2"
      
      s4<-filter(tiss, tiss$x.axis < x.s & tiss$y.axis < y.s24)
      s4$sample <-"4"      
      
      
    }}   
  
    ph.da <- rbind(s1,s2,s3,s4)
    write.csv(ph.da,file = paste0("meta/",i,".replicates_tissue_positon.csv"))
    write.csv(s1,file = paste0("meta/",i,"_meta.s1.csv" ))
    write.csv(s2,file = paste0("meta/",i,"_meta.s2.csv" ))
    write.csv(s3,file = paste0("meta/",i,"_meta.s3.csv" ))
    write.csv(s4,file = paste0("meta/",i,"_meta.s4.csv" ))
}
for (i in c("WT_sham_H_R2","WT_24h_H_R1_1mm")){  
  tiss <- meta[meta$orig.ident==i,]
  head(tiss)
  #tiss$V3<-as.numeric(tiss$V3)
  #tiss$V4<-as.numeric(tiss$V4)
  y.b<-max(tiss$y.axis)-min(tiss$y.axis)
  h<-hist(tiss$y.axis,breaks =y.b )
  n0<- which(h$counts==0)
  length(n0)
  
  if ( length(n0)==0) {
    print(paste0(i,"_y_not.divide"))
    #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_x_not.divide")
    #write.table(not,file =paste0(not,".txt") )
  } else {
    y.s<-mean(h$mids[n0])
    s34<-filter(tiss,tiss$y.axis < y.s )
    s12<-filter(tiss,tiss$y.axis > y.s)
    x.b12<-max(s12$x.axis)-min(s12$x.axis)
    h<-hist(s12$x.axis,breaks =x.b12)
    n0<- which(h$counts==0)
    
    if (length(n0)==0) {
      print(paste0(i,"_s12_x_not.divide"))
      #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_y_not.divide")
      #write.table(not,file =paste0(not,".txt") )
    } else {
      
      x.s12<-mean(h$mids[n0]) 
      
      s1<-filter(tiss, tiss$y.axis > y.s & tiss$x.axis > x.s12)
      s1$sample <-"1"
      
      s2<-filter(tiss, tiss$y.axis > y.s & tiss$x.axis < x.s12)
      s2$sample <-"2" }     
         
   x.b34<-max(s34$x.axis)-min(s34$x.axis)
   h<-hist(s34$x.axis,breaks =x.b34)
   n0<- which(h$counts==0)   
   if (length(n0)==0) {
      print(paste0(i,"_s34_x_not.divide"))
      #not <- paste0(gsub("_tissue_positon.csv","",files[i]),"_y_not.divide")
      #write.table(not,file =paste0(not,".txt") )
    } else {
      
      x.s34<-mean(h$mids[n0]) 
      
      s3<-filter(tiss, tiss$y.axis < y.s & tiss$x.axis > x.s34)
      s3$sample <-"3"
      
      s4<-filter(tiss, tiss$y.axis < y.s & tiss$x.axis < x.s34)
      s4$sample <-"4"      
      
      
    }   }
  
    ph.da <- rbind(s1,s2,s3,s4)
    write.csv(ph.da,file = paste0("meta/",i,".replicates_tissue_positon.csv"))
    write.csv(s1,file = paste0("meta/",i,"_meta.s1.csv" ))
    write.csv(s2,file = paste0("meta/",i,"_meta.s2.csv" ))
    write.csv(s3,file = paste0("meta/",i,"_meta.s3.csv" ))
    write.csv(s4,file = paste0("meta/",i,"_meta.s4.csv" ))
}

f<-dir("WT_replace_v2/meta",pattern = "replicate")
f
rep<-read.csv(paste0("WT_replace_v2/meta/",f[1]),header = T)
rownames(rep)<-rep[,1]
rep<-rep[,-1]
head(rep)

for(i in c(2:length(f))){
    tem<-read.csv(paste0("WT_replace_v2/meta/",f[i]),header = T)
    rownames(tem)<-tem[,1]
    tem<-tem[,-1]
    tem[1:5,1:5]
    rep<-rbind(rep,tem)
}
dim(rep)
table(rep$sample)

da@meta.data<-rep
head(da@meta.data)
da$sample<-paste(da$orig.ident,da$sample,sep = "_")
da$replicate<-sapply(as.character(da$sample),FUN = function(x) {substr(x,nchar(x),nchar(x))})
table(da$sample)
da$sample<-factor(da$sample,levels = c('WT_sham_H_R2_1mm_1','WT_sham_H_R2_1mm_2','WT_sham_H_R2_1mm_3','WT_sham_H_R2_1mm_4',
         'WT_sham_H_R2_1','WT_sham_H_R2_2','WT_sham_H_R2_3','WT_sham_H_R2_4',
         'WT_sham_T_210323_1','WT_sham_T_210323_2','WT_sham_T_210323_3','WT_sham_T_210323_4',
         'WT_sham_T_R2_1mm_1','WT_sham_T_R2_1mm_2','WT_sham_T_R2_1mm_3','WT_sham_T_R2_1mm_4',
         'WT_3h_H_R2_1mm_1','WT_3h_H_R2_1mm_2','WT_3h_H_R2_1mm_3','WT_3h_H_R2_1mm_4',
         'WT_3h_H_R2_1','WT_3h_H_R2_2','WT_3h_H_R2_3','WT_3h_H_R2_4',
         'WT_3h_T_R2_1','WT_3h_T_R2_2','WT_3h_T_R2_3','WT_3h_T_R2_4',
         'WT_3h_T_210330_1mm_1','WT_3h_T_210330_1mm_2','WT_3h_T_210330_1mm_3','WT_3h_T_210330_1mm_4',
         'WT_24h_H_R1_1mm_1','WT_24h_H_R1_1mm_2','WT_24h_H_R1_1mm_3','WT_24h_H_R1_1mm_4',
         'WT_24h_H_R1_1','WT_24h_H_R1_2','WT_24h_H_R1_3','WT_24h_H_R1_4',
         'WT_24h_T_201231_1','WT_24h_T_201231_2','WT_24h_T_201231_3','WT_24h_T_201231_4',
         'WT_24h_T_R1_1mm_1','WT_24h_T_R1_1mm_2','WT_24h_T_R1_1mm_3','WT_24h_T_R1_1mm_4',
         'WT_72h_H_R1_1mm_1','WT_72h_H_R1_1mm_2','WT_72h_H_R1_1mm_3','WT_72h_H_R1_1mm_4',
         'WT_72h_H_210323_1','WT_72h_H_210323_2','WT_72h_H_210323_3','WT_72h_H_210323_4',
         'WT_72h_T_R2_1','WT_72h_T_R2_2','WT_72h_T_R2_3','WT_72h_T_R2_4',
         'WT_72h_T_R1_1mm_1','WT_72h_T_R1_1mm_2','WT_72h_T_R1_1mm_3','WT_72h_T_R1_1mm_4'))

da$orig.ident<-factor(temp$orig.ident,levels = c('WT_sham_H_R2_1mm','WT_sham_H_R2','WT_sham_T_210323','WT_sham_T_R2_1mm',
                                                  'WT_3h_H_R2_1mm','WT_3h_H_R2','WT_3h_T_R2','WT_3h_T_210330_1mm',
                                                  'WT_24h_H_R1_1mm','WT_24h_H_R1','WT_24h_T_201231','WT_24h_T_R1_1mm',
                                                  'WT_72h_H_R1_1mm','WT_72h_H_210323','WT_72h_T_R2','WT_72h_T_R1_1mm'))
saveRDS(da,"WT.merge.replace_v2.raw.rds")


###Quality Control and visualization

#remove blood related genes
Hb_gene<-rownames(da)[grepl("^Hb",rownames(da))]
Hb_gene
da<-da[!rownames(da)%in%Hb_gene,colnames(da)]
da
#m.s.genes<-convertHumanGeneList(cc$s.genes)
m.s.genes<-c('Exo1','Msh2','Mcm4','Chaf1b','Rrm2','Cenpu','Mrpl36','Gmnn','Hells','Cdc6','Gins2','Uhrf1',
             'Cdc45','Slbp','Ubr7','Fen1','Rad51ap1','Mcm5','Nasp','Cdca7','Blm','Usp1','Ung','Prim1',
             'Clspn','Mcm6','Dtl','Pola1','Dscc1','Tipin','Wdr76','Casp8ap2','Tyms','Ccne2','Rrm1','Polr1b',
             'Rfc2','Rad51','E2f8','Mcm7','Pcna')
#m.g2m.genes<-convertHumanGeneList(cc$g2m.genes)
m.g2m.genes<-c('Dlgap5','Ctcf','Smc4','Kif20b','Cdc25c','Gtse1','Tpx2','Hmgb2','Cks1brt','Cdca2',
               'Top2a','Cks2','Cdca3','G2e3','Ttk','Ncapd2','Lbr','Anp32e','Ckap2','Tacc3','Nek2',
               'Cenpe','Kif11','Anln','Hjurp','Aurkb','Rangap1','Cks1b','Hmmr','Ckap5','Cdc20',
               'Psrc1','Kif23','Ect2','Kif2c','Ndc80','Nuf2','Cdca8','Birc5','Cenpf','Ube2c','Jpt1',
               'Pimreg','Nusap1','Mki67','Tubb4b','Bub1','Cenpa','Ccnb2','Aurka','Ckap2l')
da<-PercentageFeatureSet(da,"^mt-",col.name = "percent_mt")
da<-NormalizeData(da)
da<-CellCycleScoring(da,s.features = m.s.genes,g2m.features = m.g2m.genes)
da<-FindVariableFeatures(da,nfeatures = 3000,selection.method ="vst")
da<-ScaleData(da,vars.to.regress = c("S.Score","G2M.Score","nCount_Spatial","percent_mt","orig.ident"),features = rownames(da))
saveRDS(da,"WT.merge.replace_v2.regress_CC.nC.mt.ident.rds")

slide<-unique(da$orig.ident)
slide
gene_detect<-data.frame(row.names = slide)
##dir.create("QC/")
for(i in slide){
    cells<-rownames(da@meta.data)[da$orig.ident==i]
    temp<-da@assays$Spatial@counts[,cells]
    gene_detect[i,"Total gene detected"]<-sum(rowSums(temp)>0)
    gene_detect[i,"Num of Spots"]<-length(cells)
}
gene_detect

stats<-data.frame(row.names = slide)
for(i in slide){
    stats[i,paste0("Mean of ",c('nCount','nFeature','percent_mt','S.Score','G2M.Score'))]<-apply(da@meta.data[da@meta.data$orig.ident==i,c('nCount_Spatial','nFeature_Spatial','percent_mt','S.Score','G2M.Score')],2,mean)
    stats[i,paste0("Median of ",c('nCount','nFeature','percent_mt','S.Score','G2M.Score'))]<-apply(da@meta.data[da@meta.data$orig.ident==i,c('nCount_Spatial','nFeature_Spatial','percent_mt','S.Score','G2M.Score')],2,median)               
}
#stats[slide[1],c('nCount_Spatial','nFeature_Spatial','percent_mt','S.Score','G2M.Score')]<-colSums(da@meta.data[da@meta.data$orig.ident==slide[1],c('nCount_Spatial','nFeature_Spatial','percent_mt','S.Score','G2M.Score')])
stats
#write.csv(stats,"QC/WT.merge.QC_average.csv")   
stats2<-cbind(stats,gene_detect)
stats2
colMeans(stats)
write.csv(stats2,"QC/WT.merge.16sample.QC_220920.csv")

options(repr.plot.width=10,repr.plot.height=6)
VlnPlot(da,cols=colorRampPalette(c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF'))(16),
            features = "nCount_Spatial",
            group.by = "orig.ident",
            pt.size = 0)+
    NoLegend()+
    theme(axis.text.x = element_text(angle = 90))
ggsave("WT.merge.Spatial_nCount.vlnplot_220801.png",width = 10,height = 6,dpi = 300)

options(repr.plot.width=10,repr.plot.height=6)
VlnPlot(da,cols=colorRampPalette(c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF'))(16),
            features = "nFeature_Spatial",
            group.by = "orig.ident",
            pt.size = 0)+
    NoLegend()+
    theme(axis.text.x = element_text(angle = 90))
ggsave("WT.merge.Spatial_nFeature.vlnplot_220801.png",width = 10,height = 6,dpi = 300)

options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="nCount_Spatial"))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11,name = "Spectral")))(50))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
          )
ggsave("WT.merge.Spatial_nCount.spatial_220801.png",width = 20,height = 16,dpi = 300)

options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="nFeature_Spatial"))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11,name = "Spectral")))(50))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
          )
ggsave("WT.merge.Spatial_nFeature.spatial_220801.png",width = 20,height = 16,dpi = 300)


###add meta info
da@meta.data$time<-NA
da@meta.data[grepl("sham",da@meta.data$orig.ident,fixed = T),]$time<-"WT_sham"
da@meta.data[grepl("3h",da@meta.data$orig.ident,fixed = T),]$time<-"WT_3h"
da@meta.data[grepl("24h",da@meta.data$orig.ident,fixed = T),]$time<-"WT_24h"
da@meta.data[grepl("72h",da@meta.data$orig.ident,fixed = T),]$time<-"WT_72h"
da$time<-factor(da$time,levels = c('WT_sham','WT_3h','WT_24h','WT_72h'))
unique(da$time)

da@meta.data$RC<-NA
da@meta.data[grepl("_H_",da@meta.data$orig.ident,fixed = T),]$RC<-"H"
da@meta.data[grepl("_T_",da@meta.data$orig.ident,fixed = T),]$RC<-"T"
unique(da$RC)

da$distance<-"0.5mm"
da@meta.data[grepl("1mm",da$orig.ident,fixed = T),]$distance<-"1mm"

da$pos<-paste(da$RC,da$distance,sep="_")
da$pos<-factor(da$pos,levels = c('H_1mm','H_0.5mm','T_0.5mm','T_1mm'))
unique(da$pos)
##manually adjust spatial coordinates
sample_coor<-da@meta.data[da@meta.data$orig.ident==Images(da)[1],c("x.axis","y.axis")]
#head(sample_coor)
range(sample_coor$x.axis)
range(sample_coor$y.axis)
meta<-da@meta.data
range(meta$x.axis)
range(meta$y.axis)
meta$x.ad<-NA
meta$y.ad<-NA
for(i in 1:length(Images(da))){
    m=(i-1)%/%4
    x.add<-(i-4*m-1)*100
    y.add<-150*(4-m)
    meta[meta$orig.ident==Images(da)[i],"x.ad"]<-meta[meta$orig.ident==Images(da)[i],"x.axis"]+x.add
    meta[meta$orig.ident==Images(da)[i],"y.ad"]<-meta[meta$orig.ident==Images(da)[i],"y.axis"]+y.add

}
range(meta$x.ad)
range(meta$y.ad)
da@meta.data<-meta

###manually adjust the coordinates of representative samples
s<-c('WT_sham_H_R2_1mm_2','WT_sham_H_R2_4','WT_sham_T_210323_4','WT_sham_T_R2_1mm_4',
     'WT_3h_H_R2_1mm_2','WT_3h_H_R2_2','WT_3h_T_R2_1','WT_3h_T_210330_1mm_3',
     'WT_24h_H_R1_1mm_3','WT_24h_H_R1_4','WT_24h_T_201231_3','WT_24h_T_R1_1mm_3',
     'WT_72h_H_R1_1mm_2','WT_72h_H_210323_1','WT_72h_T_R2_4','WT_72h_T_R1_1mm_3')
da_s<-da[rownames(da@assays$SCT@counts),rownames(da@meta.data)[da$sample %in% s]]
da_s
me<-da_s@meta.data
me$x.ad<-0
me$y.ad<-0
for(i in s){
    temp<-me[me$sample==i,c("x.axis","y.axis")]
    xmin<-min(temp$x.axis)
    ymin<-min(temp$y.axis)
    temp$x.ad<-temp$x.axis-xmin
    temp$y.ad<-temp$y.axis-ymin
    me[me$sample==i,"x.ad"]<-temp$x.ad
    me[me$sample==i,"y.ad"]<-temp$y.ad
}
me$x.ad2<-0
me$y.ad2<-0
for(i in 1:length(s)){
    m=(i-1)%/%4
    x.add<-(i-4*m-1)*40
    y.add<-50*(3-m)
    me[me$sample==s[i],"x.ad2"]<-me[me$sample==s[i],"x.ad"]+x.add
    me[me$sample==s[i],"y.ad2"]<-me[me$sample==s[i],"y.ad"]+y.add

}
range(me$x.ad2)
range(me$y.ad2)
options(repr.plot.width=10,repr.plot.height=7)
ggplot(me,aes(x.ad2,y.ad2,group=SCT_snn_res.0.2))+geom_point(size=1.6,shape=21,stroke=0.15,aes(fill=SCT_snn_res.0.2))+
    #scale_fill_manual(values = col)+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
    )
write.csv(me,"16samples.meta.ad.csv")

###SCTransform
da<-SCTransform(da,assay = "Spatial",
                vars.to.regress = c("S.Score","G2M.Score","nCount_Spatial","percent_mt","orig.ident"),
                return.only.var.genes = F,
               do.scale=T)
da

##sample correlation based on bulk expression
mean_expr<-t(as.matrix(da@assays$SCT@counts))
mean_expr<-as.data.frame(mean_expr)
mean_expr$sample<-da$sample
mean_expr[1:5,1:5]
mean_expr %>% group_by(sample) %>%
    summarise_all("mean") -> df
df<-as.data.frame(df)
rownames(df)<-df$sample
df<-df[,-which(colnames(df)=="sample")]
dim(df)
df[1:5,1:5]
write.csv(df,"WT.replace_v2.SCT.sample.average.expression.csv")
pdf("64samples.SCTcount.average.expression.pearson.correlation.corrplot.pdf",width=40,height=40)
options(repr.plot.width=30,repr.plot.height=30)
corrplot(cor,method="color",is.corr = FALSE,
         type="upper",
         tl.col = "black",col.lim = c(0.65,1),
         col = rev(COL2("RdBu")))#,
               #lower.col = "#92a8d1",
               #col = colorRampPalette(colors = rev(col)))
dev.off()

#injury score
injury_genes<-c('Adamts1','Atf3','Ccl2','Ccnd1','Cd68','Cebpd','Cyba',
                'Fn1','Gal','Gap43','Hmox1','Hspb1','Igfbp2','Jun','Junb',
                'Fos','Lgals1','Neat1','Socs3','Tnc','S100a10','Timp1')
da<-AddModuleScore(da,features = list(injury_genes),name = "injury_score")
options(repr.plot.width=10,repr.plot.height=6)
VlnPlot(da,cols=colorRampPalette(c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF'))(16),
            features = "injury_score1",
            group.by = "orig.ident",
            pt.size = 0)+
    NoLegend()+
    theme(axis.text.x = element_text(angle = 90))
ggsave("WT.merge.SCT_22g.injuryScore.vlnplot_220801.png",width = 10,height = 6,dpi = 300)

options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="injury_score1"))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11,name = "Spectral")))(50))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
          )
ggsave("WT.merge.SCT_22g.injuryScore.spatial_220801.png",width = 20,height = 16,dpi = 300)


###extract SCT.data (gene expressed >0.1% spots) to pyscenic
SCTdata<-da@assays$SCT@data
head(apply(da@assays$SCT@counts,1,function(x) sum(x>0)))
colnames(da@assays$Spatial@meta.features)
da@assays$SCT@meta.features$expressed_spots<-apply(da@assays$SCT@counts,1,function(x) sum(x>0))
table(da@assays$SCT@meta.features$expressed_spots>20)  
scenic_ma<-da@assays$SCT@data[rownames(da)[da@assays$SCT@meta.features$expressed_spots>20],colnames(da)]
dim(scenic_ma)
write.csv(as.matrix(scenic_ma),"WT.merge.replace_v2.express_cells20.18447g.22820allspots.SCTnormalizedData.csv")


###Clustering
da<-RunPCA(da,features = VariableFeatures(da))
da<-FindNeighbors(da,reduction = "pca",dims = 1:20,k.param = 50)
for(i in c(0.1,0.15,0.2,0.25,0.3)){
    da<-FindClusters(da,resolution = i,verbose = F)
}

da<-RunUMAP(da,reduction = "pca",dims = 1:20)
da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")

col1<-c('#1C8356','#782AB6','#3283FE','#F7E1A0')
options(repr.plot.width=9,repr.plot.height=7)
DimPlot(da,reduction="umap",label=F,cols = c('#6b5b95','#eca1a6','#c4b7a6','#92a8d1'
),#alpha(c("#1CE6FF","#e6afb9","#A30059","#023fa5"),0.66                                                
        group.by = c("replicate"))+
theme(#legend.position ="bottom",
    axis.ticks = element_blank(),axis.text = element_blank(),
         axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())
ggsave("WT.merge.SCT.regress_CC.nC.mt.ident_pc20.k50.replicate.UMAP.png",width = 9,height = 7,dpi = 400)

options(repr.plot.width=9,repr.plot.height=7)
DimPlot(da,reduction="umap",label=F,cols = c('#F4A582', '#0571B0'),#alpha(c('#F4A582', '#0571B0'),0.66),
        group.by = c("RC"))+
theme(#legend.position ="bottom",
    axis.ticks = element_blank(),axis.text = element_blank(),
         axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())

ggsave("WT.merge.SCT.regress_CC.nC.mt.ident_pc20.k50.RC.UMAP.png",width = 9,height = 7,dpi = 400)

options(repr.plot.width=9,repr.plot.height=7)
DimPlot(da,reduction="umap",label=F,cols = c('#63A79C','#FB9A99'),#alpha(c('#63A79C','#DB95BB'),0.66),
        group.by = c("distance"))+
theme(#legend.position ="bottom",
    axis.ticks = element_blank(),axis.text = element_blank(),
         axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())

ggsave("WT.merge.SCT.regress_CC.nC.mt.ident_pc20.k50.distance.UMAP.png",width = 9,height = 7,dpi = 400)

##assign clusters to anatomic domains
da$domain_res02<-NA
da@meta.data[da@meta.data$SCT_snn_res.0.2 %in% c("0","3","4"),"domain_res02"]<-"WM"
da@meta.data[da@meta.data$SCT_snn_res.0.2 %in% c("1","2"),"domain_res02"]<-"MG"
da@meta.data[da@meta.data$SCT_snn_res.0.2 %in% c("5"),"domain_res02"]<-"DH"
da@meta.data[da@meta.data$SCT_snn_res.0.2 %in% c("6"),"domain_res02"]<-"VH"
unique(da$domain_res02)
table(da$domain_res02)
col_domain<-c('#BC3C29A8','#0072B5A8','#E18727A8','#20854EA8')
names(col_domain)<-c("DH","MG","VH","WM")
options(repr.plot.width=9,repr.plot.height=7)
DimPlot(da,reduction="umap",label=F,cols = col_domain,#alpha(c("#f79cd4","#1CE6FF","#FFB500","#A4E804"),0.66),
        group.by = c("domain_res02"))+
theme(#legend.position ="bottom",
    axis.ticks = element_blank(),axis.text = element_blank(),
         axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())

ggsave("WT.merge.SCT.regress_CC.nC.mt.ident_pc20.k50.domain_res02.UMAP.png",width = 9,height = 7,dpi = 400)

col<-c('#66AB94','#5D3B83','#316B9B','#E2E293','#C580AA','#C13031','#E18727A8')
show_col(col)
col
names(col)<-c("0","1","2","3","4","5","6")
da$SCT_snn_res.0.2<-factor(da$SCT_snn_res.0.2,levels = c("0","1","2","3","4","5","6"))
Idents(da)<-da$SCT_snn_res.0.2
options(repr.plot.width=9,repr.plot.height=7)
DimPlot(da,reduction="umap",label=F,cols = col,
        group.by = c("SCT_snn_res.0.2"))+
theme(#legend.position ="bottom",
    axis.ticks = element_blank(),axis.text = element_blank(),
         axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())

ggsave("WT.merge.regress_CC.nC.mt.ident_SCT.pc20.k50.res02.UMAP.png",width = 9,height = 7,dpi = 400)

pdf("WT.merge.regress_CC.nC.mt.ident_SCT.pc20.k50.res02.spatial_manually.pdf",width=20,height=16)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(meta,aes(x.ad,y.ad,group=SCT_snn_res.0.2))+geom_point(size=1.6,shape=21,stroke=0.15,aes(fill=SCT_snn_res.0.2))+
    scale_fill_manual(values = col)+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
    )
dev.off()
write.csv(da@meta.data,"WT.SCT.pc20.k50.res02.meta.data.csv")

#barplot of cluster per variables
meta<-da@meta.data
meta$n<-1
meta$domain_res02<-as.character(meta$domain_res02)
meta$domain_rep<-paste0(meta$domain_res02,"_",meta$replicate)
meta$time_rep<-paste0(meta$time,"_",meta$replicate)
meta$RC_rep<-paste0(meta$RC,"_",meta$replicate)
meta$dis_rep<-paste0(meta$distance,"_",meta$replicate)
domain_cluster<-aggregate(meta$n,list(meta$domain_res02,meta$SCT_snn_res.0.2),sum)
colnames(domain_cluster)<-c("domain","cluster","spot")
domain_rep_cluster<-aggregate(meta$n,list(meta$domain_rep,meta$SCT_snn_res.0.2),sum)
colnames(domain_rep_cluster)<-c("domain","cluster","spot")
domain_cluster$cluster<-as.character(domain_cluster$cluster)
domain_rep_cluster$cluster<-as.character(domain_rep_cluster$cluster)
time_cluster<-aggregate(meta$n,list(meta$time,meta$SCT_snn_res.0.2),sum)
colnames(time_cluster)<-c("time","cluster","spot")
time_rep_cluster<-aggregate(meta$n,list(meta$time_rep,meta$SCT_snn_res.0.2),sum)
colnames(time_rep_cluster)<-c("time","cluster","spot")
time_cluster$cluster<-as.character(time_cluster$cluster)
time_rep_cluster$cluster<-as.character(time_rep_cluster$cluster)
rc_cluster<-aggregate(meta$n,list(meta$RC,meta$SCT_snn_res.0.2),sum)
colnames(rc_cluster)<-c("RC","cluster","spot")
rc_rep_cluster<-aggregate(meta$n,list(meta$RC_rep,meta$SCT_snn_res.0.2),sum)
colnames(rc_rep_cluster)<-c("RC","cluster","spot")
rc_cluster$cluster<-as.character(rc_cluster$cluster)
rc_rep_cluster$cluster<-as.character(rc_rep_cluster$cluster)
dis_cluster<-aggregate(meta$n,list(meta$distance,meta$SCT_snn_res.0.2),sum)
colnames(dis_cluster)<-c("distance","cluster","spot")
dis_rep_cluster<-aggregate(meta$n,list(meta$dis_rep,meta$SCT_snn_res.0.2),sum)
colnames(dis_rep_cluster)<-c("distance","cluster","spot")
dis_cluster$cluster<-as.character(dis_cluster$cluster)
dis_rep_cluster$cluster<-as.character(dis_rep_cluster$cluster)

options(repr.plot.width=6,repr.plot.height=3.5)
pdf("240104_stackbar_domainPerCluster.pdf",width = 6,height = 3.5)
ggplot(domain_cluster) + 
  geom_bar(mapping = aes(x=domain, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "domain per cluster", 
       x = "domain", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=5)
pdf("240104_stackbar_domainPerCluster_replicates.pdf",width = 6,height = 5)
ggplot(domain_rep_cluster) + 
  geom_bar(mapping = aes(x=domain, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "domain per cluster", 
       x = "domain", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=3.5)
pdf("240104_stackbar_clusterPerTime.pdf",width = 6,height = 3.5)
ggplot(time_cluster) + 
  geom_bar(mapping = aes(x=time, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per time", 
       x = "cluster", 
       y = "time") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=5)
pdf("240104_stackbar_clusterPerTime_replicates.pdf",width = 6,height = 5)
ggplot(time_rep_cluster) + 
  geom_bar(mapping = aes(x=time, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per time", 
       x = "cluster", 
       y = "time") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=2.5)
pdf("240104_stackbar_clusterPerDirection.pdf",width = 6,height = 2.5)
ggplot(rc_cluster) + 
  geom_bar(mapping = aes(x=RC, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per RC", 
       x = "RC", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=5)
pdf("240104_stackbar_clusterPerDirection_replicate.pdf",width = 6,height = 5)
ggplot(rc_rep_cluster) + 
  geom_bar(mapping = aes(x=RC, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per RC", 
       x = "RC", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=2.5)
pdf("240104_stackbar_clusterPerDistance.pdf",width = 6,height = 2.5)
ggplot(dis_cluster) + 
  geom_bar(mapping = aes(x=distance, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per distance", 
       x = "distance", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

options(repr.plot.width=6,repr.plot.height=5)
pdf("240104_stackbar_clusterPerDistance_replicate.pdf",width = 6,height = 5)
ggplot(dis_rep_cluster) + 
  geom_bar(mapping = aes(x=distance, y=spot,fill=cluster), position="stack", stat="identity") + 
  scale_fill_manual(values = col)+
  labs(title = "cluster per distance", 
       x = "distance", 
       y = "cluster") + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text = element_text(size = 10, hjust = 1))
dev.off()

###cluster composition each time-domain
df<-da@meta.data[,c("time","domain_res02","SCT_snn_res.0.2")]
head(df)
df$group<-paste0(df$time,"_",df$domain_res02)
df<-df[,c(1，3,4)]
colnames(df)<-c("time","cluster","group")
df$value<-1
head(df)
temp<-df %>% group_by(group,cluster) %>% summarise(spot_num=sum(value))
head(temp)
temp$group<-factor(temp$group,levels = c(
                                         #"g1",
                                         
                                        # "g2",
                                         'WT_sham_DH',
                                         'WT_3h_DH',
                                         'WT_24h_DH',
                                         'WT_72h_DH',
                                         #"g3",
                                         'WT_sham_VH',
                                        'WT_3h_VH',
                                        'WT_24h_VH',
                                        'WT_72h_VH',
                                         'WT_sham_MG',
                                         'WT_3h_MG',
                                         'WT_24h_MG',
                                        'WT_72h_MG',
                                        'WT_sham_WM',
                                         'WT_3h_WM',
                                         'WT_24h_WM',
                                         'WT_72h_WM'
                                        )
                   #c('WT_sham_WM','WT_sham_MG','WT_sham_DH','WT_sham_VH',
                   #                     'WT_3h_WM','WT_3h_MG','WT_3h_DH','WT_3h_VH',
                   #                     'WT_24h_WM','WT_24h_MG','WT_24h_DH','WT_24h_VH',
                   #                     'WT_72h_WM','WT_72h_MG','WT_72h_DH','WT_72h_VH')
                  )
options(repr.plot.width=7,repr.plot.height=6)
ggplot(temp,aes(x=group,y=spot_num,fill=cluster))+
    geom_bar(position = "stack",stat="identity")+
    scale_fill_manual(values = c('#66AB94','#5D3B83','#316B9B','#E2E293','#C580AA','#C13031','#E18727A8',"white"))+
    theme_bw()+
    theme(panel.grid  =element_blank(),axis.text.x.bottom = element_text(angle = 90),axis.line.x = element_line(colour = "black"))#+
    #coord_cartesian(ylim = c(0,100))
ggsave("domain_cluster.barplot.png",width =7,height = 6,dpi = 300)


