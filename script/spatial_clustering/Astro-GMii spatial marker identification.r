library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
df<-read.csv("/home/jovyan/zxli_SCI/result/Seurat/reg.CC/WT_replace_v2/res02_220310/Astrocyte/GM_Gfap.intersect.astr_marker.GBA.MG_gfap.long.df.csv")
head(gba)
GMii<-read.csv("/home/jovyan/zxli_SCI/result/Seurat/sc_Astrocyte_recluster/lee/classifyBySpatialMarker/gfap.317marker.k4.m13.hclust2.fc1.pct05.pad001.pos.markers.csv")
GMii<-GMii$gene
length(GMii)
da<-AddModuleScore(da,features = list(GMii),name = "GMii_score")
options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="GMii_score1"))+
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
ggsave("WT.GMii_Score.spatial_220916.png",width = 20,height = 16,dpi = 400)

da$Gfap<-da@assays$SCT@data["Gfap",]
da$Pla2g7<-da@assays$SCT@data["Pla2g7",]
da$Slc7a10<-da@assays$SCT@data["Slc7a10",]
pdf("WT.Pla2g7.spatial_230803.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="Pla2g7"))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11,name = "Spectral")))(50))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()
dev.off()
pdf("Slc7a10.SCT.data.spatial.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="Slc7a10"))+
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
dev.off()
alldeg<-read.csv("lee/classifyBySpatialMarker/3state.alldeg.5groups.csv")
alldeg<-alldeg[,-1]
head(alldeg)
genes<-alldeg[grep("GM_Gfap",alldeg$group),"gene"]
genes

genes<-intersect(genes,rownames(da))
temp<-as.matrix(da@assays$SCT@data[genes,])
head(temp)
EM<-apply(temp,1,function(x) (summary(x)[[3]]/2+(summary(x)[[2]]+summary(x)[[5]])/4))
test_bi<-temp
for(i in rownames(temp)){
    temp_row<-temp[i,]
    temp_EM<-EM[[i]]
    temp_row[temp_row>=temp_EM]<-1
    temp_row[temp_row<temp_EM]<-0
    test_bi[i,]<-temp_row
}
ct_df<-read.csv("RCTD/testC5/WT_all/fcreg1_220503/16celltypes.dec_conf.nor.meta.csv")
rownames(ct_df)<-ct_df[,1]
ct_df<-ct_df[,-1]
head(ct_df)
ct<-colnames(ct_df)[46:61]
ct
wt_df<-ct_df[,ct]
dim(wt_df)
wt_df_bi<-wt_df
wt_df_bi<-apply(wt_df_bi,2,function(x) (ifelse(x>0,1,0)))
head(wt_df_bi)
me<-cbind(da@meta.data,wt_df_bi)
dim(me)

me$Astr_Gfap<-paste0(me$domain_res02,"_",me$Astrocyte.Gfap)
unique(me$Astr_Gfap)
me$Astr_Gfap<-ifelse(me$Astr_Gfap=="MG_1",1,0)
table(me$Astr_Gfap)

####Fisher's exact test
    #calculate LR's association with cellpair and get the LR-cp double sig spot number
    lr_cp_pval_df<-data.frame(row.names = rownames(test_bi))
    lr_cp_spot_n<-data.frame(row.names = rownames(test_bi))
    #for(i in colnames(astr_wt_bi)){
        lr_cp_pval_df[,"MG_1"]<-NA
        lr_cp_spot_n[,"MG_1"]<-NA
        for(j in rownames(test_bi)){
            #if(all(table(lr_bi_df[,j],cp_bi_df_sub[,i])==2)){
            pp1<-test_bi[j,]
            pp1<-factor(pp1,levels=c(0,1))
            pp2<-me[,"Astr_Gfap"]
            pp2<-factor(pp2,levels=c(0,1))
            f_df<-table(pp1,pp2)
            p<-fisher.test(f_df,alternative = "greater")$p.val
            lr_cp_pval_df[j,"MG_1"]<-p
            
            n<-f_df[2,2]
            lr_cp_spot_n[j,"MG_1"]<-n
            #}
            
        }
   # }
lr_cp_n_long<-lr_cp_spot_n
lr_cp_n_long$gene<-rownames(lr_cp_n_long)
lr_cp_n_long<-melt(lr_cp_n_long)
colnames(lr_cp_n_long)<-c("gene","celltype","n_spot")
lr_cp_n_long$key<-paste0(lr_cp_n_long$gene,"-",lr_cp_n_long$celltype)
head(lr_cp_n_long)
lr_cp_long_df<-lr_cp_pval_df
lr_cp_long_df$gene<-rownames(lr_cp_long_df)
lr_cp_long_df<-melt(lr_cp_long_df)
lr_cp_long_df$key<-paste0(lr_cp_long_df$gene,"-",lr_cp_long_df$variable)
head(lr_cp_long_df)
lr_cp_long<-merge(lr_cp_long_df,lr_cp_n_long[,c(3,4)],by = "key")
head(lr_cp_long)
dim(lr_cp_long)
write.csv(lr_cp_long,"Astrocyte/GM_Gfap.deg.GBA.MG_gfap.long.df.csv")

genes<-unique(lr_cp_long[lr_cp_long$value<0.05 & lr_cp_long$n_spot>1000,"gene"])
meta<-da@meta.data
meta[,genes]<-da@assays$SCT@data[genes,]
options(repr.plot.width=20,repr.plot.height=16)
for(i in genes){
    ggplot(meta,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill=i))+
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
ggsave(paste0("lee/GM_Gfap/",i,".GBA.spatial_manually.png"),width=20,height=16,dpi=400)
}


