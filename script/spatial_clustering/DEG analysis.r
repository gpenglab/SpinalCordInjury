library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(future)
library(grid)
library(RColorBrewer)
library(stringr)
library(biomaRt)
library(scales)
library(colorspace)
library(clusterProfiler)
options(connectionObserver = NULL)
library(reshape2)

da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
da

###domains marker display including well-known and novel in our sham data
da_sham<-da[rownames(da),rownames(da@meta.data[da@meta.data$time=="WT_sham",])]
da_sham
unique(da_sham$time)
domain_marker<-FindAllMarkers(da_sham,logfc.threshold = 0.5,min.pct = 0.5,only.pos = T)
domain_marker<-domain_marker[domain_marker$p_val_adj<0.01,]
write.csv(domain_marker,"sham_domain.marker.fc05.pct05.p001.pos.csv")
domain_marker %>% 
    group_by(cluster) %>%
    #filter(pct.1-pct.2>0.15) %>%
    top_n(n=20,wt = avg_log2FC) ->top10
top10
sf<-"domain_marker/"
if(!dir.exists(sf))
    dir.create(sf)
options(repr.plot.width=9,repr.plot.height=7)
for(i in unique(top10$gene)){
    p<-FeaturePlot(da,reduction="umap",features = i,slot = "scale.data",min.cutoff = median(da@assays$SCT@scale.data[i,]))+
        theme(#legend.position ="bottom",
            axis.ticks = element_blank(),axis.text = element_blank(),
            axis.line = element_line(size = 0.8,arrow = arrow()) )#+scale_x_discrete(breaks=c())
    plot(p)
    ggsave(paste0(sf,"WT.merge.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.sham_domain_marker.scale.data.cutMedian.",i,".UMAP_220614.png"),width = 9,height = 7,dpi = 400)
}
options(repr.plot.width=20,repr.plot.height=18)
for(i in unique(top10$gene)){
    p<-SpatialFeaturePlot(da,features = i,slot = "scale.data",combine = F,min.cutoff = median(da@assays$SCT@scale.data[i,])
                  )

    do.call(ggpubr::ggarrange,c(p,list(ncol=4,nrow=4,common.legend=T,legend="right")))
    ggsave(paste0(sf,"WT.merge.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.sham_domain_marker.scale.data.cutMedian.",i,".spatial_220614.png"),
           width = 20,height = 18,dpi = 400)
}

###cluster marker spatial display
Idents(da)<-da$SCT_snn_res.0.2
cluster_marker<-FindAllMarkers(da,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.25 ,only.pos = T)
write.csv(cluster_marker,"pc20_k50.res02.fc025.pct025.pos.clusterMarker.csv")
cluster_marker %>% 
    group_by(cluster) %>%
    #filter(pct.1-pct.2>0.15) %>%
    top_n(n=10,wt = avg_log2FC) ->top10
top10
me<-da@meta.data
temp<-t(as.matrix(da@assays$SCT@data[unique(top10$gene),]))
me<-cbind(me,temp)
sf<-"cluster_marker/"
if(!dir.exists(sf))
    dir.create(sf)
options(repr.plot.width=20,repr.plot.height=16)
for(i in unique(top10$gene)){
    ggplot(me,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill=i))+
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
ggsave(paste0(sf,i,".spatial_manually.png"),width=20,height=16,dpi=400)
}

#cluster markers' functional pathway
temp<- cluster_marker %>% filter(p_val_adj<0.01 )
backid<-bitr(rownames(da),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
backgenes<-backid[,2]
write.table(backgenes,"WT.GOterm.SCT.backgeneids.txt",row.names = T,col.names = T)
library(org.Mm.eg.db)
library(ggplot2)
grlabs<-split(cluster_marker
cluster)
gcSample = lapply(grlabs, function(gr) as.numeric(bitr(gr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID))
pvalueCutoff = 0.01
qvalueCutoff = 0.01
xx.mus.go <- compareCluster(gcSample, OrgDb='org.Mm.eg.db', fun='enrichGO', 
                            pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ont = "BP", 
                            readable=T,universe=backgenes
                           )
saveRDS(xx.mus.go,"pc20.k50.res02.cluster_fc025.pct025.pos.markergene_pqv001.go_220622.rds")
df<-xx.mus.go@compareClusterResult
df1 = df %>% group_by(.,Cluster) %>% filter(.,Count>=5) %>% top_n(.,-5,p.adjust)  
df1$Description = factor(df1$Description,levels = rev(unique(df1$Description)))
options(repr.plot.width=7.5, repr.plot.height=7)
p=ggplot(df1,aes(Cluster,Description))+geom_point(aes(size=Count,color=p.adjust))+        
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+
    theme(axis.title = element_blank(),axis.text = element_text(size = 13))
p
ggsave("WT.SCT.reg.pc20.k50.res02.cluster_fc025.pct025.pos.markergene_pqv001.count5.pval.top5.dotplot_220622.png",width = 7.5,height = 7,dpi = 400)

###time DEG according to domain
Idents(da)<-"time"
wmda<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="WM",])]
mgda<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="MG",])]
vhda<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="VH",])]
dhda<-da[rownames(da),rownames(da@meta.data[da@meta.data$domain_res02=="DH",])]
HT<-c(wmda,mgda,dhda,vhda)
name<-c("WM","MG","DH","VH")
sd<-"time_domain.marker/"
if(!dir.exists(sd))
    dir.create(sd)
for(i in c(1:4)){
    temp<-HT[[i]]
    mark<-FindAllMarkers(temp,logfc.threshold = 0.5,min.pct = 0.25,only.pos = F)
    mark$up_down<-NA
    mark[mark$avg_log2FC<0,"up_down"]<-"down"
    mark[mark$avg_log2FC>0,"up_down"]<-"up"
    mark$group<-paste(mark$cluster,name[i],mark$up_down,sep = "-")
    mark$group<-factor(mark$group,
                        levels = c(paste0('WT_sham-',name[i],'-up'),paste0('WT_3h-',name[i],'-up'),paste0('WT_24h-',name[i],'-up'),paste0('WT_72h-',name[i],'-up'),
                                   paste0('WT_sham-',name[i],'-down'),paste0('WT_3h-',name[i],'-down'),paste0('WT_24h-',name[i],'-down'),paste0('WT_72h-',name[i],'-down')    
                                                  ))
    mark <- filter(mark,p_val_adj<0.01)
    write.csv(mark,paste0(sd,"WT.SCT.",name[i],".time_Marker.p001.csv"))
    mark_up<-mark[mark$up_down=="up",]
    mark_up$group<-factor(mark_up$group,levels = levels(mark$group)[1:4])
    write.csv(mark_up,paste0(sd,"WT.SCT.",name[i],".time_upMarker.p001.csv"))
    mark_up %>%
        group_by(group) %>% #filter(pct.1>0.5) %>%
        top_n(n = 10, wt = avg_log2FC) -> top10
    options(repr.plot.width=12,repr.plot.height=12)
    DoHeatmap(temp, features = top10$gene) + NoLegend()
    ggsave(paste0(sd,name[i],".time_Marker.pos.pct025.top10_fc.heatmap.png"),width = 12,height = 12,dpi = 400)
    mark_up %>%
        group_by(group) %>% #filter(pct.1>0.5) %>%
        top_n(n = 30, wt = avg_log2FC) -> top30
    options(repr.plot.width=10,repr.plot.height=15)
    DoHeatmap(temp, features = top30$gene) + NoLegend()
    ggsave(paste0(sd,name[i],".time_Marker.pos.pct025.top30_fc.heatmap.png"),width = 10,height = 15,dpi = 400)           
}

for(j in c("WM","MG","VH","DH")){
    temp2<-HT[[j]]
    #name<-names(HT[j])
    for(a in c(levels(Idents(temp2))[2:4])){
        mark<-FindMarkers(temp2,logfc.threshold = 0.5,min.pct = 0.25,only.pos = F,ident.1=a,ident.2="WT_sham")
            
        mark <- filter(mark,p_val_adj<0.01)
        dim(mark)
        mark$group<-paste0(a,"_",j)
        table(mark$group)
        write.csv(mark,paste0(sd,"WT.SCT.",j,"_",a,".timevs.sham_Marker.np.p001.csv"))
    }
}

##merge domain deg and number barplot
sf<-"time_domain.marker/"
f<-dir(sf,pattern = "WT.SCT.*.np.*csv")
ht<-read.csv(paste0(sf,f[1]))
for(i in 2:length(f)){
    temp<-read.csv(paste0(sf,f[i]))
    ht<-rbind(ht,temp)
}
write.csv(ht,"time_domain.marker/WT.domain.timevs.sham.fc05.pct025.p001.np.marker.csv")
#domain bulk
ht$change<-"specific"
deg_3h_com<-Reduce(intersect,list(ht[ht$group=='WT_3h_WM',"X"],ht[ht$group=='WT_3h_MG',"X"],
                      ht[ht$group=='WT_3h_VH',"X"],ht[ht$group=='WT_3h_DH',"X"]))
ht[ht$group %in% c('WT_3h_WM','WT_3h_MG','WT_3h_VH','WT_3h_DH')& ht$X %in% deg_3h_com,"change"]<-"common"
deg_24h_com<-Reduce(intersect,list(ht[ht$group=='WT_24h_WM',"X"],ht[ht$group=='WT_24h_MG',"X"],
                      ht[ht$group=='WT_24h_VH',"X"],ht[ht$group=='WT_24h_DH',"X"]))
ht[ht$group %in% c('WT_24h_WM','WT_24h_MG','WT_24h_VH','WT_24h_DH')& ht$X %in% deg_24h_com,"change"]<-"common"
deg_72h_com<-Reduce(intersect,list(ht[ht$group=='WT_72h_WM',"X"],ht[ht$group=='WT_72h_MG',"X"],
                      ht[ht$group=='WT_72h_VH',"X"],ht[ht$group=='WT_72h_DH',"X"]))
ht[ht$group %in% c('WT_72h_WM','WT_72h_MG','WT_72h_VH','WT_72h_DH')& ht$X %in% deg_72h_com,"change"]<-"common"
ht$up_down<-NA
ht$up_down<-ifelse(ht$avg_log2FC>0,"up","down")
ht$time<-"sham"
ht[grepl("3h",ht$group,fixed = T) ,"time"]<-"3h"
ht[grepl("24h",ht$group,fixed = T) ,"time"]<-"24h"
ht[grepl("72h",ht$group,fixed = T) ,"time"]<-"72h"
ht$time<-factor(ht$time,levels = c("3h","24h","72h"))
write.csv(ht,"time_domain.marker/WT.domain.timevs.sham.fc05.pct025.p001.np.marker_ad.csv")

#barplot to show the specific and common number of DEGs
temp<-ht
temp$value<-ifelse(temp$up_down=="up",1,-1)
temp2<-temp %>% group_by(group,time,up_down,change) %>% summarise(total_count=sum(value))
temp2$HT<-"WM"
temp2[grepl("VH",temp2$group),"HT"]<-"VH"
temp2[grepl("DH",temp2$group),"HT"]<-"DH"
temp2[grepl("MG",temp2$group),"HT"]<-"MG"
temp2$change<-factor(temp2$change,levels = c("specific","common"))
temp2$up_down<-factor(temp2$up_down,levels = c("up","down"))
temp2$group<-factor(temp2$group,levels = c('WT_3h_DH','WT_3h_MG','WT_3h_VH','WT_3h_WM',
                               'WT_24h_DH','WT_24h_MG','WT_24h_VH','WT_24h_WM',
                               'WT_72h_DH','WT_72h_MG','WT_72h_VH','WT_72h_WM'))
temp2$fill<-paste0(temp2$up_down,"_",temp2$change)
temp2$fill<-factor(temp2$fill,levels = c('up_specific','up_common','down_specific','down_common'))
write.csv(temp2,"time_domain.marker/domain.time.vs.sham.barplot.number.csv")
temp2$group<-factor(temp2$group,levels = c('WT_3h_WM','WT_3h_MG','WT_3h_DH','WT_3h_VH',
                               'WT_24h_WM','WT_24h_MG','WT_24h_DH','WT_24h_VH',
                                           'WT_72h_WM','WT_72h_MG',
                               'WT_72h_DH','WT_72h_VH'))
temp2$fill<-factor(temp2$fill,levels = c('down_specific','down_common','up_specific','up_common'))                               
pdf("time_domain.marker/domain_time.vs.sham.np.deg.barplot_221029.pdf",width = 7,height = 6)
options(repr.plot.width=7,repr.plot.height=6)
ggplot(temp2,aes(x=group,y=total_count,fill=fill))+
    geom_bar(position = "stack",stat="identity",alpha=0.8)+
    scale_fill_manual(values = c('#2B83BA','#ABDDA4','#D7191C','#FDAE61'))+
    theme_bw()+
    theme(panel.grid  =element_blank(),axis.text.x.bottom = element_text(angle = 90),axis.line.x = element_line(colour = "black"))+
    coord_cartesian(ylim = c(-600,600))
dev.off()

#GO term
temp<-read.csv("time_domain.marker/WT.domain.timevs.sham.fc05.pct025.p001.np.marker.csv")
temp<-temp[,-1]
temp$time<-paste0(str_split(temp$group,"_",simplify = T)[,1],"_",str_split(temp$group,"_",simplify = T)[,2])
temp$up_down<-""
temp$up_down[temp$avg_log2FC>0]<-"up"
temp$up_down[temp$avg_log2FC<0]<-"down"
temp$group<-paste0(temp$group,"-",temp$up_down)
table(temp$group)
write.csv(temp,paste0(sd,"domain.time.bulk.vs.sham.fc05.pct025.pad001.rbind.csv"))
grlabs<-split(temp$X,temp$group)
gcSample = lapply(grlabs, function(gr) as.numeric(bitr(gr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID))
pvalueCutoff = 0.01
qvalueCutoff = 0.01
xx.mus.go <- compareCluster(gcSample, OrgDb='org.Mm.eg.db', fun='enrichGO',#maxGSSize=150,
                            pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                            ont = "BP", readable=T,universe=backgenes
                           )
EnrichResGO_mm = xx.mus.go@compareClusterResult
saveRDS(xx.mus.go,paste0(sd,"WT.domain.time.vs.sham.fc05.pct025.p001.updown_pqval001.background.go.rds"))
df<-xx.mus.go@compareClusterResult
df$time<-paste0(str_split(df$Cluster,"_",simplify = T)[,1],"_",str_split(df$Cluster,"_",simplify = T)[,2])
df$time<-factor(df$time,levels = c("WT_3h","WT_24h","WT_72h"))
df<-df[order(df$time),]
df$up_down<-str_split(df$Cluster,"-",simplify = T)[,2]
df$up_down<-factor(df$up_down,levels = c("up","down"))
df<-df[order(df$up_down),]
write.csv(df,"WT.domain.time.vs.sham.fc05.pct025.p001.updown_pqval001.background.go.adjust.dataframe.csv")
df2<-df[grep("DH",df$Cluster),]
df1 = df1 %>% group_by(.,Cluster) %>% filter(Count>=5) %>% top_n(.,-5,p.adjust)  
df1<-df2[df2$Description %in% unique(df1$Description),]
df1$Cluster<-factor(df1$Cluster,levels = c('WT_3h_DH-up','WT_3h_MG-up','WT_3h_VH-up','WT_3h_WM-up',
                                           'WT_24h_DH-up','WT_24h_MG-up','WT_24h_VH-up','WT_24h_WM-up',
                                           'WT_72h_DH-up','WT_72h_MG-up','WT_72h_VH-up','WT_72h_WM-up',
                                           'WT_3h_DH-down','WT_3h_MG-down','WT_3h_VH-down','WT_3h_WM-down',
                                           'WT_24h_DH-down','WT_24h_MG-down','WT_24h_VH-down','WT_24h_WM-down',
                                           'WT_72h_DH-down','WT_72h_MG-down','WT_72h_VH-down','WT_72h_WM-down'
                                          ))
df1<-df1[order(df1$time),]
df1<-df1[order(df1$up_down),]
df1$Description = factor(df1$Description,
                         levels = unique(c(unique(df1[df1$up_down=="up",]$Description),rev(unique(df1[df1$up_down=="down",]$Description)))))  
options(repr.plot.width=6, repr.plot.height=4,font.size=15)
library(ggplot2)
p=ggplot(df1,aes(time,Description))+geom_point(aes(size=Count,color=p.adjust))+        
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+facet_wrap(~up_down,nrow=2,scales = "free")+
    theme(axis.title = element_blank(),axis.text = element_text(size = 13),axis.text.x.bottom = element_text(angle = 90))
p
ggsave(paste0(sd,"WT.DH.time.Marker.fc05.pct025.p001_GO.pqval001.C5.pval.top5.go_updown.png"),width = 6,height = 8.5,dpi = 400)

df2<-df[grep("MG",df$Cluster),]
df1 = df1 %>% group_by(.,Cluster) %>% filter(Count>=5) %>% top_n(.,-5,p.adjust) 
df1<-df2[df2$Description %in% unique(df1$Description),]
df1<-df1[order(df1$time),]
df1<-df1[order(df1$up_down),]
df1$Description = factor(df1$Description,
                         levels = unique(c(unique(df1[df1$up_down=="up",]$Description),rev(unique(df1[df1$up_down=="down",]$Description))))) 
options(repr.plot.width=6.2, repr.plot.height=8.5,font.size=15)
library(ggplot2)
p=ggplot(df1,aes(time,Description))+geom_point(aes(size=Count,color=p.adjust))+
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+facet_wrap(~up_down,nrow=2,scales = "free")+
    theme(axis.title = element_blank(),axis.text = element_text(size = 13),axis.text.x.bottom = element_text(angle = 90))
p
ggsave(paste0(sd,"WT.MG.time.Marker.fc05.pct025.p001_GO.pqval001.C5.pval.top5.go_updown.png"),width = 6.2,height = 8.5,dpi = 400)                         
df2<-df[grep("VH",df$Cluster),]
df1 = df1 %>% group_by(.,Cluster) %>% filter(Count>=5) %>% top_n(.,-5,p.adjust) 
df1<-df2[df2$Description %in% unique(df1$Description),]
df1<-df1[order(df1$time),]
df1<-df1[order(df1$up_down),]
df1$Description = factor(df1$Description,
                         levels = unique(c(unique(df1[df1$up_down=="up",]$Description),rev(unique(df1[df1$up_down=="down",]$Description))))) 
options(repr.plot.width=6, repr.plot.height=8.5,font.size=15)
library(ggplot2)
p=ggplot(df1,aes(time,Description))+geom_point(aes(size=Count,color=p.adjust))+ 
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+facet_wrap(~up_down,nrow=2,scales = "free")+
    theme(axis.title = element_blank(),axis.text = element_text(size = 13),axis.text.x.bottom = element_text(angle = 90))
p
ggsave(paste0(sd,"WT.VH.time.Marker.fc05.pct025.p001_GO.pqval001.C5.pval.top5.go_updown.png"),width = 6,height = 8.5,dpi = 400)

df2<-df[grep("WM",df$Cluster),]
df1 = df1 %>% group_by(.,Cluster) %>% filter(Count>=5) %>% top_n(.,-5,p.adjust) 
df1<-df2[df2$Description %in% unique(df1$Description),]
df1<-df1[order(df1$time),]
df1<-df1[order(df1$up_down),]
df1$Description = factor(df1$Description,
                         levels = unique(c(unique(df1[df1$up_down=="up",]$Description),rev(unique(df1[df1$up_down=="down",]$Description))))) 
options(repr.plot.width=6, repr.plot.height=8.5,font.size=15)
library(ggplot2)
p=ggplot(df1,aes(time,Description))+geom_point(aes(size=Count,color=p.adjust))+ 
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+facet_wrap(~up_down,nrow=2,scales = "free")+
    theme(axis.title = element_blank(),axis.text = element_text(size = 13),axis.text.x.bottom = element_text(angle = 90))
p
ggsave(paste0(sd,"WT.WM.time.Marker.fc05.pct025.p001_GO.pqval001.C5.pval.top5.go_updown.png"),width = 6,height = 8.5,dpi = 400)



###HT DEGs
Hda<-da[rownames(da),rownames(da@meta.data[da@meta.data$RC=="WT_H",])]
Tda<-da[rownames(da),rownames(da@meta.data[da@meta.data$RC=="WT_T",])]
HT<-c(Hda,Tda)
names(HT)<-c("H","T")
for(j in c("H","T")){
    temp2<-HT[[j]]
    #name<-names(HT[j])
    for(a in c(levels(Idents(temp2))[2:4])){
        mark<-FindMarkers(temp2,logfc.threshold = 0.5,min.pct = 0.25,only.pos = F,ident.1=a,ident.2="WT_sham")
            
        mark <- filter(mark,p_val_adj<0.01)
        dim(mark)
        mark$group<-paste0(a,"_",j)
        table(mark$group)
        write.csv(mark,paste0(sd,j,"_",a,".timevs.sham_Marker.fc05.pct025.p001.csv"))
    }
}
f<-dir(sd,pattern = "^[H-T]")
temp<-read.csv(paste0(sd,f[1]))
for(i in 2:length(f)){
    temp2<-read.csv(paste0(sd,f[i]))
    temp<-rbind(temp,temp2)
}
temp$up_down<-""
temp$up_down[temp$avg_log2FC>0]<-"up"
temp$up_down[temp$avg_log2FC<0]<-"down"
temp$time<-temp$group
temp$group<-paste0(temp$time,"-",temp$up_down)
write.csv(temp,paste0(sd,"HT.time.bulk.vs.sham.fc05.pct025.pad001.rbind.csv"))

###Distance
Hda<-da[rownames(da),rownames(da@meta.data[da@meta.data$distance=="1mm",])]
Tda<-da[rownames(da),rownames(da@meta.data[da@meta.data$distance=="0.5mm",])]
HT<-c(Hda,Tda)
names(HT)<-c("1mm","0.5mm")
for(j in c("1mm","0.5mm")){
    temp2<-HT[[j]]
    for(a in c(levels(Idents(temp2))[2:4])){
        mark<-FindMarkers(temp2,logfc.threshold = 0.5,min.pct = 0.25,only.pos = F,ident.1=a,ident.2="WT_sham")
            
        mark <- filter(mark,p_val_adj<0.01)
        dim(mark)
        mark$group<-paste0(a,"_",j)
        table(mark$group)
        write.csv(mark,paste0(sd,j,"_",a,".timevs.sham_Marker.fc05.pct025.p001.csv"))
    }
}
f<-dir(sd,pattern = "^[0-1]")
temp<-read.csv(paste0(sd,f[1]))
for(i in 2:length(f)){
    temp2<-read.csv(paste0(sd,f[i]))
    temp<-rbind(temp,temp2)
}
temp$up_down<-""
temp$up_down[temp$avg_log2FC>0]<-"up"
temp$up_down[temp$avg_log2FC<0]<-"down"
temp$time<-temp$group
temp$group<-paste0(temp$time,"-",temp$up_down)
write.csv(temp,paste0(sd,"distance.time.bulk.vs.sham.fc05.pct025.pad001.rbind.csv"))
