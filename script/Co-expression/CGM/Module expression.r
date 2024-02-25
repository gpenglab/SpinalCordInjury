library(Seurat)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
options(connectionObserver=NULL)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(ggridges)
library(ggsci)
library(scales)
source("script/self_function/save_pheatmap_pdf.R")

### Module seurat score correlation with trait
da<-readRDS("WT_replace_v2/res02_220310/WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
modgene<-read.csv("thr100.sp18.hclust_average.min30.deep2.h015.34modules.full.gene.list.csv")
rownames(modgene)<-modgene$gene
dim(modgene)
modgene_filt<-modgene[!modgene$module=="grey",]
##calculate module gene score
for(i in unique(modgene_filt$module)){
    genes<-modgene_filt[modgene_filt$module==i,"gene"]
    da<-AddModuleScore(da,features = list(genes),name = i)
}
module_score<-da@meta.data[,c(which(colnames(da@meta.data)=="coral1"):ncol(da@meta.data))]
colnames(module_score)<-gsub("1$","",colnames(module_score))
write.csv(module_score,"WT.domain.thr100.sp18.hclust_average.min30.deep2.h015.33modules.SCT.score.ex_grey.csv")
name<-colnames(da@meta.data)[which(colnames(da@meta.data)=="coral1"):ncol(da@meta.data)]
meta<-da@meta.data
sf<-"modulescore_spatial/"
if(!dir.exists(sf))
    dir.create(sf)
options(repr.plot.width=20,repr.plot.height=16)
for(i in name){
    ggplot(meta,aes_string("x.ad","y.ad",color=i))+geom_point(size=1)+
        #scale_colour_gradientn(colours=colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
        ggplot2::scale_colour_gradientn(colours = c('#0000ff','#3232ff','#6666ff','#9898ff','#ccccff',
                                              '#fffefe','#ffcccc','#ff9898','#ff6666','#ff3232','#ff0000'))+
        xlab(paste0("")) +
        ylab(paste0("")) + 
        theme(panel.background = element_blank(),
              panel.grid.major =element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank()
          # axis.ticks.x = element_blank(),
    )#+
    #plot(p)
    #p<-SpatialFeaturePlot(mergda,features = c(i),min.cutoff=mean(mergda@meta.data[,i]),combine = F)
    #do.call(ggpubr::ggarrange,c(p,list(ncol=4,nrow=8,common.legend=T,legend="right")))
    ggsave(paste0(sf,"WT_",i,".module.score.spatial.png"),width=20,height=16,dpi=300)
    }

##ridgeplot of each module at time_domain level
df<-da@meta.data
colnames(df)<-gsub("1$","",colnames(df))
df$time<-factor(df$time,levels = rev(c("WT_sham","WT_3h","WT_24h","WT_72h")))
col<-c('#BC3C29A8','#0072B5A8','#E18727A8','#20854EA8')
show_col(col)
names(col)<-c('DH','MG','VH','WM')
name<-gsub("1$","",name)
name
options(repr.plot.width=6,repr.plot.height=4)
fig1<-"moduleScore_ridgePlot"
if(!dir.exists(fig1))
    dir.create(fig1)
for(i in name){
    ggplot(df,aes_string(x=i,y="time",color="domain",fill="domain"))+
    geom_density_ridges(color="white")+
    scale_y_discrete(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0),name = paste0(i," module score")#,limits=c(min(df[,"lightcyan1"]),max(df[,"lightcyan1"]))
                      )+
    scale_fill_manual(values = c(alpha(col,0.66)),labels=c("DH","MG","VH","WM"))+
    coord_cartesian(clip = "off")+
    guides(fill=guide_legend(override.aes=list(fill=col,
                                               color=NA,points_color=NA)))+
    #ggtitle(paste0(i," module score at domain after SCI"))+
    theme_ridges(center=TRUE)+
    theme_ridges(center=TRUE)+
    theme(panel.grid.major = element_blank(),axis.title.y.left = element_blank())
    ggsave(paste0(fig1,"/",i,"_moduleScore.ridgePlot.png"),width = 6,height = 4,dpi = 300)
}

#1.time_domain correlation and pval
mod_dom_mean<-module_score
mod_dom_mean$time_domain<-paste0(meta$time,"_",meta$domain_res02)
mod_dom_mean<-aggregate(mod_dom_mean[,-which(colnames(mod_dom_mean)=="time_domain")],list(mod_dom_mean$time_domain),mean)
rownames(mod_dom_mean)<-mod_dom_mean[,1]
mod_dom_mean<-mod_dom_mean[,-1]

trait<-data.frame(row.names = rownames(mod_dom_mean))
for(i in rownames(trait)){
    #trait[,i]<-0
    for(j in 1:nrow(trait)){
        trait[j,i]<-ifelse(rownames(trait)[j]==i,1,0)
    }
}
trait
c<-cor(mod_dom_mean,trait)
pval<-corPvalueStudent(c,16)
pval
c_long<-melt(c)
c_long$Var2<-factor(c_long$Var2,
                    levels = c('WT_sham_WM','WT_3h_WM','WT_24h_WM','WT_72h_WM',
        'WT_sham_MG','WT_3h_MG','WT_24h_MG', 'WT_72h_MG',
        'WT_sham_DH','WT_3h_DH','WT_24h_DH', 'WT_72h_DH',
        'WT_sham_VH','WT_3h_VH','WT_24h_VH', 'WT_72h_VH')
                    #c('WT_sham_DH','WT_sham_MG','WT_sham_VH','WT_sham_WM',
                             #  'WT_3h_DH','WT_3h_MG','WT_3h_VH','WT_3h_WM',
                             #  'WT_24h_DH','WT_24h_MG','WT_24h_VH','WT_24h_WM',
                             #  'WT_72h_DH','WT_72h_MG','WT_72h_VH','WT_72h_WM'
                             # )
                   )
pval_long<-melt(pval)
pval_long$c<-c_long$value
pval_long$l<-paste0("(",as.character(round(pval_long$c,2)),")","\n",as.character(round(pval_long$value,2)))
c<-c[,c('WT_sham_WM','WT_3h_WM','WT_24h_WM','WT_72h_WM',
        'WT_sham_MG','WT_3h_MG','WT_24h_MG', 'WT_72h_MG',
        'WT_sham_DH','WT_3h_DH','WT_24h_DH', 'WT_72h_DH',
        'WT_sham_VH','WT_3h_VH','WT_24h_VH', 'WT_72h_VH'
        )]
c_df<-as.data.frame(c)
c_df$module<-rownames(c_df)
c_df<-melt(c_df)
df1 = c_df %>% group_by(.,variable) %>% top_n(.,5,value)
c_ad<-c[rev(unique(df1$module)),]
row_anno<-data.frame(row.names = rownames(c_ad),"module"=rownames(c_ad))
row_anno
row_color<-row_anno$module
names(row_color)<-row_anno$module
col_anno<-data.frame(row.names = colnames(c_ad),"domain"=c(rep("WM",4),rep("MG",4),rep("DH",4),rep("VH",4)),
                    "time"=c(rep(c("sham","3h","24h","72h"),4)))
col_anno
domain_col<-c('#20854EA8','#0072B5A8','#BC3C29A8','#E18727A8')
names(domain_col)<-c("WM", "MG", "DH", "VH")
time_col<-c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF')
names(time_col)<-c("sham","3h","24h","72h")
col_list<-list(module=row_color,domain=domain_col,time=time_col)
gap_col<-c(4,8,12)
#pdf("module-time_domain relationships.pheatmap_220530.pdf")
options(repr.plot.width=8,repr.plot.height=12)
p<-pheatmap(c_ad,
            col=colorRampPalette(colors = c('#0000ff',
                                            #'#440154FF',
                                            '#fffefe',
                                            '#ff0000'
                                            #'#FDE725FF'
            ))(50),
         scale="none",
         cluster_cols=F,cluster_rows=F,#treeheight_row=10,cutree_rows=7,
         annotation_row=row_anno,
         annotation_col=col_anno,
         border_color="NA",
         fontsize=6,
         annotation_colors=col_list,
         legend=T,gaps_col=gap_col,
         #cellwidth=25,cellheight=25,
         show_rownames=T)
save_pheatmap_pdf(p,filename = "module-time_domain relationships.pheatmap_220721.pdf",width = 8,height = 12)
#dev.off()

###module go term
WGCNA.modules.full.gene.list<-read.csv("thr100.sp18.hclust_average.min30.deep2.h015.34modules.full.gene.list.csv")
modgene<-WGCNA.modules.full.gene.list[WGCNA.modules.full.gene.list$module!="grey",]
length(unique(modgene$module))
backid<-read.table("WT_replace_v2/WT.GOterm.SCT.backgeneids.txt")
backid<-as.character(backid$x)
grlabs<-split(modgene$gene,modgene$module)
gcSample<-lapply(grlabs,function(gr) as.numeric(bitr(gr,fromType="SYMBOL",toType = "ENTREZID",OrgDb="org.Mm.eg.db")$ENTREZID))
pvalueCutoff=0.1
qvalueCutoff=0.1
xx.mus.go<-compareCluster(gcSample,OrgDb="org.Mm.eg.db",fun="enrichGO",universe=backid,
                          pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff,
                         ont="MF",readable=T)
saveRDS(xx.mus.go,"thr100.sp18.min30.deep2.h015.33modules.module.gene.go.pqval01_MF.220805.rds")
df<-xx.mus.go@compareClusterResult
df<-df %>% filter(p.adjust<0.05 & qvalue<0.05)
write.csv(df,"25module.pqval005.ALL.Goterm.dataframe_220805.csv")
df1<-df %>% filter(p.adjust<0.05 & qvalue<0.05) %>% group_by(.,Cluster) %>% top_n(.,-3,p.adjust)%>% top_n(.,3,Count) # order by cluster and select each top10 terms
df1$Description<-factor(df1$Description,levels=rev(unique(df1$Description))) #order Description by cluster 
col1<-unique(as.character(df1$Cluster))
col_darken<-darken(col1,0.15,space = "HLS")
### barplot
options(repr.plot.width=12,repr.plot.height=15)
ggplot(df1,aes(x=Description,y=-log10(p.adjust),
              fill=Cluster))+
    geom_bar(stat="identity",position = "dodge")+
    coord_flip()+
    scale_fill_manual(values = alpha(col_darken,0.88))+
    theme_bw()+
    #geom_text(aes(label=Count),hjust=1,size=3.5,color="white")+
    theme(axis.text.y.left = element_text(size = 10),panel.background = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor=element_blank())
ggsave("25module.pqval005.top3.Goterm.barplot_220810.png",width = 12,height = 15,dpi = 300)

##representative go terms
temp<-df[df$Cluster=="cyan" & df$Description %in% c("positive regulation of angiogenesis","response to wounding",
                                                    #"myeloid leukocyte migration",
                                                    "leukocyte chemotaxis","response to interferon-gamma"),]
temp<-temp[order(temp$p.adjust),]
temp$Description<-factor(temp$Description,levels = rev(unique(temp$Description)))
options(repr.plot.width=6,repr.plot.height=3)
ggplot(temp,aes(x=Description,y=-log10(p.adjust),
              fill=Cluster))+
    geom_bar(stat="identity",position = "dodge")+
    coord_flip()+
    scale_fill_manual(values = alpha(darken("cyan",0.15,space = "HLS"),0.88))+
    theme_bw()+
    #geom_text(aes(label=Count),hjust=1,size=3.5,color="white")+
    theme(axis.text.y.left = element_text(size = 10),panel.background = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor=element_blank())
ggsave("cyan.selected.Goterm.barplot_220822.png",width = 6,height = 3,dpi = 300)

# darkseagreen4(13)
temp<-df[df$Cluster=="darkseagreen4" & df$Description %in% c("ERK1 and ERK2 cascade",#"regulation of tumor necrosis factor production",
                                                    "response to lipopolysaccharide",
                                                    "cellular response to interleukin-1","tumor necrosis factor production"),]
temp<-temp[order(temp$p.adjust),]
temp$Description<-factor(temp$Description,levels = rev(unique(temp$Description)))
options(repr.plot.width=6,repr.plot.height=3)
ggplot(temp,aes(x=Description,y=-log10(p.adjust),
              fill=Cluster))+
    geom_bar(stat="identity",position = "dodge")+
    coord_flip()+
    scale_fill_manual(values = alpha(darken("darkseagreen4",0.15,space = "HLS"),0.88))+
    theme_bw()+
    #geom_text(aes(label=Count),hjust=1,size=3.5,color="white")+
    theme(axis.text.y.left = element_text(size = 10),panel.background = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor=element_blank())
ggsave("darkseagreen4.selected.Goterm.barplot_220822.png",width = 6,height = 3,dpi = 300)

# darkseagreen4(13)
temp<-df[df$Cluster=="darkorange" & df$Description %in% c("synapse organization",#"regulation of tumor necrosis factor production",
                                                    "oxidative phosphorylation",
                                                    "synaptic vesicle cycle","neurotransmitter secretion"),]
temp<-temp[order(temp$p.adjust),]
temp$Description<-factor(temp$Description,levels = rev(unique(temp$Description)))
options(repr.plot.width=6,repr.plot.height=3)
ggplot(temp,aes(x=Description,y=-log10(p.adjust),
              fill=Cluster))+
    geom_bar(stat="identity",position = "dodge")+
    coord_flip()+
    scale_fill_manual(values = alpha(darken("darkorange",0.15,space = "HLS"),0.88))+
    theme_bw()+
    #geom_text(aes(label=Count),hjust=1,size=3.5,color="white")+
    theme(axis.text.y.left = element_text(size = 10),panel.background = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor=element_blank())
ggsave("darkorange.selected.Goterm.barplot_220822.png",width = 6,height = 3,dpi = 300)






