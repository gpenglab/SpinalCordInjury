library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(pheatmap)
library(venn)
library(VennDiagram)
library(clusterProfiler)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
# For some of the plots:
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
source("script/self_function/pheatmap_add_flag.R")
source("script/self_function/save_pheatmap_pdf.R")
source("CSI.R")
source("script/self_function/aux_rss.R")

###hclust 317 astrocyte deconvolution marker in spot level
#Astrocyte-Gfap normalized spot expression data
da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
meta<-read.csv("fcreg1.32cp.conf_nor.astr_state.meta.csv")
rownames(meta)<-meta[,1]
meta<-meta[,-1]
## add specific state label to each spot including WM_Gfap, GM_Gfap, GM_Slc7a10 (exclude each other and Svep1 and AstroEpendymal)
meta$Astr_state2<-""
meta[meta$Astro.Svep1==0&meta$Astrocyte.Gfap>0&meta$Astrocyte.Slc7a10==0&meta$Astroependymal==0,"Astr_state2"]<-"Astr_Gfap"
meta[meta$Astro.Svep1==0&meta$Astrocyte.Gfap==0&meta$Astrocyte.Slc7a10>0&meta$Astroependymal==0,"Astr_state2"]<-"Astr_Slc"
cells<-rownames(meta)[meta$Astr_state2 %in% c("Astr_Gfap","Astr_Slc")]
ma<-CreateSeuratObject(counts = da@assays$SCT@scale.data[astr_mark,cells],meta.data = meta[cells,])
saveRDS(ma,"astr.prop0.317markers.SCTscale.data.meta.rds")
# subset sample and Astrocyte_gfap
meta<-ma@meta.data
meta_sub<-subset(meta,sample %in% c('WT_sham_H_R2_1mm_2','WT_sham_H_R2_4','WT_sham_T_210323_4','WT_sham_T_R2_1mm_4',
                                                  'WT_3h_H_R2_1mm_2','WT_3h_H_R2_2','WT_3h_T_R2_1','WT_3h_T_210330_1mm_3',
                                                  'WT_24h_H_R1_1mm_3','WT_24h_H_R1_4','WT_24h_T_201231_3','WT_24h_T_R1_1mm_3',
                                                  'WT_72h_H_R1_1mm_2','WT_72h_H_210323_1','WT_72h_T_R2_4','WT_72h_T_R1_1mm_3')& Astr_state2=="Astr_Gfap")                                                  
cells<-rownames(meta_sub)
exp<-as.matrix(ma@assays$RNA@counts[,cells])
astr_marker<-rownames(exp)
exp_bi<-ifelse(exp>0,1,-1)
meta_sub$domain2<-ifelse(meta_sub$domain_res02=="WM","WM","GM")
### time_distance_domain bulk
sp_df<-data.frame()
for(i in unique(meta_sub$orig.ident)){
    if(grepl("sham",i) | grepl("3h",i)){
        temp_me<-meta_sub[meta_sub$Astr_state%in%c("WM_Astr_Gfap"#,"MG_Astr_Slc"
                                                  )&
                            meta_sub$orig.ident==i,]
        #table(temp_me$Astr_state)
    
        temp_ma<-as.data.frame(t(exp_bi[,rownames(temp_me)]))
        temp_ma$group<-temp_me$domain_res02
        temp_mean<-temp_ma %>% group_by(group) %>% summarise_all("mean")
        temp_mean$section<-i
        temp_mean<-as.data.frame(temp_mean)
    }
    if(grepl("24",i) | grepl("72h",i)){
        temp_me<-meta_sub[meta_sub$Astr_state%in%c("WM_Astr_Gfap",#"MG_Astr_Slc",
                                           "MG_Astr_Gfap")
                           & meta_sub$orig.ident==i,]
        table(temp_me$Astr_state)
    
        temp_ma<-as.data.frame(t(exp_bi[,rownames(temp_me)]))
        temp_ma$group<-temp_me$domain_res02
        temp_mean<-temp_ma %>% group_by(group) %>% summarise_all("mean")
        temp_mean$section<-i
        temp_mean<-as.data.frame(temp_mean)
    }
    sp_df<-rbind(sp_df,temp_mean)
    
}

dim(sp_df)
sp_df$key<-paste0(sp_df$section,"_",sp_df$group)
unique(sp_df$key)
sp_df$section<-factor(sp_df$section,levels = c('WT_sham_H_R2_1mm','WT_sham_H_R2','WT_sham_T_210323','WT_sham_T_R2_1mm',
                                              'WT_3h_H_R2_1mm','WT_3h_H_R2','WT_3h_T_R2','WT_3h_T_210330_1mm',
                                              'WT_24h_H_R1_1mm','WT_24h_H_R1','WT_24h_T_201231','WT_24h_T_R1_1mm',
                                              'WT_72h_H_R1_1mm','WT_72h_H_210323','WT_72h_T_R2','WT_72h_T_R1_1mm'))

sp_df$time<-NA
sp_df$time[grepl("sham",sp_df$section)]<-"sham"
sp_df$time[grepl("3h",sp_df$section)]<-"3h"
sp_df$time[grepl("24h",sp_df$section)]<-"24h"
sp_df$time[grepl("72h",sp_df$section)]<-"72h"
table(sp_df$time)
sp_df$RC<-NA
sp_df$RC[grepl("_H_",sp_df$section)]<-"rostral"
sp_df$RC[grepl("_T_",sp_df$section)]<-"caudal"
table(sp_df$RC)
sp_df$distance<-NA
sp_df$distance<-ifelse(grepl("_1mm",sp_df$section),"1mm","0.5mm")
table(sp_df$distance)
sp_df$domain<-sp_df$group
sp_df$time<-factor(sp_df$time,levels = c("sham","3h","24h","72h"))
sp_df$RC<-factor(sp_df$RC,levels = c("rostral","caudal"))
sp_df$distance<-factor(sp_df$distance,levels = c("0.5mm","1mm"))
sp_df$domain<-factor(sp_df$domain,levels = c("WM","MG")#c("WM","MG","DH","VH")
                    )

write.csv(sp_df,"astr_gfap.317marker.scale.express.binary.mean.in.time_Ct2domain.meta.csv")

ha1<-columnAnnotation(time=sp_df$time,
                       RC=sp_df$RC,
                       distance=sp_df$distance,
                       domain=sp_df$domain,
                     col=list(time=structure(names=c("sham","3h","24h","72h"),c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF")),
                             RC=structure(names=c("rostral","caudal"),c('#F4A582','#0571B0')),
                             distance=structure(names=c("1mm","0.5mm"),c('#63A79C','#FB9A99')),
                             domain=structure(names=c("WM","MG"),#c("WM","MG","DH","VH"),
                                              c("#20854EA8","#0072B5A8")#c("#20854EA8","#0072B5A8","#BC3C29A8","#E18727A8")
                                             )))

temp<-t(sp_df[,astr_marker])
options(repr.plot.width=10,repr.plot.height=13)
set.seed(220609)
pdf("gfap.317.markers.scale.binary.mean.time_ct2domain.clustered.k5.heatmap.pdf",width = 8,height = 10)
Heatmap(temp,cluster_columns = T,show_row_names = F,km = 5,#row_names_max_width = unit(2,"cm"),
        cluster_rows = T,
        col = c("#05445E","#92a8d1","#F6E6E8","#E56997","#B91646"),
        #clustering_distance_rows = "pearson",
        show_column_names = F,
        #left_annotation = ha,
        top_annotation = ha1#,show_column_names = F
       )
dev.off()

## get the cluster annotation dataframe
r.dend<-row_dend(p)
rcl.list<-row_order(p)
lapply(rcl.list,function(x) length(x)) #check/confirm size of clusters
#loop to extract genes for each cluster
for(i in 1:length(rcl.list)){
    if(i==1){
        clu<-t(t(astr_mark[rcl.list[[i]]]))
        out<-cbind(clu,names(rcl.list)[i])
        colnames(out)<-c("gene","Cluster")
    } else {
        clu<-t(t(astr_mark[rcl.list[[i]]]))
        clu<-cbind(clu,names(rcl.list)[i])
        out<-rbind(out,clu)
    }
}
out<-as.data.frame(out)
rownames(out)<-out$gene
write.csv(out,"gfap.astr.317marker.scaled.express.binary.mean.in.time_ct2domain.heatmap.km4.geneclusters.csv")

c1<-out[out$Cluster==1,"gene"]
c2<-out[out$Cluster==2,"gene"]
c3<-out[out$Cluster==3,"gene"]
c4<-out[out$Cluster==4,"gene"]

###hclust and kmeans based on two genelist score dimension
set.seed(220609)
hc2<-cutree(hclust(dist(da@meta.data[,c('module33',
                                        'module11'
                                       )],
                        method = "maximum")),
            k = 2)
options(repr.plot.width=5,repr.plot.height=5)
ggplot(da@meta.data,aes(x=module33,y=module11,col=as.factor(hc2)))+
    #scale_color_manual(values = )+
    geom_point()+
    theme_classic()+
     scale_color_manual(values = c('#911eb4',"#ff7b25"))#+
    #geom_vline(xintercept = 0)+
    #geom_hline(yintercept = 0)
ggsave("sc.Astrocyte_Gfap.hclust2.bySpatialM13Score.scaterplot.png",width = 5,height = 5,dpi = 300)

da$hcluster<-as.factor(hc2)
Idents(da)<-da$hcluster
options(repr.plot.width=5,repr.plot.height=5)
VlnPlot(da,features = c("Gfap","Slc7a10"),cols = c('#911eb4',"#ff7b25"),slot = "data",pt.size = 0)
ggsave("sc.Astrocyte_Gfap.hclust2.bySpatialM13Score.Gfap_Slc7a10.data.Vlnplot.png",width = 5,height = 5,dpi = 300)

##marker and its bulk expression
test1<-da
Idents(test1)<-test1$hcluster
mark<-FindAllMarkers(test1,logfc.threshold = 0.5,only.pos = T,return.thresh = 0.01)
mark_s<-mark[mark$pct.1>0.5&mark$avg_log2FC>1,]
write.csv(mark_s,"gfap.317marker.k4.m13.hclust2.fc1.pct05.pad001.pos.markers.csv")
write.csv(mark,"gfap.317marker.k4.m13.hclust2.fc05.pad001.pos.markers.csv")
genes<-mark[,#mark$cluster==2,
            "gene"]

all_m<-unique(genes)
p_genes<-intersect(c(c1,c3),all_m)
p_genes
bulk_exp<-test1@assays$SCT@scale.data[all_m,]
df<-as.data.frame(t(bulk_exp))
df$cluster<-as.character(test1$hcluster)
df_mean<-aggregate(df[,colnames(df)!="cluster"],list(df$cluster),mean)

rownames(df_mean)<-df_mean[,1]
df_mean<-df_mean[,-1]
row_col<-data.frame(row.names = colnames(df_mean))
row_col$color<-"grey"
row_col[p_genes,"color"]<-"#ff7b25"#"#911eb4"

temp<-as.matrix(t(df_mean))
options(repr.plot.width=5,repr.plot.height=10)
p<-pheatmap::pheatmap(temp,show_rownames = T,cluster_cols = T,scale = "none",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "average",#labels_row = lab_row,
                   cellwidth = 20,
                   cellheight = 1,
                   fontsize = 8,
                   #na_col = NULL,
                   border_color = NA
                   #,scale = "none",col = viridis_pal()(50)
       )
cols<-row_col[order(match(rownames(temp),p$gtable$grobs[[5]]$label)),]
p$gtable$grobs[[5]]$gp=gpar(col=cols)
add.flag(p,kept.labels = p_genes,repel.degree = 1)

pdf("gfap.m13.hclust2.cluster.fc05.pad001.marker.bulk.pheatmap.pdf",width = 5,height = 10)
add.flag(p,kept.labels = p_genes,repel.degree = 1)
dev.off()

saveRDS(test1,"gfap.317_bi_m13.hclust2.cluster.rds")
write.csv(test1@meta.data,"gfap.317_bi_m13.hclust2.cluster.metadata.csv")

###FInd DEG between Gfap_WM, Gfap_GM and Slc7a10 astrocyte
da1<-readRDS("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.rds")
Idents(da1)<-da1$SCT_snn_res.0.35
da1$state<-"GM_Slc7a10"
da1$state[da1$gfap_hclust==1]<-"WM_Gfap"
da1$state[da1$gfap_hclust==2]<-"GM_Gfap"

###extract normalized data to do scenic analysis
meta<-da1@meta.data
write.csv(meta,"L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.metadata.csv")
norm_ma<-da1@assays$SCT@data
write.csv(norm_ma,"L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40.SCT.normalized.data.csv")
saveRDS(da1,"L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.rds")

###top enriched regulons of 3 states
regulon<-read.csv("lee_astr/regulons.csv")
regulonAUC<-read.csv("lee_astr/auc_mtx.csv")
rownames(regulonAUC)<-regulonAUC[,1]
regulonAUC<-regulonAUC[,-1]
#check all cells match
table(rownames(regulonAUC)==rownames(meta))
#split the cells by time
cellsPerTime<-split(rownames(meta),meta$state)
#regulonAUC2<-regulonAUC[cell,]
regulonAUC2<-regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
#Calculate average expression
regulonActivity_mean<-sapply(cellsPerTime,function(cells) colMeans(regulonAUC2[cells,]))

#scale expression
regulonActivity_mean_Scaled<-t(scale(t(regulonActivity_mean),center = T,scale = T))
regulonActivity_mean_Scaled_rm<-regulonActivity_mean_Scaled[complete.cases(regulonActivity_mean_Scaled),]
options(repr.plot.width=4, repr.plot.height=7) # To set the figure size in Jupyter
png("astr.3states.all.regulons.mean_scaled.activityscore.heatmap.png",width=4,height=7,units="in",res=300)
draw(ComplexHeatmap::Heatmap(regulonActivity_mean_Scaled_rm#[unique(top10$Regulon),]
                             , name="Regulon activity",
                       row_names_gp=grid::gpar(fontsize=0),cluster_columns=FALSE,cluster_rows=TRUE)) 
dev.off()
#To see the exact values
topRegulators <- reshape2::melt(regulonActivity_mean_Scaled_rm)
colnames(topRegulators) <- c("Regulon", "state", "RelativeActivity")
topRegulators$state <- factor(as.character(topRegulators$state))

topRegulators %>% 
    group_by(state) %>%
    top_n(.,n=10,wt=RelativeActivity) ->top10

options(repr.plot.width=4, repr.plot.height=8) # To set the figure size in Jupyter
png("astr.3states.top10.regulons.mean_scaled.activityscore.heatmap.png",width=4,height=8,units="in",res=300)
draw(ComplexHeatmap::Heatmap(regulonActivity_mean_Scaled_rm[unique(top10$Regulon),]
                             , name="Regulon activity",
                       row_names_gp=grid::gpar(fontsize=13),cluster_columns=FALSE,cluster_rows=FALSE)) 
dev.off()

### calculate RSS
rss<-calcRSS(AUC = t(regulonAUC),cellAnnotation = meta$state)
write.csv(rss,"Astrcoyte.3state.rrs.dataframe.csv")
rss_rm<-rss[complete.cases(rss),]
topRegulators <- reshape2::melt(rss_rm)
colnames(topRegulators) <- c("Regulon", "state", "RSS")
topRegulators$state <- factor(as.character(topRegulators$state))
topRegulators %>% 
    group_by(state) %>%
    top_n(.,n=10,wt=RSS) ->top10
temp<-regulonActivity_mean_Scaled_rm[unique(top10$Regulon),]
options(repr.plot.width=3, repr.plot.height=7) # To set the figure size in Jupyter
pdf("Astrocyte.3states.RSS.top10_regulons.mean.scaled.activityscore.heatmap_220728.pdf",width=3,height=7)
pheatmap(regulonActivity_mean_Scaled_rm[unique(top10$Regulon),], 
         name="Activity",
         main="Regulon activity",
         border_color="NA",
         cluster_cols=FALSE,
         cluster_rows=FALSE)
#save_pheatmap_pdf(p@layout,filename = "Astrocyte.3states.RSS.top10_regulons.mean.scaled.activityscore.heatmap_220728.pdf",width=3,height=7)
dev.off()




###DEG of each state
#GM_Gfap vs. WM_Gfap
deg1<-FindMarkers(da1,logfc.threshold = 0.5,ident.1="GM_Gfap",ident.2="WM_Gfap")
deg1<-subset(deg1,p_val_adj<0.01)
deg1$gene<-rownames(deg1)
deg1$group<-ifelse(deg1$avg_log2FC>0,"GM_Gfap","WM_Gfap")
deg1$group<-factor(deg1$group,levels = c("GM_Gfap","WM_Gfap"))
deg1<-deg1[order(deg1$group),]
deg1<-deg1[order(deg1$avg_log2FC,decreasing = TRUE),]
write.csv(deg1,"GM_Gfap.vs.WM_Gfap.fc05.pad001.DEG.csv")

#GM_Gfap vs. GM_Slc7a10
deg2<-FindMarkers(da1,logfc.threshold = 0.5,ident.1="GM_Gfap",ident.2="GM_Slc7a10")
deg2<-subset(deg2,p_val_adj<0.01)
deg2$gene<-rownames(deg2)
deg2$group<-ifelse(deg2$avg_log2FC>0,"GM_Gfap","GM_Slc7a10")
deg2$group<-factor(deg2$group,levels = c("GM_Gfap","GM_Slc7a10"))
deg2<-deg2[order(deg2$group),]
deg2<-deg2[order(deg2$avg_log2FC,decreasing = TRUE),]
write.csv(deg2,"GM_Gfap.vs.GM_Slc7a10.fc05.pad001.DEG.csv")

### WM_GFap vs. GM_Slc7a10
deg3<-FindMarkers(da1,logfc.threshold = 0.5,ident.1="WM_Gfap",ident.2="GM_Slc7a10")
deg3<-subset(deg3,p_val_adj<0.01)
deg3$gene<-rownames(deg3)
deg3$group<-ifelse(deg3$avg_log2FC>0,"WM_Gfap","GM_Slc7a10")
deg3$group<-factor(deg3$group,levels = c("WM_Gfap","GM_Slc7a10"))
deg3<-deg3[order(deg3$group),]
deg3<-deg3[order(deg3$avg_log2FC,decreasing = TRUE),]

deg1_p<-rownames(deg1)[deg1$avg_log2FC>0]
deg1_n<-rownames(deg1)[deg1$avg_log2FC<0]

deg2_p<-rownames(deg2)[deg2$avg_log2FC>0]
deg2_n<-rownames(deg2)[deg2$avg_log2FC<0]

deg3_p<-rownames(deg3)[deg3$avg_log2FC>0]
deg3_n<-rownames(deg3)[deg3$avg_log2FC<0]

deg_dp<-intersect(deg1_p,deg2_p)
write.table(deg_dp,"hclust.GM_Gfap.specific.marker.fc05.pad001.txt")
deg_wm<-intersect(deg1_n,deg3_p)
deg_gm<-intersect(deg2_n,deg3_n)


alldeg<-rbind(deg1,deg2,deg3)
alldeg$group<-as.character(alldeg$group)
alldeg[alldeg$gene %in% union(deg1_p,deg2_p),"group"]<-"GM_Gfap.specific"
alldeg[alldeg$gene %in% union(deg3_p,deg1_n),"group"]<-"WM_Gfap.specific"
alldeg[alldeg$gene %in% union(deg3_n,deg2_n),"group"]<-"GM_Slc7a10.specific"
alldeg[alldeg$gene %in% dp_wm_deg,"group"]<-"GM_Gfap.WM_Gfap"
alldeg[alldeg$gene %in% dp_gm_deg,"group"]<-"GM_Gfap.GM_Slc7a10"
alldeg<-alldeg[!duplicated(alldeg$gene),]
alldeg$group<-factor(alldeg$group,
                     levels = c("WM_Gfap.specific","GM_Gfap.WM_Gfap",
                                             "GM_Slc7a10.specific",
                                             "GM_Gfap.GM_Slc7a10","GM_Gfap.specific"))

alldeg<-alldeg[order(alldeg$group,decreasing = FALSE),]
write.csv(alldeg,"3state.alldeg.5groups.csv")




###Display the surface marker genes of Astro-GMii
sf_gene<-read_excel("/S2_File.xlsx",sheet = "Table B")
# dp state surface markers
cand_sg<-intersect(deg_dp,sf_gene$`ENTREZ gene symbol`)
pdf("Astr_GMii.bothpositive.surface.genes.Vlnplot.220818.pdf",width = 8,height = 8)
options(repr.plot.width=8,repr.plot.height=8)
VlnPlot(da1,features = cand_sg,ncol = 5,
        cols = c('#911eb4',"#ff7b25","#e6beff"),
        #slot = "scale.data",
        pt.size = 0)
dev.off()
# WM_Gfap surface markers
cand_sg2<-intersect(deg_wm,sf_gene$`ENTREZ gene symbol`)
# GM_Slc7a10 surface markers
cand_sg3<-intersect(deg_gm,sf_gene$`ENTREZ gene symbol`)


###Overlap of SA, scGMii_DEG and trajectory genes. Check the correlation with GMii state spots
decov<-readRDS("WT.merge.SCT_RCTD.c5.rds")
deg<-mark_s[mark_s$cluster==2,'gene']
deg<-intersect(deg,rownames(decov@assays$SCT@data))
traj<-read.csv("slingshot.dynamic.gene.logfc1.5.csv")
traj_g<-intersect(traj$X,rownames(decov@assays$SCT@data))
cross_gene<-Reduce(intersect,list(deg,c(c1,c3),traj_g))

venn_list<-list(deg,c(c1,c3),traj_g)
names(venn_list)<-c('sc_GMii.vs.WMrs_markers','spatial SA markers','Trajectory genes')

pdf("GMii.spatial.sc.marker.overlap.pdf")
venn(venn_list)
dev.off()

# 1. get the deconv proportion within GM region
gmii_decov<-decov$Astrocyte_Gfap
gmii_decov[decov$time %in% c('WT_sham','WT_3h')]<-0
gmii_decov[!decov$domain %in% c('DH','MG','VH')]<-0
gmii_decov[gmii_decov>0]<-1
expr<-decov@assays$SCT@scale.data[cross_gene,]
SA_bi<-apply(expr,1,function(x) ifelse(x<(mean(x)+sd(x)),0,1))

# define function
jaccard<-function(a,b){
    a<-names(a)[a==1]
    b<-names(b)[b==1]
    int<-length(intersect(a,b))
    union<-length(a)+length(b)-int
    return(int/union)
}


simpson<-function(a,b){
    a<-names(a)[a==1]
    b<-names(b)[b==1]
    int<-length(intersect(a,b))
    mins<-min(length(a),length(b))
    return(int/mins)
}

test<-apply(SA_bi,2,function(x) jaccard(x, gmii_decov))
test[order(test,decreasing = TRUE)]


###Neuron loss, SA_score and injury_score
da<-readRDS("gfap.317_bi_m13.hclust2.cluster.rds")
injury_genes<-c('Adamts1','Atf3','Ccl2','Ccnd1','Cd68','Cebpd','Cyba',
                'Fn1','Gal','Gap43','Hmox1','Hspb1','Igfbp2','Jun','Junb',
                'Fos','Lgals1','Neat1','Socs3','Tnc','S100a10','Timp1')
da<-AddModuleScore(da,features = list(injury_genes),name = "injury_score")
me<-da@meta.data
neuron<-read.csv("15celltypes.dec_conf_nor.meta.csv")
rownames(neuron)<-neuron[,1]
neuron<-neuron[,-1]

meta<-cbind(me[,c('orig.ident','distance',"time","RC","pos","domain_res02","x.ad","y.ad","injury_score1",
                  "astr_marker_module13_score1","Gfap","Slc7a10","Pla2g7","Igfbp2")],neuron[,"Neuron"])
colnames(meta)[ncol(meta)]<-"Neuron"
write.csv(meta,"Neuronloss_SAscore_Injuryscore.meta.csv")
meta$pos<-factor(meta$pos,levels = unique(meta$pos))
#subset 24hpi and 72hpi
meta_sub<-meta[meta$time %in% c('WT_24h','WT_72h'),]
meta_sub<-meta_sub[order(meta_sub$pos,decreasing = TRUE),]
me3<-meta_sub[meta_sub$time=="WT_72h",]
me4<-me3[me3$domain_res02 %in% c("MG","DH","VH"),]
me4$dis2<-ifelse(me4$distance=="1mm",1,0.5)
spc<-cor(me4[,c("injury_score1","astr_marker_module13_score1",'Neuron',"dis2")],method = "spearman")
spc<-melt(spc)
# normalize
me4$injury_score1<-(me4$injury_score1-min(me4$injury_score1))/(max(me4$injury_score1)-min(me4$injury_score1))
me4$astr_marker_module13_score1<-(me4$astr_marker_module13_score1-min(me4$astr_marker_module13_score1))/(max(me4$astr_marker_module13_score1)-min(me4$astr_marker_module13_score1))
me4$spot<-seq(1,nrow(me4))
temp<-melt(me4[,c("spot","Neuron",#"injury_score1",
                  "astr_marker_module13_score1","pos")],
           id.vars = c("spot","pos"),variable.name = "group")
pdf("Neuronloss_SAscore.loessfitting.72h.GM.position_221024.pdf",height = 5,width = 9)
options(repr.plot.width=9,repr.plot.height=5)
ggplot(temp, aes(x = spot, y = value,color=group)) +
    geom_point(size=0.01) +
    scale_color_manual(values = c('#034f84',
                                  #'#5D3B83',
                                  '#C580AA',
                                  '#316B9B'#,'#E2E293'
                                 ))+
    stat_smooth(method = "loess",se=FALSE,size=2#,
        #col = "#C42126",
        #se = FALSE,
        #size = 1
               )+
    theme_classic()
dev.off()


###bulk expression heatmap and Go terms of all deg across three spatial astrocyte states
bulk_exp<-as.matrix(da1@assays$SCT@data[alldeg$gene,])
df<-as.data.frame(t(bulk_exp))
df$state<-as.character(da1$state)
df_mean<-aggregate(df[,colnames(df)!="state"],list(df$state),mean)
rownames(df_mean)<-df_mean[,1]
df_mean<-df_mean[,-1]
df_mean<-df_mean[c(3,2,1),]
write.csv(df_mean,"3states.5groups.deg.bulk.data.csv")
col_ha<-as.data.frame(alldeg$group)
colnames(col_ha)<-"group"
rownames(col_ha)<-alldeg$gene

## go term
library(org.Mm.eg.db)
library(ggplot2)
#grlabs<-split(gene.nodle.df$gene,gene.nodle.df$module)
grlabs<-split(alldeg$gene,alldeg$group)
gcSample = lapply(grlabs, function(gr) as.numeric(bitr(gr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID))
pvalueCutoff = 0.05
qvalueCutoff = 0.05
xx.mus.go <- compareCluster(gcSample, OrgDb='org.Mm.eg.db', fun='enrichGO', 
                            pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ont = "BP", 
                            readable=T,universe=backgenes,maxGSSize=1000
                           )
saveRDS(xx.mus.go,"3state.fc05.pad005.alldeg.5group.pqval001.max1000.go.rds")
df<-xx.mus.go@compareClusterResult
df1 = df %>% group_by(.,Cluster) %>% filter(.,Count>=5) %>% top_n(.,-20,p.adjust) 
df1$Description = factor(df1$Description,levels = rev(unique(df1$Description))) 
options(repr.plot.width=10, repr.plot.height=18)
library(ggplot2)
p=ggplot(df1,aes(Cluster,Description))+geom_point(aes(size=Count,color=p.adjust))+        #基因数目表示SIZE，p.adjust代表颜色深浅
    scale_color_gradientn(colors=brewer.pal(9,"Spectral"),guide = guide_colorbar(reverse = T))+
    theme_bw()+
    theme(axis.title = element_blank(),axis.text.y.left = element_text(size=13),
          axis.text.x.bottom = element_text(size = 13,angle = 45,hjust = 1,vjust = 1))
p
ggsave("3state.fc05.pad001.alldeg.5group.pqval005.max1000.go.C5.pad_top20.point.png",width = 10,height = 18,dpi = 300)

##correlate to other studies and add associated gene into heatmap
## proflammatory-factor
infm<-c("Angpt2","Tek","Ccl2","Ackr1","Cd34","Sell","Cd44","Sele","Csf1","Csf1r","Glg1","Selplg",
        "Tnfsf10","Tnfrsf10b","Alox5ap","Alox5","C3","C3ar1","Ccl5","Esam","C3ar1","Ccr1",
        "Hebp1","Lgals9","Havcr2")

infm_s<-infm[infm %in% alldeg$gene]

### anti-infalmmatory genes
anti_inf<-c("Sirpa","Tyro3","Cx3cl1","Cx3cr1","Il6","Gas6","Anxa1","Pros1","Pomc",
            "Tgfb2","Il6r","Il6st","Tyro3","Mc1r","Jag1","Grn","Cd47",
            "Cd46","Mc1r","Axl","Tgfbr1","Tgfbr2",
            "Itgax","Itgb2","Fpr1","Tnfrsf1b","Tnfrsf11b")
anti_c<-alldeg[alldeg$gene %in% anti_inf,"gene"]
### trophic factors
tr_f<-c("Bmp","Bmpr1b","Bmpr2","Enho","Gpr19","Fgf1","Fgfr3","Fgf9","Nectin1","Psap","Bmp4",
                         "Gpr37l1","Flt1","Vegfa","Nrp1","Nrp2","Fgf7","Fgfr1","Gdf11","Tgfbr1",
                         "Acvr2b","Hgf","Met","Pdgfb","Pdgfra","Pdgfc","Pdgfd","Rspo3","Lgr5","Wnt2","Pdgfrb",
                         "Fzd4","Wnt5a","Fzd1","Fzd5")

tr_f_s<-tr_f[tr_f %in% alldeg$gene]
prot_c<-c("Csf1","Pros1",tr_f_s)
cand_df<-alldeg[alldeg$gene %in% c(#pig,
                                   anti_c,tr_f_s,infm_s),]
cand_df$label<-""
#cand_df$label[cand_df$gene %in% pig]<-"PIGs"
cand_df$label[cand_df$gene %in% anti_c]<-"Anti-inflammatory"
cand_df$label[cand_df$gene %in% tr_f_s]<-"Trophic.factor"
cand_df$label[cand_df$gene %in% infm_s]<-"Pro-inflammatory"
cand_df$protective<-ifelse(cand_df$gene %in% prot_c,"y","n")
# add state bulk normalized expression
bulk_exp<-as.matrix(da1@assays$SCT@data[cand_df$gene,])
df<-as.data.frame(t(bulk_exp))
df$state<-as.character(da1$state)
df_mean<-aggregate(df[,colnames(df)!="state"],list(df$state),mean)
rownames(df_mean)<-df_mean[,1]
df_mean<-scale(df_mean[,-1])
cand_df<-cbind(cand_df,t(df_mean))
df1<-cand_df[,c("gene","group","label","protective","GM_Gfap","GM_Slc7a10","WM_Gfap")]
df1_l<-melt(df1)
df1_l$label<-factor(df1_l$label,levels = c('Trophic.factor','Anti-inflammatory','Pro-inflammatory'#,'PIGs'
                                          ))
df1_l<-df1_l[order(df1_l$label),]
write.csv(df1_l,"associate.gene_spatial.state.average.expression.df_220630.csv")

temp<-as.matrix(t(df_mean))
col_ha<-as.data.frame(alldeg$group)
colnames(col_ha)<-"group"
rownames(col_ha)<-alldeg$gene
row_col<-data.frame(row.names = colnames(df_mean))
row_col$color<-"grey"
row_col[unique(df1_l$gene),"color"]<-"#ff7b25"
group<-c('#ff7b25','#eea29a','#911eb4','#d5f4e6','#e6beff')
names(group)<-c("GM_Gfap.specific","GM_Gfap.WM_Gfap","WM_Gfap.specific","GM_Gfap.GM_Slc7a10","GM_Slc7a10.specific")
an_col<-list(group=group)

options(repr.plot.width=6,repr.plot.height=12)
p<-pheatmap::pheatmap(temp,show_rownames = T,show_colnames = T,
                      annotation_colors = an_col,
                       annotation_row = col_ha,
                      cluster_cols = F,
                      cluster_rows = F,scale = "row",
                      clustering_distance_rows = "euclidean",
                      clustering_method = "average",#labels_row = lab_row,
                   cellwidth = 20,
                   cellheight = 0.5,
                   fontsize = 8,
                   #na_col = NULL,
                   border_color = NA
                   #,scale = "none",col = viridis_pal()(50)
       )
cols<-row_col[order(match(rownames(temp),p$gtable$grobs[[5]]$label)),]
p$gtable$grobs[[5]]$gp=gpar(col=cols)
add.flag(p,kept.labels = unique(df1_l$gene),repel.degree = 1)
pdf("3state.alldeg.5groups.fc05.pad001.marker.bulk.data.mean.scale_associategenes.heatmap_220630.pdf",width = 5,height = 10)
add.flag(p,kept.labels = unique(df1_l$gene),repel.degree = 1)
dev.off()




