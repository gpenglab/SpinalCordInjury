library(heatmaply)
library(scFunctions)
library(tidyverse)
library(ggridges)
library(pheatmap)
library(RColorBrewer)
library(dendextend)
library(reshape2)
library(plotly)
library(htmlwidgets)
source("CSI.R")
source("script/self_function/save_pheatmap_pdf.R")

binary_auc<-read.csv("../../WT.SCT.binarizedbyPyscenic.auc_mtx.csv")
rownames(binary_auc)<-binary_auc$Cell
binary_auc<-binary_auc[,-1]
#transposition the matrix
binary_regulons_trans<-as.matrix(t(binary_auc))
regulons_trans<-as.matrix(t(auc_mtx))

###Calculate regulon specificity score (RSS)
meta<-read.csv("WT_replace_v2/res02_220310/WT.SCT.pc20.k50.res02.meta.data.csv")
rownames(meta)<-meta[,1]
meta<-meta[,-1]
meta$time_domain<-paste0(meta$time,"_",meta$domain_res02)
meta$cell_type<-meta$time_domain
table(colnames(binary_regulons_trans)==rownames(meta))
# calculate the RSS for all regulons over all spatial clusters
rrs_df<-calculate_rrs(meta,binary_regulons = binary_regulons_trans,cell_type_column = "cell_type")
rss<-melt(rrs_df)
colnames(rss)<-c("regulon","cell_type","RSS")
write.csv(rss,"WT.regulons.res02_time-domain.rss.dataframe.csv")
write.csv(rrs_df,"WT.binary.regulons.res02_time-domain.rrs.dataframe.csv")
rrs_df %>% 
    group_by(cell_type) %>%
    top_n(.,n=20,wt=RSS) ->rss_sub
subR<-unique(rss_sub$regulon)
# visualize all regulons over all cell types using heatmaps.
rrs_df_wide<-rrs_df %>% spread(cell_type,RSS)
rownames(rrs_df_wide)<-rrs_df_wide$regulon
rrs_df_wide<-rrs_df_wide[,2:ncol(rrs_df_wide)]
## subset all regulons that don't have at least an RSS of 0.2 for one cell type
rrs_df_wide_specific<-rrs_df_wide[subR,]
# order col
rrs_df_wide_specific<-rrs_df_wide_specific[,c('WT_sham_DH','WT_3h_DH','WT_24h_DH','WT_72h_DH',
                                              'WT_sham_MG','WT_3h_MG','WT_24h_MG','WT_72h_MG',
                                              'WT_sham_VH','WT_3h_VH','WT_24h_VH','WT_72h_VH',
                                              'WT_sham_WM','WT_3h_WM','WT_24h_WM','WT_72h_WM')
                                          ]
col_anno<-data.frame(row.names = colnames(rrs_df_wide_specific),"domain"=c(rep("WM",4),rep("MG",4),rep("DH",4),rep("VH",4)),
                    "time"=c(rep(c("sham","3h","24h","72h"),4)))
domain_col<-c('#20854EA8','#0072B5A8','#BC3C29A8','#E18727A8')
names(domain_col)<-c("WM", "MG", "DH", "VH")
time_col<-c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF')
names(time_col)<-c("sham","3h","24h","72h")
col_list<-list(domain=domain_col,time=time_col)
gap_col<-c(4,8,12)
d_hc<-hclust(dist(rrs_df_wide_specific))
den<-as.dendrogram(hclust(dist(rrs_df_wide_specific)))
options(repr.plot.width=6,repr.plot.height=7)
#pdf("WT.regulon.RSS.thrs01.cluster.heatmaply.pdf",width=8,height=8)
p<-pheatmap(rrs_df_wide_specific,col=colorRampPalette(rev(brewer.pal(10,"Spectral")))(50),scale="none",
          cluster_cols=F,cluster_rows=T,treeheight_row=10,cutree_rows=7,annotation_col=col_anno,border_color="NA",
         annotation_colors=col_list,legend=T,gaps_col=gap_col,
            show_rownames=T#,legend_labels=c("RSS")
          #file = "WT.SCT.pct01.data.regulons.01wide.RSS.of.clusters.pheatmap.html",
          ) 
save_pheatmap_pdf(p,"WT.res02.time-domain.binary.regulon.RSS.top3.pheatmap.pdf",width = 6,height = 7)

### Calculate connection specificity index (CSI) for all regulons/selected top regulons
regulonAUC<-t(binary_auc)
regulonAUC_sub<-regulonAUC[subR,]
regulons_csi <- calculate_csi(regulonAUC_sub,
                              calc_extended = FALSE)
test_sub<-t(test)
pearson_cor <- cor(test)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    gather(regulon_2,pcc,-regulon_1) %>%
    mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))

pcc$PN<-NA
pcc[pcc$pcc<0,"PN"]<-"0"
pcc[pcc$pcc>0,"PN"]<-"1"
pcc$PN<-as.factor(pcc$PN)
pearson_cor_long<-merge(pcc,regulons_csi[,c("regulon_pair","CSI")],by = "regulon_pair")
write.csv(pearson_cor_long,"time_domain.RSS.top20.regulons.genescaleexp_correlation.and.csi.df.csv")
pearson_cor_long<-pearson_cor_long[abs(pearson_cor_long$pcc)>0.65#&pearson_cor_long$CSI>0
                                   ,]
pearson_cor_long<-pearson_cor_long[!pearson_cor_long$regulon_1==pearson_cor_long$regulon_2,]
pcc_df_filt<-pearson_cor_long
cluster_df_sub<-clusters_df[clusters_df$regulon %in% unique(c(pcc_df_filt$regulon_1,pcc_df_filt$regulon_2)),]
reg<-unique(c(pcc_df_filt$regulon_1,pcc_df_filt$regulon_2))
write.csv(pcc_df_filt,"rss.top20.auc.regulon.auc_pearson.correlation065.long.CSI.csv",quote=FALSE)
write.csv(cluster_df_sub,"rss.top20.auc.regulon.auc_pearson.correlation065.long.CSI.nodes.csv",quote=FALSE)
write.csv(pearson_cor_long,"rss.top20.regulon.pearson.correlation.long.csv",quote=FALSE)
write.csv(regulons_csi,"WT.SCT.binary.regulons.SCT.res02.time-domain.RSS.top20.csi.csv")
regulons_csi$regulon_1<-gsub("...$","",regulons_csi$regulon_1)
regulons_csi$regulon_2<-gsub("...$","",regulons_csi$regulon_2)
regulons_csi$regulon_pair<-paste0(regulons_csi$regulon_1,"_",regulons_csi$regulon_2)

# finnally, export the SCI cluster modules by the following code with the same cluster number to create the heatmap
csi_csi_wide<- regulons_csi %>% 
    spread(regulon_2,CSI)

future_rownames<-csi_csi_wide$regulon_1
csi_csi_wide<-as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
rownames(csi_csi_wide)<-future_rownames
write.csv(csi_csi_wide,"WT.SCT.res02.time-domain.RSS_top20.binary_csi.wide.matrix.csv")
subR<-rownames(csi_csi_wide)
regulons_hclust<-hclust(dist(csi_csi_wide,method="euclidean"),method = "complete")
clusters<-cutree(regulons_hclust,k=7)
clusters_df<-data.frame("regulon"=names(clusters),
                       "csi_cluster"=clusters)
anno_row<-data.frame(row.names = rownames(clusters_df),"CSI_module"=as.character(clusters_df[,2]))
anno_col<-brewer.pal(7,"Set3")
names(anno_col)<-unique(anno_row$CSI_module)
anno_col<-list(CSI_module=anno_col)
p<- pheatmap(csi_csi_wide,annotation_row = anno_row,annotation_colors = anno_col,
           show_colnames = FALSE,
           show_rownames=FALSE,
           display_numbers=F,
           color = viridis(n = 25),
           cutree_cols = 7,
           cutree_rows = 7,
           #fontsize_row = 0,#font_size_regulons,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           treeheight_row = 5,
           treeheight_col = 5,
           clustering_method="complete",
           clustering_distance_ro0ws = "euclidean",
           clustering_distance_cols = "euclidean",
           width = 300,
           height = 800)
  #dev.off()
save_pheatmap_pdf(p,"WT.res02_time_domain.rss_top20.csi.complete_nclust7.hcluster.withAnno.pheatmap.pdf")

## export module node dataframe to cytoscape for regulatory network visualization
rownames(clusters_df)<-gsub("...$","",rownames(clusters_df))
clusters_df$regulon<-gsub("...$","",clusters_df$regulon)
write.csv(clusters_df,"rss.top20.regulons.csi.7hcluster.df.csv",quote=FALSE)

## Finally, calculate activity scores for each csi module based on the specifity scores of all the regulons in that module.
csi_cluster_activity_wide<-calc_csi_module_activity(clusters_df,regulonAUC_sub,metadata = meta,cell_type_column = "cell_type")
csi_cluster_activity_wide<-csi_cluster_activity_wide[,c('WT_sham_WM','WT_3h_WM','WT_24h_WM','WT_72h_WM',
                                                       'WT_sham_MG','WT_3h_MG','WT_24h_MG','WT_72h_MG',
                                                       'WT_sham_DH','WT_3h_DH','WT_24h_DH',"WT_72h_DH",
                                              'WT_sham_VH','WT_3h_VH','WT_24h_VH','WT_72h_VH'
                                              )]
csi_long<-melt(t(scale(t(csi_cluster_activity_wide),center = T,scale = T)),value.name = "csi",varnames = "cluster_csi")#,variable.name = "cluster_exp"
colnames(csi_long)<-c("cluster_csi","cluster_exp","csi")
top3<-csi_long %>% group_by(cluster_exp) %>% top_n(.,wt = csi,n = 1)
top3$cluster_csi = factor(top3$cluster_csi,levels = rev(unique(top3$cluster_csi)))
csi_cluster_wide<-csi_cluster_activity_wide[c(2,3,1,6,4,5,7),
                                           ]
col_anno<-data.frame(row.names = colnames(csi_cluster_wide),"domain"=c(rep("WM",4),rep("MG",4),rep("DH",4),rep("VH",4)),
                    "time"=c(rep(c("sham","3h","24h","72h"),4)))
anno_row<-data.frame(row.names = unique(clusters_df$csi_cluster),"CSI_module"=as.character(unique(clusters_df$csi_cluster)))
anno_col<-brewer.pal(7,"Set3")
names(anno_col)<-anno_row$CSI_module
col_list<-list(CSI_module=anno_col,domain=domain_col,time=time_col)  
options(repr.plot.width=6,repr.plot.height=5)
p<-pheatmap(csi_cluster_wide,border_color = "NA",
           show_colnames = TRUE,
           color = viridis(n = 100),#colorRampPalette(rev(brewer.pal(6,"RdYlBu")))(50),#viridis(n = 100),
           cluster_cols = F,
           cluster_rows = F,
           annotation_col=col_anno,
           annotation_colors=col_list,
            legend=T,
            gaps_col=gap_col,
            annotation_row=anno_row,
            
           #scale="row",
        
           #clustering_distance_rows = "euclidean"#,
           #clustering_distance_cols = "euclidean"
            )
save_pheatmap_pdf(p,"WT.SCT.binary.regulon.res02_time-domain.rss_top20.7csi.heatmap.viridis_220720.pdf")

### CSI module average activity in space
auc_mtx<-read.csv("/home/jovyan/zxli_SCI/result/pyscenic/WT.merge.replace_v2/SCT/auc_mtx.csv")
rownames(auc_mtx)<-auc_mtx[,1]
auc_mtx<-auc_mtx[,-1]
colnames(auc_mtx)<-gsub("...$","",colnames(auc_mtx))
meta<-cbind(meta,auc_mtx)
temp<-auc_mtx[,clusters_df$regulon]
temp<-as.data.frame(t(temp))
temp$module<-clusters_df$csi_cluster
temp_m<-temp
temp_m$module<-as.character(temp_m$module)
temp_m<-aggregate(temp_m[,colnames(temp_m)!="module"],list(temp_m$module),mean)
rownames(temp_m)<-temp_m[,1]
temp_m<-temp_m[,-1]
write.csv(temp_m,"7csi.module.average.auc.activity.eachSpot.df.csv")
rownames(temp_m)<-paste0("module",rownames(temp_m))
meta_m<-cbind(meta,t(temp_m))
saveRDS(meta_m,"7csi.module.average.auc.activity.eachSpot.meta.rds")
sf<-"module.spatial_220720/"
if(!dir.exists(sf))
    dir.create(sf)
for(i in paste0("module",seq(1:7))){
    options(repr.plot.width=20,repr.plot.height=16)
    ggplot(meta_m,aes_string("x.ad","y.ad",color=i))+geom_point(size=1)+
        scale_color_gradientn(colours = c('#0000ff','#3232ff','#6666ff','#9898ff','#ccccff',
                                              '#fffefe','#ffcccc','#ff9898','#ff6666','#ff3232','#ff0000'))+
        #scale_fill_gradientn(colours = colorRampPalette(viridis(n = 10,option = "B"))(50))+
        xlab(paste0("")) +
        ylab(paste0("")) + 
        theme(panel.background = element_blank(),
              panel.grid.major =element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
              axis.text.y = element_blank()#element_text(size = 20,colour = "black")
              # axis.ticks.x = element_blank(),
        )
    ggsave(paste0(sf,"WT.merge.regress_CC.nC.mt.ident_",i,".spatial_manually_220720.png"),width=20,height=16,dpi=300)
}



