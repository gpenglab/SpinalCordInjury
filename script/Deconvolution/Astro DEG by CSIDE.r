library(dplyr)
library(spacexr)
library(stringr)
library(clusterProfiler)
options(connectionObserver = NULL)
library(biomaRt)
library(Matrix)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(data.table)
library(tidyverse)
library(pheatmap)
library(scales)
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

###Load reference
reference<-readRDS("Reference.renamed.merged13celltypes.220507.rds")
meta<-read.csv("../16celltypes.dec_conf.nor.meta.csv")
rownames(meta)<-meta[,1]
meta<-meta[,-1]
celltypes<-c('Astro.Svep1','Astrocyte.Gfap','Astrocyte.Slc7a10',
             'Dendritic','Fibroblast','Homeostatic.Microglia',
             'Monocyte','Neuron','Neutrophil','Oligodendrocyte',
             'Vascular','OPC','Macrophage','R.Microglia',
             'Myeloid','Ependymal')
weights_df<-meta[,celltypes]

###check astrocyte subtypes spatial distribution specificity
wm_cells<-rownames(subset(meta,domain=="WM"&orig.ident=="WT_sham_H_R2"))
gm_cells<-rownames(subset(meta,domain%in%c("MG","DH","VH")&orig.ident=="WT_sham_H_R2"))
mg_cells<-rownames(subset(meta,domain%in%c("MG")&orig.ident=="WT_sham_H_R2"))
wdv_cells<-rownames(subset(meta,domain%in%c("WM","DH","VH")&orig.ident=="WT_sham_H_R2"))

gfap_cells<-rownames(subset(meta,Astrocyte.Gfap>0&orig.ident=="WT_sham_H_R2"))
slc_cells<-rownames(subset(meta,Astrocyte.Slc7a10>0&orig.ident=="WT_sham_H_R2"))

wm_gfap_cells<-intersect(wm_cells,gfap_cells)
gm_slc_cells<-intersect(gm_cells,slc_cells)

sub_cells<-unique(c(wm_gfap_cells,gm_slc_cells))

##merge into one celltype
weights_ad<-weights_df
weights_ad$Astrocyte<-rowSums(weights_ad[,c("Astro.Svep1","Astrocyte.Gfap","Astrocyte.Slc7a10")])
weights_ad<-weights_ad[,-which(colnames(weights_ad)%in%c("Astro.Svep1","Astrocyte.Gfap","Astrocyte.Slc7a10"))]
for(i in unique(meta$sample)){
    cells<-meta[meta$sample==i,"X"]
    temp_wt<-weights_ad[cells,]
    write.csv(temp_wt,paste0(datadir,i,".conf_nor.13cp.weights.df.csv"))
}

###load sham_H_1mm 4 replicates
puck1<-readRDS("WT_sham_H_R2_1.puck.rds"
                     )

puck2<-readRDS("WT_sham_H_R2_2.puck.rds"
                     )

puck3<-readRDS("WT_sham_H_R2_3.puck.rds"
                     )

puck4<-readRDS("WT_sham_H_R2_4.puck.rds"
                     )

###Create CSIDE object
spatial.replicates <- list(puck1, puck2, puck3, puck4)
replicate_names <- c('1','2','3','4')
group_ids <- c(1,1,1,1)
myRCTD.reps <- create.RCTD.replicates(spatial.replicates,
                                      reference, 
                                      replicate_names, 
                                      group_ids = group_ids,
                                      UMI_min=0,
                                      UMI_max=2000000,
                                      CELL_MIN_INSTANCE=10,
                                      keep_reference = T)

###get 13 celltypes' deconvolution markers
sc_df<-get_marker_data(myRCTD.reps@RCTD.reps[[1]]@cell_type_info$info[[2]],myRCTD.reps@RCTD.reps[[1]]@cell_type_info$info[[1]],myRCTD.reps@RCTD.reps[[1]]@internal_vars$gene_list_bulk)
sc_df$gene<-rownames(sc_df)
#marker expression matrix
sc_exp<-myRCTD.reps@RCTD.reps[[1]]@cell_type_info$info[[1]]
sc_exp_norm<-sc_exp[sc_df$gene,]
sc_exp_norm<-sc_exp_norm / rowSums(sc_exp_norm)
sc_exp_norm$gene<-rownames(sc_exp_norm)

sc_df<-cbind(sc_df,sc_exp_norm[sc_df$gene,])
write.csv(sc_df,"Astr_merged.total13celltypes.normalized.df.csv")

### visualize Astrocyte markers in deconvolution reference matrix
p_df<-sc_exp_norm[sc_df$gene[sc_df$cell_type=="Astrocyte"],]
options(repr.plot.width=5,repr.plot.height=14)
pdf("Astrocyte.markers.express.in.10Xsc.celltype.marker.heatmap.pdf",width = 5,height = 14)
Heatmap(p_df[,-which(colnames(p_df)=="gene")],show_row_names = F,
        cluster_rows = F ,col=c("#05445E",#"#54BAB9",
                                                                      #"white",
                                                                      "#F6E6E8",
                                                                      #"#E9DAC1",
                                                                      "#E56997",#"#F51720",
                                                                      "#B91646"))
dev.off()

cells1<-rownames(puck1@coords)
cells2<-rownames(puck2@coords)
cells3<-rownames(puck3@coords)
cells4<-rownames(puck4@coords)

weight_df1<-weights_ad[cells1,]
weight_df2<-weights_ad[cells2,]
weight_df3<-weights_ad[cells3,]
weight_df4<-weights_ad[cells4,]

myRCTD.reps@RCTD.reps[[1]]<-import_weights(myRCTD.reps@RCTD.reps[[1]],weight_df1)
myRCTD.reps@RCTD.reps[[2]]<-import_weights(myRCTD.reps@RCTD.reps[[2]],weight_df2)
myRCTD.reps@RCTD.reps[[3]]<-import_weights(myRCTD.reps@RCTD.reps[[3]],weight_df3)
myRCTD.reps@RCTD.reps[[4]]<-import_weights(myRCTD.reps@RCTD.reps[[4]],weight_df4)

saveRDS(myRCTD.reps,"13cp.sham_H.myRCTD.reps.rds")

###set variable that deviate two regional astrocyte
exp.var1<-as.integer(cells1%in%gfap_cells)
names(exp.var1)<-cells1
exp.var2<-as.integer(cells2%in%gfap_cells)
names(exp.var2)<-cells2
exp.var3<-as.integer(cells3%in%gfap_cells)
names(exp.var3)<-cells3
exp.var4<-as.integer(cells4%in%gfap_cells)
names(exp.var4)<-cells4

exvar_list <- list(exp.var1, exp.var2, exp.var3, exp.var4)
saveRDS(exvar_list,"sham_H.Astr.subtype.exvar.list.rds")

##only retain major cell types at uninjury tissue
cell_types<-c('Astrocyte','Neuron','Oligodendrocyte',"Homeostatic.Microglia","Vascular")

###run CSIDE
myRCTD.reps <- run.CSIDE.replicates(myRCTD.reps, exvar_list,population_de = T, 
                                    cell_type_threshold = 1,
                                    cell_types,
                                    doublet_mode = F,#gene_threshold=0.01,
                                    #fdr=0.05,
                                    weight_threshold = 0.9
                                    #,normalize_expr=T
                                   )
saveRDS(myRCTD.reps,"13cp.sham_H_astrocyte.gfap.vs.slc.4reps.cpt1.wt09.CSIDE.rds")

#get the DEG of regional astrocyte and enriched go terms
de_result<-myRCTD.reps@population_sig_gene_df
test<-de_result$Astrocyte
test<-arrange(test,log_fc_est)
write.csv(test,"13cp.sham_H_astr.gfap.vs.slc7a10_4reps.CSIDE.cpt1.wt09.sig_gene_df.csv")
backid<-read.table("WT.GOterm.SCT.backgeneids.txt")
backid<-backid$x
sf1<-"goterm_result/"
if(!dir.exists(sf1))
    dir.create(sf1)
pvalueCutoff=0.05
qvalueCutoff=0.05

#df<-data.frame()
for(i in names(de_result)){
    temp<-de_result[[i]]
    temp<-subset(temp,q_val<0.05&p<0.05)
    if(nrow(temp)>0){
        temp$gene<-rownames(temp)
        temp$group<-NA
        temp[temp$log_fc_est<0,"group"]<-paste0(i,"_sham.MG")
        temp[temp$log_fc_est>0,"group"]<-paste0(i,"_sham.WM")
        if(nrow(temp)>30){
            grlabs<-split(temp$gene,temp$group)
            gcSample<-lapply(grlabs,function(gr) as.numeric(bitr(gr,fromType="SYMBOL",toType = "ENTREZID",OrgDb="org.Mm.eg.db")$ENTREZID))
        
            xx.mus.go<-compareCluster(gcSample,OrgDb="org.Mm.eg.db",fun="enrichGO",
                                      pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff,
                                 ont="BP",readable=T,universe=backid)
            if(nrow(xx.mus.go@compareClusterResult)>0){
                saveRDS(xx.mus.go,paste0(sf1,i,"_H.4reps.13cp.sham.WM.vs.MG.cpt1.wt09.deg.pq005.go.rds"))
        }          
    }  
    } 
}

temp<-temp[temp$p.adjust<0.05&temp$Count>=5,]
df1<-subset(temp,Cluster=="Astrocyte_sham.WM")
df2<-subset(temp,Cluster=="Astrocyte_sham.MG")

df_ad<-rbind(df1,df2#,df3
            )
df_ad$Cluster<-factor(df_ad$Cluster,levels = c("Astrocyte_sham.WM","Astrocyte_sham.MG"))
df_ad<-df_ad[order(df_ad$Cluster),]
top10<-df_ad %>%group_by(.,Cluster) %>% top_n(.,-10,p.adjust) %>% top_n(.,10,Count)
top10<-arrange(top10,Cluster)
top10$Description<-factor(top10$Description,levels=rev(unique(top10$Description))) #order Description by cluster 
# barplot
options(repr.plot.width=8,repr.plot.height=8)
ggplot(top10,aes(x=Description,y=-log10(p.adjust),
              fill=Cluster))+
    geom_bar(stat="identity",position = "dodge")+
    coord_flip()+
    scale_fill_manual(values = rev(c("#DE8BE6","#911eb4")))+
    theme_bw()+
    theme(axis.text.y.left = element_text(size = 13),panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor=element_blank())
ggsave("13cp.Astrocyte_sham_H.WM.vs.MG.Go.withMetabolic.top10.barplot_220701.png",width = 8,height = 8,dpi = 300)

#vocano plot
test$celltype<-ifelse(test$log_fc_est>0,"Astrocyte.Gfap","Astrocyte.Slc7a10")
### reverse the direction 
test$log_fc_est<--test$log_fc_est
pdf("13cp.sham_H.Astr.DEG.vocalno.withlabel_reversed220701.plot.pdf",width = 15,height = 15)
ggplot(test,aes(x=log_fc_est,y=-log10(p),col=celltype,label=rownames(test)
               ))+
    geom_point(size=2)+
    geom_text_repel(max.overlaps = 18)+
    scale_color_manual(values = c("#911eb4","#DE8BE6"))+
    theme_minimal()+
    theme(panel.grid =  element_blank())+
    geom_vline(xintercept = c(-0.5,0.5),col="grey")+
    geom_hline(yintercept = -log10(0.05),col="grey")
dev.off()
#options(repr.plot.width=9,repr.plot.height=5)
#ggplot(test,aes(x=log_fc_est,y=-log10(q_val),col=celltype,label=rownames(test)))+
#    geom_point()+
#    geom_text_repel(max.overlaps = 18)+
#    scale_color_manual(values = c("#911eb4","#e6beff"))+
#    theme_minimal()+
#    theme(panel.grid =  element_blank())+
#    geom_vline(xintercept = c(-1,1),col="grey")+
#    geom_hline(yintercept = -log10(0.01),col="grey")
#ggsave("sham_H.Astrocyte_Gfap.vs.Slc7a10.DEG.vocalno.plot.png",width = 9,height = 5,dpi = 400)

###enriched regulons
gfap_deg<-rownames(test[test$celltype=="Astrocyte.Gfap",])
slc_deg<-rownames(test[test$celltype=="Astrocyte.Slc7a10",])

mods<-list(gfap=gfap_deg,slc=slc_deg)
regulon_list<-readRDS('515.regulon.genes.list.rds')
r_genes<-unique(unlist(regulon_list))
backgenes<-readRDS("WT.SCT.data.rds")
backgenes<-rownames(backgenes)
total_genes<-backgenes
t_l<-length(total_genes)
### create hypergeometric dataframe for each GOI(gene of interest)
mod_regulon_df_list<-list()
#inter_mod_reg_df<-data.frame(row.names = names(mods),module=names(mods))
#c<-2
f<-"astr_deg.regulon.intersect_genes/"
if(!dir.exists(f))
    dir.create(f)

for(i in names(regulon_list)){
    inter_list<-list()
    r_genes<-regulon_list[[i]]
    k<-length(r_genes)
    # find the module success gene number in target regulon,regulon number, remaining regulons size and module size
    hyper_df <- data.frame(matrix(nrow=length(mods),ncol = 9),row.names = names(mods))
    

    colnames(hyper_df) <- c("q",  #regulon in module
                            "m",  #module size
                            "n",  #regulon in background 
                            "k",  #regulon size
                            "p.val",
                            "p.val.ad",
                            "sig",
                            "ratio",
                            "score")

    for(j in names(mods)){
        m_genes<-mods[[j]]
        suc_genes<-intersect(m_genes,r_genes)
        suc_genes_l<-length(suc_genes)
        m<-length(m_genes)
        n<-t_l-m
        p<-round(phyper(suc_genes_l,m,n,k,lower.tail = F,log.p = F),4)
        r<-round(suc_genes_l/k,4)
        
        if(suc_genes_l==0){
            inter_list[[j]]<-""
        }else{
            inter_list[[j]]<-c(suc_genes)
        }
        
        
        hyper_df[j,c("q","m","n","k","p.val")]<-c(suc_genes_l,m,n,k,p)
        hyper_df$ratio<-round(hyper_df$q/hyper_df$k,4) #calculate the regulon enriched ratio in each module

    }
    
    saveRDS(inter_list,paste0(f,i,".intersect_genes.in.astr_deg.rds"))
    
    hyper_df$p.val.ad<-round(p.adjust(hyper_df$p.val,"bonferroni"),4)
    hyper_df$sig<-round((1-hyper_df$p.val.ad)/0.05,4)
    hyper_df[hyper_df$p.val.ad>=0.01,"sig"]<-0
    hyper_df$score<-hyper_df$sig*hyper_df$ratio
    mod_regulon_df_list[[i]]<-hyper_df
}
### create module-regulon p.val matrix
p_ma<-do.call(cbind.data.frame,mod_regulon_df_list)
p_ma<-as.matrix(p_ma[,grepl(".sig",colnames(p_ma))])
colnames(p_ma)<-gsub(".sig","",colnames(p_ma))
p_ma<-t(p_ma)
# create ratio matrix
r_ma<-do.call(cbind.data.frame,mod_regulon_df_list)
r_ma<-as.matrix(r_ma[,grepl("ratio",colnames(r_ma))])
colnames(r_ma)<-gsub(".ratio","",colnames(r_ma))
r_ma<-t(r_ma)
s_ma<-do.call(cbind.data.frame,mod_regulon_df_list)
s_ma<-as.matrix(s_ma[,grepl("score",colnames(s_ma))])
colnames(s_ma)<-gsub(".score","",colnames(s_ma))
s_ma<-t(s_ma)
write.csv(p_ma,"astrDEG_regulon.hyper.pad001.pval.ad_significance.matrix.csv")
write.csv(s_ma,"astrDEG-regulon.hyper.pad001.pval-ratio.score.matrix.csv")
saveRDS(mod_regulon_df_list,"astrDEG.regulon.hypergeomic.test.list.pad001.rds")
#dotplot
df<-do.call(rbind.data.frame,mod_regulon_df_list)
df$module<-str_split(rownames(df),"\\.",simplify = T)[,2]
df$regulon<-str_split(rownames(df),"\\.",simplify = T)[,1]
df_sub<-df[df$p.val.ad<0.01&df$q>=5&df$ratio>0.01,]
top5<-df_sub %>% group_by(.,module) %>% top_n(15,wt = -p.val.ad) #%>% top_n(15,wt = ratio)
top5<- top5 %>% group_by(.,module) %>% arrange(.,module)
top5$regulon<-factor(top5$regulon,levels = unique(top5$regulon))
options(repr.plot.width=4,repr.plot.height=5)
ggplot(top5,aes(x=regulon,y=module))+
    geom_point(aes(size=q,color=sig#,shape=sc_shared
                  ))+
    scale_color_gradientn(colours = c("darkblue","pink","red"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60,size = 12,vjust = 1,hjust = 1),
          axis.text.y.left = element_text(angle = 30,size = 12),
          axis.title = element_blank()
         )
ggsave("astrDEG.regulon.pad001.num5.pointplot.png",width = 4,height = 2.6,dpi = 300)




