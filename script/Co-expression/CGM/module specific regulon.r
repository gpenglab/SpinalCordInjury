library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
#library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(scales)
source("script/self_function/save_pheatmap_pdf.R")

### load module gene
modgene<-read.csv("thr100.sp18.hclust_average.min30.deep2.h015.34modules.full.gene.list.csv")
rownames(modgene)<-modgene$gene
modgene<-modgene[!modgene$module=="grey",]
#load regulon auc matrix
auc<-read.csv("auc_mtx.csv")
#load spatial meta data
meta<-read.csv("WT_replace_v2/res02_220310/WT.SCT.pc20.k50.res02.meta.data.csv")
rownames(meta)<-meta[,1]
meta<-cbind(meta,auc)
write.csv(meta,"WT.SCT.reg.pc20.k50.res02.auc.meta.csv")
mods<-split(modgene$gene,modgene$module)
#load regulon
regulon<-read.csv("regulons.csv")
rownames(regulon)<-regulon$tfs
regulon$targets<-as.character(regulon$targets)
# get regulon genes list
regulon_list<-list()
for(i in regulon$tfs){
    regulon_list[[i]]<-c(i,unlist(strsplit(regulon[regulon$tfs==i,"targets"],",")))
}
saveRDS(regulon_list,"WT.regulon.list.rds")
r_genes<-unique(unlist(regulon_list))
##background genes
backgenes<-readRDS("/home/jovyan/result/cellchat/time_domain_220325/WT.SCT.data.rds")
backgenes<-rownames(backgenes)
total_genes<-backgenes
t_l<-length(total_genes)
### check if all module genes are contained in background genes
intgenes<-intersect(r_genes,total_genes)

### create hypergeometric dataframe for each module
mod_regulon_df_list<-list()
#inter_mod_reg_df<-data.frame(row.names = names(mods),module=names(mods))
#c<-2
f<-"intersect_genes/"
if(!dir.exists(f))
    dir.create(f)

for(i in names(regulon_list)){
    inter_list<-list()
    #inter_mod_reg_df[,i]<-NA
    r_genes<-regulon_list[[i]]
    #r_genes<-intersect(r_genes,intgenes) ### filter regulon genes that not exist in all module gene background
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
    
    saveRDS(inter_list,paste0(f,i,".intersect_genes.in.modules.rds"))
    
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
write.csv(p_ma,"module_regulon.hyper.pad001.pval.ad_significance.matrix.csv")
write.csv(s_ma,"module-regulon.hyper.pad001.pval-ratio.score.matrix.csv")
saveRDS(mod_regulon_df_list,"module.regulon.hypergeomic.test.list.pad001.rds")

df<-do.call(rbind.data.frame,mod_regulon_df_list)
df$module<-str_split(rownames(df),"\\.",simplify = T)[,2]
df$regulon<-str_split(rownames(df),"\\.",simplify = T)[,1]
#remove regulon size <=10
df<-subset(df,k>10)
top3<- df %>% subset(.,score>0) %>% group_by(module)  %>% top_n(.,wt = sig,n = 10) %>% top_n(.,wt = score,n = 3) #
top3<-arrange(top3,module)
df3<-p_ma[unique(top3$regulon),unique(top3$module)]
col_anno<-data.frame(row.names = colnames(df3),"module"=colnames(df3))
module_col<-colnames(df3)
names(module_col)<-colnames(df3)
col_list<-list(module=module_col)
options(repr.plot.width=7, repr.plot.height=12) 
png("selected.module.pad001.selfdefine.top3.score.regulons.heatmap.png",width=7,height=12,res = 300)
p<-pheatmap(df3, name="Score",cluster_rows = F,cluster_cols = F,border_color = NA,scale = "row",
         col=viridisLite::viridis(n = 50),annotation_col=col_anno,annotation_colors=col_list
        )
dev.off()   


