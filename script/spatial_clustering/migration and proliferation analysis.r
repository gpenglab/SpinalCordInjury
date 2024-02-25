library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
source("script/self_function/save_pheatmap_pdf.R")

da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
###positive regulation of cell migration score
cm<-read_excel("lee/classifyBySpatialMarker/GO_term_positiveregulationofcellmigration_20220620_052409.xlsx")
cmgenes<-unique(cm$Symbol)
cmgenes<-intersect(cmgenes,rownames(da@assays$SCT@data))
da<-AddModuleScore(da,features = list(cmgenes),name = "positive_regulation_of_cell_migration")
summary(da$positive_regulation_of_cell_migration1)
da$pos_mig_bi<-as.character(ifelse(da$positive_regulation_of_cell_migration1>0.07204,1,0)) #binary based on mean
table(da$pos_mig_bi)
options(repr.plot.width=20,repr.plot.height=16)
ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes(fill= positive_regulation_of_cell_migration1))+
    #scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7,name = "Spectral")))(50))+
    scale_fill_gradientn(colours = c('#034f84','#92a8d1',#'#b8a9c9', 
    
        '#d6d4e0',#"#deeaee",
    #"#f0f0f0",
    #'#fff2df',
       
        #"#eeac99",
        '#f4a688',
        "#ED797B",                             
        #"#f7786b",                             
        
        "#d64161",
        "#c94c4c"#,
        #"#c83349"
                                    ))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
    )
ggsave("WT.merge.regress_CC.nC.mt.ident_positive_regulation_of_cell_migration.spatial_manually_220701.png",width=20,height=16,dpi=400)

### calculate the colocalization pval of cell type each timepoint -220627
ct_df<-read.csv("RCTD/testC5/WT_all/fcreg1_220503/16celltypes.dec_conf.nor.meta.csv")
rownames(ct_df)<-ct_df[,1]
ct_df<-ct_df[,-1]
ct<-colnames(ct_df)[46:61]
# remove Svep1 and Dendritic
ct<-ct[!ct%in%c("Astro.Svep1","Dendritic")]
wt_df<-ct_df[,ct]
wt_df_bi<-wt_df
wt_df_bi<-apply(wt_df_bi,2,function(x) (ifelse(x>0,1,0)))
me<-cbind(da@meta.data,wt_df_bi) 
## Fisher test and jaccard simmilarity
pval_df<-data.frame(row.names = unique(me$time))
spot_n<-data.frame(row.names = unique(me$time))
jc_df<-data.frame(row.names = unique(me$time))
#pcc_df<-data.frame(row.names = unique(me$time))
for(i in unique(me$time)){
    temp_me<-me[me$time==i,]
    cells<-rownames(temp_me)
    temp_wt<-wt_df_bi[cells,ct]
    pp1<-temp_me[,"pos_mig_bi"]
    lr_ind<-rownames(temp_me)[pp1==1]
    for(j in ct){
        if(sum(temp_wt[,j])==0){
            pval_df[i,j]<-1
            spot_n[i,j]<-0
            jc_df[i,j]<-0
        }
        else{
            # unbinary wt
            #cc2<-wt_df[cells,j]
            # binary wt
            pp2<-temp_me[,j]
            cp_ind<-rownames(temp_me)[pp2==1]
            # pcc<-cor(pp1,pp2)
            pp1<-factor(pp1,levels=c(0,1))
            pp2<-factor(pp2,levels=c(0,1))
            f_df<-table(pp1,pp2)
            p<-fisher.test(f_df)$p.val
            pval_df[i,j]<-p
            
            n<-f_df[2,2]
            #nn<-f_df[1,1]
            spot_n[i,j]<-n
            
            JC<-length(intersect(cp_ind,lr_ind))/length(unique(c(cp_ind,lr_ind)))
            jc_df[i,j]<-JC
            
            #cc<-cor(cc1,cc2)
            #pcc_df[i,j]<-cc
        }
    }
}
pval_l<-pval_df
pval_l$time<-rownames(pval_l)
pval_l<-melt(pval_l)
colnames(pval_l)<-c("time","celltype","pval")
pval_l$key<-paste0(pval_l$time,"-",pval_l$celltype)
nspot_l<-spot_n
nspot_l$time<-rownames(nspot_l)
nspot_l<-melt(nspot_l)
colnames(nspot_l)<-c("time","celltype","n_spot")
nspot_l$key<-paste0(nspot_l$time,"-",nspot_l$celltype)
jc_l<-jc_df
jc_l$time<-rownames(jc_l)
jc_l<-melt(jc_l)
colnames(jc_l)<-c("time","celltype","JC")
jc_l$key<-paste0(jc_l$time,"-",jc_l$celltype)
df_l<-merge(pval_l,nspot_l,by = "key")
dim(df_l)
df_l<-merge(df_l,jc_l,by = "key")
dim(df_l)
write.csv(df_l,"pos_mig_bi.celltype.colocalization.pval.nspot.csv")
df_l_ad<-df_l
pval_ma<-df_l_ad[,c("time","celltype","JC")]
pval_ma<-reshape(pval_ma,idvar = "time",timevar = "celltype",direction = "wide")
rownames(pval_ma)<-pval_ma[,1]
pval_ma<-pval_ma[,-1]
colnames(pval_ma)<-gsub("JC.","",colnames(pval_ma))
pval_ma<-as.matrix(pval_ma[c(4,2,1,3),])
## remove decreasing cell types: Neuron, Oligodendrocyte
pval_ma_s<-pval_ma[,!colnames(pval_ma)%in% c("Neuron","Oligodendrocyte")]
dim(pval_ma_s)
options(repr.plot.width=7,repr.plot.height=4)
p<-pheatmap(pval_ma_s,scale = "none",cluster_rows = F,cluster_cols = T
           )
save_pheatmap_pdf(p,"pos_mig_mean_bi.14ct.Fisher_test.pval.jc005.pheatmap_220727.pdf",width = 7,height = 4)

###proliferation score
genes<-c("Mki67","Pcna","Gab2","Mcm6","Ccnd1","E2f1","Top2a")
da<-AddModuleScore(da,features = list(genes),name = "proliferation_score")
options(repr.plot.width=20,repr.plot.height=16)

    ggplot(da@meta.data,aes(x.ad,y.ad))+geom_point(size=1.6,shape=21,stroke=0.15,aes_string(fill="proliferation_score1"))+
    scale_fill_gradientn(colours = c('#034f84','#92a8d1',#'#b8a9c9', 
    
        '#d6d4e0',#"#deeaee",
    #"#f0f0f0",
    #'#fff2df',
       
        #"#eeac99",
        '#f4a688',
        "#ED797B",                             
        #"#f7786b",                             
        
        "#d64161",
        "#c94c4c"#,
        #"#c83349"
                                    ))+
    xlab(paste0("")) +
    ylab(paste0("")) + 
    theme(panel.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),#element_text(size = 15,colour = "black"),
          axis.text.y = element_blank()#element_text(size = 20,colour = "black")
          # axis.ticks.x = element_blank(),
          )
ggsave("WT.merge.SCT_7g.proliferationScore.spatial_240203.png",width = 20,height = 16,dpi = 300)
saveRDS(da,"WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
summary(da$proliferation_score1)[[5]]
da$pos_poliferate_bi<-as.character(ifelse(da$proliferation_score1>0.058,1,0)) #binary based on Q3 value
table(da$pos_poliferate_bi)
### calculate the colocalization pval of cell type each timepoint      
## Fisher test and jaccard simmilarity
pval_df<-data.frame(row.names = unique(me$time))
spot_n<-data.frame(row.names = unique(me$time))
jc_df<-data.frame(row.names = unique(me$time))
for(i in unique(me$time)){
    temp_me<-me[me$time==i,]
    cells<-rownames(temp_me)
    temp_wt<-wt_df_bi[cells,ct]
    pp1<-temp_me[,"pos_poliferate_bi"]
    lr_ind<-rownames(temp_me)[pp1==1]
    for(j in ct){
        if(sum(temp_wt[,j])==0){
            pval_df[i,j]<-1
            spot_n[i,j]<-0
            jc_df[i,j]<-0
        }
        else{
            pp2<-temp_me[,j]
            cp_ind<-rownames(temp_me)[pp2==1]
            pp1<-factor(pp1,levels=c(0,1))
            pp2<-factor(pp2,levels=c(0,1))
            f_df<-table(pp1,pp2)
            p<-fisher.test(f_df)$p.val
            pval_df[i,j]<-p
            
            n<-f_df[2,2]
            spot_n[i,j]<-n
            
            JC<-length(intersect(cp_ind,lr_ind))/length(unique(c(cp_ind,lr_ind)))
            jc_df[i,j]<-JC
        }
    }
}

pval_l<-pval_df
pval_l$time<-rownames(pval_l)
pval_l<-melt(pval_l)
colnames(pval_l)<-c("time","celltype","pval")
pval_l$key<-paste0(pval_l$time,"-",pval_l$celltype)
head(pval_l)

nspot_l<-spot_n
nspot_l$time<-rownames(nspot_l)
nspot_l<-melt(nspot_l)
colnames(nspot_l)<-c("time","celltype","n_spot")
nspot_l$key<-paste0(nspot_l$time,"-",nspot_l$celltype)
head(nspot_l)

jc_l<-jc_df
jc_l$time<-rownames(jc_l)
jc_l<-melt(jc_l)
colnames(jc_l)<-c("time","celltype","JC")
jc_l$key<-paste0(jc_l$time,"-",jc_l$celltype)
head(jc_l)

df_l<-merge(pval_l,nspot_l,by = "key")
dim(df_l)
df_l<-merge(df_l,jc_l,by = "key")
dim(df_l)
write.csv(df_l,"pos_poliferate_bi_byq3.14celltype.colocalization_eachtime.pval.nspot_240125.csv")
df_l_ad<-df_l
pval_ma<-df_l_ad[,c("time","celltype","JC")]
head(pval_ma)
pval_ma<-reshape(pval_ma,idvar = "time",timevar = "celltype",direction = "wide")
rownames(pval_ma)<-pval_ma[,1]
pval_ma<-pval_ma[,-1]
colnames(pval_ma)<-gsub("JC.","",colnames(pval_ma))
pval_ma<-as.matrix(pval_ma[c(4,2,1,3),])
head(pval_ma)
## remove decreasing cell types: Neuron, Oligodendrocyte
pval_ma_s<-pval_ma[,!colnames(pval_ma)%in% c("Neuron","Oligodendrocyte")]
dim(pval_ma_s)
options(repr.plot.width=6.3,repr.plot.height=3.8)
p<-pheatmap(pval_ma_s,scale = "none",cluster_rows = F,cluster_cols = T#,color = brewer.pal(n = 7,name = "Reds")
           )
save_pheatmap_pdf(p,"240125_poliferation_q3_bi.14ct.Fisher_test.pval.jc005.pheatmap.pdf",width = 6.3,height = 3.8)

###barplot
## spot number of colocalization for each cell type
celltypes<-c("Astrocyte.Gfap","OPC","Fibroblast","R.Microglia","Vascular")
temp<-ct_df[,celltypes]
temp_bi<-apply(temp,2,FUN = function(x) ifelse(x>0,1,0))
temp_df<-cbind(da@meta.data[,c("pos_poliferate_bi","time")],temp_bi)
temp_df$pos_poliferate_bi<-as.numeric(temp_df$pos_poliferate_bi)
temp_df$time<-as.character(temp_df$time)
co_list<-list()
for (j in celltypes){
    co_df<-data.frame(row.names = unique(temp_df$time))
    co_df$yes<-0
    co_df$no<-0
    for (i in unique(temp_df$time)){
        temp2<-temp_df[temp_df$time==i,]
        s1<-rownames(temp2)[temp2$pos_poliferate_bi==1]
        s2<-rownames(temp2)[temp2[,j]==1]
        n1<-length(unique(intersect(s1,s2)))
        n2<-length(s2)-n1
        co_df[i,]<-c(n1,n2)
    }
    co_list[[j]]<-co_df
}
saveRDS(co_list,"240223.celltype.colocalized.with.poliferationScore.rds")
for (i in celltypes){
    test<-co_list[[i]]
    test$time<-rownames(test)
    test<-gather(test,group,num,yes:no)
    test$time<-factor(test$time,levels = unique(test$time))
    test$group<-factor(test$group,levels = c("yes","no"))
    options(repr.plot.width=4,repr.plot.height=6)
    ggplot(test,aes(x=time,y=num,fill=group))+
        geom_bar(stat="identity")+
        theme_classic()+
        scale_fill_manual(values = c("#c94c4c","#b2b2b2"))+
        ylim(0,4000)
    ggsave(paste0("240223_barplot_",i,"_colocalization_with_poliferateScore.png"),dpi = 400,width = 4,height = 6)
}





###colocalization of poliferation between Gfap and Igfbp2
me2<-da@meta.data
me2$Igfbp2<-da@assays$SCT@scale.data['Igfbp2',]
me2$Igfbp2_bi<-ifelse(me2$Igfbp2>0.1624,1,0) #binary based on Q3
me2$Gfap<-da@assays$SCT@scale.data['Gfap',]
me2$Gfap_bi<-ifelse(me2$Gfap>0.4738,1,0)
## Fisher test and jaccard simmilarity
pval_df<-data.frame(row.names = unique(me$time))
spot_n<-data.frame(row.names = unique(me$time))
jc_df<-data.frame(row.names = unique(me$time))
#pcc_df<-data.frame(row.names = unique(me$time))
gene<-c("Igfbp2_bi","Gfap_bi")
for(i in unique(me2$time)){
    temp_me<-me2[me2$time==i,]
    cells<-rownames(temp_me)
    temp_wt<-temp_me[cells,gene]
    pp1<-temp_me[,"pos_poliferate_bi"]
    lr_ind<-rownames(temp_me)[pp1==1]
    for(j in gene){
        if(sum(temp_wt[,j])==0){
            pval_df[i,j]<-1
            spot_n[i,j]<-0
            jc_df[i,j]<-0
        }
        else{
            pp2<-temp_me[,j]
            cp_ind<-rownames(temp_me)[pp2==1]
            # pcc<-cor(pp1,pp2)
            pp1<-factor(pp1,levels=c(0,1))
            pp2<-factor(pp2,levels=c(0,1))
            f_df<-table(pp1,pp2)
            p<-fisher.test(f_df)$p.val
            pval_df[i,j]<-p
            
            n<-f_df[2,2]
            #nn<-f_df[1,1]
            spot_n[i,j]<-n
            
            JC<-length(intersect(cp_ind,lr_ind))/length(unique(c(cp_ind,lr_ind)))
            jc_df[i,j]<-JC
        }
    }
}

pval_l<-pval_df
pval_l$time<-rownames(pval_l)
pval_l<-melt(pval_l)
colnames(pval_l)<-c("time","gene","pval")
pval_l$key<-paste0(pval_l$time,"-",pval_l$gene)
head(pval_l)

nspot_l<-spot_n
nspot_l$time<-rownames(nspot_l)
nspot_l<-melt(nspot_l)
colnames(nspot_l)<-c("time","gene","n_spot")
nspot_l$key<-paste0(nspot_l$time,"-",nspot_l$gene)
head(nspot_l)

jc_l<-jc_df
jc_l$time<-rownames(jc_l)
jc_l<-melt(jc_l)
colnames(jc_l)<-c("time","gene","JC")
jc_l$key<-paste0(jc_l$time,"-",jc_l$gene)
head(jc_l)
df_l<-merge(pval_l,nspot_l,by = "key")
dim(df_l)
df_l<-merge(df_l,jc_l,by = "key")
dim(df_l)
write.csv(df_l,"pos_poliferate_bi_byq3.Igfbp2_Gfap.colocalization_eachtime.pval.nspot_240125.csv")
df_l_ad<-df_l
pval_ma<-df_l_ad[,c("time","gene","JC")]
head(pval_ma)
pval_ma<-reshape(pval_ma,idvar = "time",timevar = "gene",direction = "wide")
rownames(pval_ma)<-pval_ma[,1]
pval_ma<-pval_ma[,-1]
colnames(pval_ma)<-gsub("JC.","",colnames(pval_ma))
pval_ma<-as.matrix(pval_ma[c(4,2,1,3),])
head(pval_ma)
pval_ma_s<-pval_ma
dim(pval_ma_s)
options(repr.plot.width=4.5,repr.plot.height=2)
p<-pheatmap(t(pval_ma),scale = "none",cluster_rows = F,cluster_cols = F
           )
save_pheatmap_pdf(p,"240125_pos_poliferate_bi_byq3.Igfbp2_Gfap_biByq3.colocalization.Fisher_test.pval.jc005.pheatmap.pdf",width = 4.5,height = 2)


