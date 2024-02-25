library(patchwork)
library(CellChat)
options(stringsAsFactors = FALSE)
library(Seurat)
library(igraph)
library(stringr)
library(stringi)
library(reshape2)
library(RColorBrewer)

#Get cellchat signaling gene list

### extract cellchat LR interaction to our data
CellChatDB<-CellChatDB.mouse
interaction_input<-CellChatDB$interaction

cc_lr<-as.data.frame(str_split(interaction_input$interaction_name_2," - ",simplify = T))
cc_r<-as.data.frame(str_split(gsub("\\)","",gsub("\\(","",cc_lr$V2)),"\\+",simplify = T))
colnames(cc_r)<-c("R1","R2")
cc_r$R1<-gsub("[[:space:]]","",cc_r$R1)
cc_r$R2<-gsub("[[:space:]]","",cc_r$R2,fixed = TRUE)
cc_lr<-cbind(cc_lr[,1],cc_r)
cc_lr<-apply(cc_lr,2,function(x) gsub("[[:space:]]","",x))

rownames(cc_lr)<-rownames(interaction_input)
colnames(cc_lr)<-c("ligand","receptor1","receptor2")

write.csv(cc_lr,"cellchat.all.lr.interaction.df.csv")

lr_genes<-unique(cc_lr[cc_lr!=""])

#Load spatial data
ma<-readRDS("WT.SCT.data.rds")
ma_bi<-ma

meta<-read.csv("WT.SCT.pc20.k50.res02.meta.data.csv")
rownames(meta)<-meta[,1]

#Load spatial cell type dataframe
wt_df<-read.csv("16celltypes.dec_conf.nor.meta.csv")
rownames(wt_df)<-wt_df[,1]
#colnames(wt_df)
celltypes<-colnames(wt_df)[c(47:ncol(wt_df))]

weight_df<-wt_df[,celltypes]
#remove Neuron and Oligodendrocyte noise through EM and then get new wt_df
summary(wt_df$Neuron[wt_df$Neuron!=0])
summary(wt_df$Oligodendrocyte[wt_df$Oligodendrocyte>0])
test<-wt_df
test$Neuron_bi<-ifelse(test$Neuron>0.34405,1,0)
test$Oligo_bi<-ifelse(test$Oligodendrocyte>0.44393,1,0)
options(repr.plot.width=20,repr.plot.height=18)
ggplot(test,aes(x=x.ad,y=y.ad,color=as.factor(Neuron_bi)))+
    geom_point(size=1.2)+
    scale_color_manual(values = c("darkblue","yellow"))+
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
ggsave("Neuron.mean.binary.spatial.png",width = 20,height = 18,dpi = 300)
options(repr.plot.width=20,repr.plot.height=18)
ggplot(test,aes(x=x.ad,y=y.ad,color=as.factor(Oligo_bi)))+
    geom_point(size=1.2)+
    scale_color_manual(values = c("darkblue","yellow"))+
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
ggsave("Oligodendrocyte.mean.binary.spatial.png",width = 20,height = 18,dpi = 300)
test$Neuron<-test$Neuron*test$Neuron_bi
test$Oligodendrocyte<-test$Oligodendrocyte*test$Oligo_bi
test<-test[,!colnames(test)%in%c("Neuron_bi","Oligo_bi")]
wt_df_ad<-test

weight_df_ad<-wt_df_ad[,celltypes]
write.csv(weight_df_ad,"16celltypes.dec_conf.nor_NOrmMean.csv")

sub_lr_c<-intersect(lr_genes,rownames(ma))

#Define LR interaction at spot by EM

### calculate lr average expression in each time 
me_exp<-readRDS("cellchat.lr.express.and.gm_express.meta.rds")
lr_time_mean<-data.frame(row.names = colnames(lr_bi_df))
for(i in unique(meta$time)){
    lr_time_mean[,i]<-0
    cells<-rownames(meta)[meta$time==i]
    for(j in rownames(lr_time_mean)){
        temp_m<-me_exp[cells,]
        s<-summary(temp_m[,j])
        em<-0.5*(s[[3]])+(s[[2]]+s[[5]])/4
        lr_time_mean[j,i]<-em
    }
}
lr_time_sum<-data.frame(row.names = colnames(lr_bi_df))
for(i in unique(meta$time)){
    lr_time_sum[,i]<-0
    cells<-rownames(meta)[meta$time==i]
    for(j in rownames(lr_time_sum)){
        temp_m<-me_exp[cells,]
        s<-sum(temp_m[,j])
        lr_time_sum[j,i]<-s
    }
}
write.csv(lr_time_mean,"lr.4times.gm.csv")

test_me<-meta
test_me<-cbind(test_me,lr_bi_df)

me_exp<-readRDS("cellchat.lr.express.and.gm_express.meta.rds")

###Firstly calculate the lR geometric mean (gm) across 16 samples
#subset represent samples
meta_sub<-subset(meta,sample %in% c('WT_sham_H_R2_1mm_2','WT_sham_H_R2_4','WT_sham_T_210323_4','WT_sham_T_R2_1mm_4',
                                                  'WT_3h_H_R2_1mm_2','WT_3h_H_R2_2','WT_3h_T_R2_1','WT_3h_T_210330_1mm_3',
                                                  'WT_24h_H_R1_1mm_3','WT_24h_H_R1_4','WT_24h_T_201231_3','WT_24h_T_R1_1mm_3',
                                                  'WT_72h_H_R1_1mm_2','WT_72h_H_210323_1','WT_72h_T_R2_4','WT_72h_T_R1_1mm_3'))
cells<-rownames(meta_sub)
ma_sub<-ma[sub_lr_c,cells]

lr_bi_df<-readRDS("SCT.data.LR.binary.RDS")

###Test celltype pair with LR at each time point
#1.firstly check the LR binary spatial distribution
test_me<-cbind(meta_sub,temp_lr_bi_df)
options(repr.plot.width=20,repr.plot.height=18)
ggplot(test_me,aes(x=x.ad,y=y.ad,color=as.factor(PSAP_GPR37L1)))+
    geom_point(size=1.2)+
    #scale_color_manual(values = rev(c("orange","#911eb4","#e6beff","blue","grey")))+
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())+
    scale_color_manual(values = c("#440154FF","#FDE725FF"))
ggsave("PSAP_GPR37L1.binary.spatial.png",width = 20,height =18,dpi = 300 )

#2.Fisher exact test on LR and CP
num_talk<-apply(lr_bi_df,2,sum)
# remove num_talk <20                           
lr_bi_df_sub<-lr_bi_df[,num_talk>20]
cal_jc<-function(weight_df){
    df<-data.frame(row.names = colnames(weight_df))
    for(i in colnames(weight_df)){
        df[,i]<-0
        for(j in colnames(weight_df)){
            jc<-sum(weight_df[,i]>0 & weight_df[,j]>0)/sum(weight_df[,i]>0 | weight_df[,j]>0)
            df[j,i]<-jc
        }
    }
    return(df)
}

for(s in unique(meta_sub$time)){
   ### construct expression matrix for each section
    cell<-meta_sub[meta_sub$time==s,"barcode"]
    temp_lr_bi_df<-lr_bi_df_sub[cell,]                           
    #cell-pair matrix
    temp_wt<-weight_df_ad[cell,]#weight_df[cell,]
    temp_wt<-temp_wt[,colSums(temp_wt)>1] #remove un-existed cell types
    wt_bi_df<-data.frame(row.names = rownames(temp_wt))
    cor_df<-cor(temp_wt,method = "spearman")
    #jc_df<-cal_jc(temp_wt)
    # construct cell type pairs
    ij <- combn(colnames(temp_wt), 2,simplify = T)
    length(ij)
    for(i in 1:ncol(ij)){
        pair<-ij[,i]
        cn<-paste(pair[1],pair[2],sep = "_")
        if(cor_df[pair[1],pair[2]]<(-0.4) & s %in% c("WT_sham","WT_3h")){
            wt_bi_df[,cn]<-0
            }
        #if(jc_df[pair[1],pair[2]]<0.06){
        #    wt_bi_df[,cn]<-0
        #    }
        else{
            inter_spots1<-rownames(temp_wt[temp_wt[,pair[1]]>0 & temp_wt[,pair[2]]>0,])
            wt_bi_df[,cn]<-0
            if(length(inter_spots1)>10)
                wt_bi_df[inter_spots1,cn]<-1  
            }
              
        }
    head(wt_bi_df)
    f2<-"cell-pair.colocalize.binary_result2/"
    if(!dir.exists(f2))
        dir.create(f2)
    write.csv(wt_bi_df,paste0(f2,s,".binaried.celltype_pair.csv"))
    cp_bi_df_sub<-wt_bi_df[,colSums(wt_bi_df)!=0]
    
    ###Fisher's exact test
    #calculate LR's association with cellpair and get the LR-cp double sig spot number
    lr_cp_pval_df<-data.frame(row.names = colnames(temp_lr_bi_df))
    lr_cp_spot_n<-data.frame(row.names = colnames(temp_lr_bi_df))
    lr_cp_jc_df<-data.frame(row.names = colnames(temp_lr_bi_df))
    #lr_cp_pcc_df<-data.frame(row.names = colnames(temp_lr_bi_df))                          
    for(i in colnames(cp_bi_df_sub)){
        lr_cp_pval_df[,i]<-NA
        lr_cp_spot_n[,i]<-NA
        for(j in colnames(temp_lr_bi_df)){
            #if(all(table(temp_temp_lr_bi_df[,j],cp_bi_df_sub[,i])==2)){
            pp1<-temp_lr_bi_df[,j]
            lr_ind<-rownames(temp_lr_bi_df)[pp1==1]
            pp2<-cp_bi_df_sub[,i]
            cp_ind<-rownames(cp_bi_df_sub)[pp2==1]
            # pcc<-cor(pp1,pp2)
            pp1<-factor(pp1,levels=c(0,1))
            pp2<-factor(pp2,levels=c(0,1))
            f_df<-table(pp1,pp2)
            p<-fisher.test(f_df,alternative = "greater")$p.val
            lr_cp_pval_df[j,i]<-p
            
            n<-f_df[2,2]
            #nn<-f_df[1,1]
            lr_cp_spot_n[j,i]<-n
            
        }
    }
    
    lr_cp_n_long<-lr_cp_spot_n
    lr_cp_n_long$lr<-rownames(lr_cp_n_long)
    lr_cp_n_long<-melt(lr_cp_n_long)
    colnames(lr_cp_n_long)<-c("lr","cp","n_spot")
    lr_cp_n_long$key<-paste0(lr_cp_n_long$lr,"-",lr_cp_n_long$cp)
    head(lr_cp_n_long)
    ## add the most possible cp of each LR and its pval
    f3<-"LR.cp.Fisher.pval_result2/"
    if(!dir.exists(f3))
        dir.create(f3)
    f3a<-"LR.cp.Fisher.pval_result2/pval.df/"
    if(!dir.exists(f3a))
        dir.create(f3a)                       
    write.csv(lr_cp_pval_df,paste0(f3a,s,".cellpair.LR.fishers.pval.df"))
    
    lr_cp_long_df<-lr_cp_pval_df#[,-which(colnames(lr_cp_pval_df)%in%c("top1.cp","top1.cp.pval"))]
    lr_cp_long_df$lr<-rownames(lr_cp_long_df)
    lr_cp_long_df<-melt(lr_cp_long_df)
    lr_cp_long_df$key<-paste0(lr_cp_long_df$lr,"-",lr_cp_long_df$variable)
    colnames(lr_cp_long_df)<-c("lr","cp","pval","key")
    head(lr_cp_long_df)                      
    
    lr_cp_long<-merge(lr_cp_long_df,lr_cp_n_long[,c(3,4)],by = "key")                  
    head(lr_cp_long)
    dim(lr_cp_long)
    
    lr_cp_long_sub<-subset(lr_cp_long,n_spot>10&pval<0.01)
    #calculate cp's lr average expression
    lr_cp_long_sub$EM<-0
    lr_cp_long_sub$cp_lr_mean<-0
    for(i in unique(lr_cp_long_sub$cp)){
        cp_spots<-rownames(cp_bi_df_sub)[cp_bi_df_sub[,i]==1]
        cp_lr<-unique(lr_cp_long_sub[lr_cp_long_sub$cp==i,"lr"])
        for(j in cp_lr){
            cp_lr_MExp<-mean(me_exp[cp_spots,j])
            lr_cp_long_sub[lr_cp_long_sub$cp==i & lr_cp_long_sub$lr==j,"EM"]<-lr_bi_gm[[j]]
            lr_cp_long_sub[lr_cp_long_sub$cp==i & lr_cp_long_sub$lr==j,"cp_lr_mean"]<-cp_lr_MExp
        }
    }
    f3b<-"LR.cp.Fisher.pval_result2/pval.long.df/"
    if(!dir.exists(f3b))
        dir.create(f3b) 
    write.csv(lr_cp_long,paste0(f3b,s,".cellpair.LR.pval.spots.long.all.csv"))
    write.csv(lr_cp_long_sub,paste0(f3b,s,".cellpair.LR.pval.spots.long.csv"))            
}  


#Just calculate sample1 and directly get the result of each timepoint
f3b<-"LR.cp.Fisher.pval_result2/pval.long.df/"
f<-dir(f3b)
f<-f[c(2,4,6,8)]
temp1<-read.csv(paste0(f3b,f[4]))
temp1<-temp1[,-1]
temp1<-subset(temp1,pval<1e-3&n_spot>10)
# load sham data
temp1<-read.csv(paste0(f3b,f[4]))
temp1<-temp1[,-1]
temp1<-subset(temp1,pval<=1e-3&n_spot>10)
#colnames(temp1)<-c("key","lr","cp","pval","n_spot","EM")
temp1$autocrine<-c(str_split(temp1$lr,"_",simplify = T)[,1]==str_split(temp1$lr,"_",simplify = T)[,2])
head(temp1)
temp2<-read.csv(paste0(f3b,f[2]))
temp2<-temp2[,-1]
temp2<-subset(temp2,pval<=1e-3&n_spot>10)
#colnames(temp2)<-c("key","lr","cp","pval","n_spot","JC")
temp2$autocrine<-c(str_split(temp2$lr,"_",simplify = T)[,1]==str_split(temp2$lr,"_",simplify = T)[,2])
head(temp2)
temp3<-read.csv(paste0(f3b,f[1]))
temp3<-temp3[,-1]
temp3<-subset(temp3,pval<=1e-3&n_spot>10)
#colnames(temp3)<-c("key","lr","cp","pval","n_spot","JC")
temp3$autocrine<-c(str_split(temp3$lr,"_",simplify = T)[,1]==str_split(temp3$lr,"_",simplify = T)[,2])
head(temp3)
temp4<-read.csv(paste0(f3b,f[3]))
temp4<-temp4[,-1]
temp4<-subset(temp4,pval<=1e-3&n_spot>10)
#colnames(temp4)<-c("key","lr","cp","pval","n_spot","JC")
temp4$autocrine<-c(str_split(temp4$lr,"_",simplify = T)[,1]==str_split(temp4$lr,"_",simplify = T)[,2])
head(temp4)

sf3<-"4time.spmc-04.pval1e-4_spot10.lr.cp_sc2/"
    if(!dir.exists(sf3))
        dir.create(sf3)

# add crosstalk shared with single cell result at corresponding timepoint
sc_df<-read.csv("uninj.fc05.pc05.min10.interaction.net.csv")
sc_df$source[sc_df$source %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$target[sc_df$target %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$key1<-paste0(sc_df$interaction_name,"-",sc_df$source,"_",sc_df$target)
sc_df$key2<-paste0(sc_df$interaction_name,"-",sc_df$target,"_",sc_df$source)
head(sc_df)
sc_int<-intersect(unique(temp1$key),unique(c(sc_df$key1,sc_df$key2)))
sc_int                   
temp1$sc_shared<-"N"
temp1$sc_shared<-ifelse(temp1$key %in% sc_int,"Y","N")
head(temp1)
table(temp1$sc_shared)

write.csv(temp1,paste0(sf3,"sham.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))     

temp2$sc_shared<-"N"
#temp2$sc_shared<-ifelse(temp2$key %in% sc_int,"Y","N")
head(temp2)
table(temp2$sc_shared)

write.csv(temp2,paste0(sf3,"3h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))     

# add crosstalk shared with single cell result at corresponding timepoint
sc_df<-read.csv("d1.fc05.pc05.min10.interaction.net.csv")
sc_df$source[sc_df$source %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$target[sc_df$target %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$key1<-paste0(sc_df$interaction_name,"-",sc_df$source,"_",sc_df$target)
sc_df$key2<-paste0(sc_df$interaction_name,"-",sc_df$target,"_",sc_df$source)
head(sc_df)
sc_int<-intersect(unique(temp3$key),unique(c(sc_df$key1,sc_df$key2)))
sc_int                   
temp3$sc_shared<-"N"
temp3$sc_shared<-ifelse(temp3$key %in% sc_int,"Y","N")
head(temp3)
table(temp3$sc_shared)

write.csv(temp3,paste0(sf3,"24h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))     

# add crosstalk shared with single cell result at corresponding timepoint
sc_df<-read.csv("d3.fc05.pc05.min10.interaction.net.csv")
sc_df$source[sc_df$source %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$target[sc_df$target %in% c("WM.Gfap","GM.Gfap")]<-"Astrocyte.Gfap"
sc_df$key1<-paste0(sc_df$interaction_name,"-",sc_df$source,"_",sc_df$target)
sc_df$key2<-paste0(sc_df$interaction_name,"-",sc_df$target,"_",sc_df$source)
head(sc_df)
sc_int<-intersect(unique(temp4$key),unique(c(sc_df$key1,sc_df$key2)))
sc_int                   
temp4$sc_shared<-"N"
temp4$sc_shared<-ifelse(temp4$key %in% sc_int,"Y","N")
head(temp4)
table(temp4$sc_shared)

write.csv(temp4,paste0(sf3,"72h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))     

#4 timepoints merged 
temp1<-read.csv(paste0(sf3,"sham.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))
temp1<-temp1[,-1]
temp2<-read.csv(paste0(sf3,"3h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))
temp2<-temp2[,-1]
temp3<-read.csv(paste0(sf3,"24h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))
temp3<-temp3[,-1]
temp4<-read.csv(paste0(sf3,"72h.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))
temp4<-temp4[,-1]

temp1$time<-"sham"
temp2$time<-"3h"
temp3$time<-"24h"
temp4$time<-"72h"

temp<-Reduce(rbind,list(temp1,temp2,temp3,temp4))

### centralize each lr
temp_scale<-temp
temp_scale$cp_lr_mean_scaled<-0

for(i in unique(temp$lr)){
    lr_temp<- temp_scale[temp_scale$lr==i,]
    lr_min<-min(lr_temp$cp_lr_mean)
    lr_max<-max(lr_temp$cp_lr_mean)
    lr_em<-unique(temp_scale[temp_scale$lr==i,"EM"])
    temp_scale[temp_scale$lr==i,"cp_lr_mean_scaled"]<-(temp_scale[temp_scale$lr==i,"cp_lr_mean"]-lr_em)/(lr_max-lr_em)
}

temp<-temp_scale
write.csv(temp,paste0(sf3,"4times.merged.spmc-04_2.pval1e-3_spot10.lr.cp.long.df_sc.csv"))

###Visualize
sf2<-"geompoint.cp_spmc-04_2.all_top3.lr.png2/"
if(!dir.exists(sf2))
    dir.create(sf2)
#remove Dendritic and Astro.Svep1
temp<-temp[!grepl("Astro.Svep1",temp$cp),]
### adjust
unique(temp$cp)
temp$cp<-factor(temp$cp,levels = c("Astrocyte.Slc7a10_Neuron","Astrocyte.Gfap_Neuron","Neuron_Vascular",
                                  "Homeostatic.Microglia_Neuron","Neuron_R.Microglia","Neuron_OPC","Neuron_Oligodendrocyte",
                                   "Astrocyte.Gfap_Oligodendrocyte","Oligodendrocyte_Vascular",'Oligodendrocyte_Ependymal',"Fibroblast_Oligodendrocyte",
                                   "Homeostatic.Microglia_Oligodendrocyte","Oligodendrocyte_R.Microglia","Oligodendrocyte_Myeloid",
                                   "Monocyte_Oligodendrocyte",'Neutrophil_Oligodendrocyte',"Oligodendrocyte_OPC","Oligodendrocyte_Macrophage","Astrocyte.Gfap_Vascular",
                                   "Astrocyte.Slc7a10_Vascular","Astrocyte.Gfap_Fibroblast","Astrocyte.Gfap_Homeostatic.Microglia",
                                   "Astrocyte.Gfap_R.Microglia","Astrocyte.Gfap_Myeloid","Astrocyte.Slc7a10_OPC","Astrocyte.Gfap_OPC",
                                   "Astrocyte.Gfap_Neutrophil","Astrocyte.Gfap_Monocyte","Astrocyte.Gfap_Ependymal","Astrocyte.Gfap_Macrophage","Fibroblast_R.Microglia",
                                   "Vascular_Macrophage","Fibroblast_Macrophage","OPC_Macrophage","Fibroblast_OPC",'Fibroblast_Monocyte','Vascular_OPC','OPC_R.Microglia'
                                  )
               )


temp$time<-factor(temp$time,levels = c("sham","3h","24h","72h"))
temp<-arrange(temp,time)
temp$group<-paste0(temp$time,"-",temp$cp)
temp$group<-factor(temp$group,levels = unique(temp$group))
top5<-temp %>% filter(autocrine==FALSE)%>%
    group_by(group) %>% top_n(2,wt = -pval) %>%  top_n(2,wt=n_spot) %>% top_n(2,wt=EM) 

top5<-temp[temp$lr%in%unique(top5$lr),]

#clustering LR
test<-lr_time_sum[unique(top5$lr),]
test2<-t(test)
# check the variance of each lr to define its temporal pattern
lr_sd<-apply(as.data.frame(test2),2,sd)             
test2<-as.data.frame(apply(test2,2,function(x) order(x)))
sim<-cor(test2,method = "spearman")
h<-hclust(dist(x=sim),method = "complete")
lr_order<-rev(rownames(sim)[h$order])

### directly look at the spatial pattern
fd<-"cp.top2.lr.spatial/"
if(!dir.exists(fd))
    dir.create(fd)
for(i in lr_order){
    options(repr.plot.width=20,repr.plot.height=18)
    ggplot(me_exp,aes_string(x='x.ad',y='y.ad',color=i))+
        geom_point(size=1.2)+
        #scale_color_manual(values = rev(c("orange","#911eb4","#e6beff","blue","grey")))+
        theme(panel.background = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())+
        scale_color_gradientn(colours =c("#DCDCDC","#A52A2A"))
        #scale_color_manual(values = c("#d6d4e0","#c83349"))
        #scale_color_manual(values = c('#f2e394','#d5f4e6','#80ced6','#eea29a','#f9d5e5'#'#5b9aa0'#
         #                             ,'#618685'#,
                                      
         #                             ,'#d6cbd3','#a2836e'
         #                             ))
    ggsave(paste0(fd,i,".SCT.data.GM.spatial.png"),width = 20,height = 18,dpi = 400)
}

### manually check throught time trend
lr_order<-c('EFNB3_EPHB1','PSAP_GPR37','SEMA4D_PLXNB3','SEMA4D_PLXNB1','PTN_PTPRZ1','PTN_SDC3',
            'FGF1_FGFR1','FGF1_FGFR3','FGF1_FGFR2','FGF9_FGFR3','GAS6_TYRO3','NFASC_CNTN1_CNTNAP1','PENK_OPRL1',#decrease
            'NCAM1_L1CAM','NRXN1_NLGN3','NRXN3_NLGN1','NRXN3_NLGN2',
            'NRXN3_NLGN3','NRXN2_NLGN3','FGF9_FGFR1', #switch
            'PSAP_GPR37L1','CXCL12_ACKR3','ADCYAP1_ADCYAP1R1','IL1B_IL1R2',
            'GRN_SORT1','CSF1_CSF1R','PDGFA_PDGFRA','C3_C3AR1',
            'SEMA6D_PLXNA1_TREM2','BDNF_NTRK2','SEMA6D_PLXNA1',
            'COL1A1_SDC4','COL6A2_ITGAV_ITGB8',
            'COL1A2_CD44','COL1A1_CD44','COL1A2_SDC4','COL6A2_SDC4','COL4A1_SDC4',
            'COL9A3_CD44','COL6A3_CD44','COL6A3_SDC1','COL6A3_SDC4','COL4A2_SDC4',
            'LAMB2_SV2C','LAMB2_CD44',
            'FN1_SDC1','FN1_CD44','FN1_ITGAV_ITGB1','FN1_ITGA5_ITGB1',
            'ANGPTL4_CDH5','ANGPTL4_SDC4',
            'SPP1_CD44','SPP1_ITGA5_ITGB1','SPP1_ITGAV_ITGB5',
            'THBS1_CD36','THBS1_SDC4','THBS1_CD47','THBS1_SDC1','THBS2_SDC4',
            'TNC_SDC4','TNC_SDC1',
            'TNR_SDC1','TNR_SDC4' #increase
            )

top5$lr<-factor(top5$lr,levels = lr_order)

options(repr.plot.width=25,repr.plot.height=13)
p<-ggplot(top5,aes(x=cp,y=lr))+
    geom_point(aes(size=-log10(pval),color=cp_lr_mean_scaled,shape=sc_shared
                  ))+
    scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7,name = "Spectral")))(20)#c("#8227d8","#f7c5c4","#e85957"
                                    #  ,"#e00512")
                         )+ 
    theme_bw()+#ggtitle(label = "            sham                                   3h                                   24h                                             72h")+
    theme(axis.text.x = element_text(angle = 60,size = 12,vjust = 1,hjust = 1),
          axis.text.y.left = element_text(angle = 30,size = 12),axis.title = element_blank(),
         strip.text=element_text(color="white")          
         )+
    facet_grid(~time)

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

pdf(paste0(sf2,"alltime.pval1e-3.n10.cp.top2pval.nspot_220825_sc.pdf"),width = 25,height = 13)
grid::grid.draw(g)
#ggsave(paste0(sf2,"alltime.pval005.n10.jc005.cp.top2pval.nspot_220725.png"),width = 25,height = 13,dpi = 300)
dev.off()

#### just show switch LR -221011
lr<-c('NCAM1_L1CAM','NRXN1_NLGN3','NRXN3_NLGN1','NRXN3_NLGN2','NRXN3_NLGN3','NRXN2_NLGN3','FGF9_FGFR1')
#cp<-c("Astrocyte.Slc7a10_Neuron","Astrocyte.Gfap_Neuron","Neuron_Vascular","Homeostatic.Microglia_Neuron","Neuron_R.Microglia","Neuron_OPC","Neuron_Oligodendrocyte")
top5_sub<-top5[ top5$lr %in% lr,]
dim(top5_sub)

options(repr.plot.width=7,repr.plot.height=4)
p<-ggplot(top5_sub,aes(x=cp,y=lr))+
    geom_point(aes(size=-log10(pval),color=cp_lr_mean_scaled#,shape=sc_shared
                  ))+
    scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7,name = "Spectral")))(20)#c("#8227d8","#f7c5c4","#e85957"
                                    #  ,"#e00512")
                         )+ 
    theme_bw()+#ggtitle(label = "            sham                                   3h                                   24h                                             72h")+
    theme(axis.text.x = element_text(angle = 60,size = 12,vjust = 1,hjust = 1),
          axis.text.y.left = element_text(angle = 30,size = 12),axis.title = element_blank(),
         strip.text=element_text(color="white")          
         )+
    facet_grid(~time)

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

pdf(paste0(sf2,"alltime.pval1e-3.n10.cp.top2pval.nspot_221011_switchLR.pdf"),width = 7,height = 4)
grid::grid.draw(g)
#ggsave(paste0(sf2,"alltime.pval005.n10.jc005.cp.top2pval.nspot_220725.png"),width = 25,height = 13,dpi = 300)
dev.off()

###Interaction network
sf<-"igraph.cellinteraction.number2/"
if(!dir.exists(sf))
    dir.create(sf)

geometric.mean<-function(x) {x=exp(mean(log(x)))
                            return(x)
                            }
temp1<-temp[temp$time=="sham",]
head(temp1)
temp2<-temp[temp$time=="3h",]
head(temp2)
temp3<-temp[temp$time=="24h",]
head(temp3)
temp4<-temp[temp$time=="72h",]
head(temp4)

sn1<-as.data.frame(table(temp1$cp))
sn1$source<-str_split(sn1$Var1,"_",simplify = T)[,1]
sn1$target<-str_split(sn1$Var1,"_",simplify = T)[,2]
rownames(sn1)<-sn1[,1]
sn1<-sn1[,-1]
colnames(sn1)<-c("value","from","to")
sn1<-sn1[,c(2,3,1)]
sn1$pval<- -log10(aggregate(temp1$pval,list(temp1$cp),FUN=geometric.mean)[,2])
sn1$n_spot<- round(aggregate(temp1$n_spot,list(temp1$cp),FUN=mean)[,2],2)
sn1$exp<-round(aggregate(temp1$cp_lr_mean_scaled,list(temp1$cp),FUN=mean)[,2],2)
head(sn1)
write.csv(sn1,paste0(sf,"sham.spmc-04.pval1e-3_spot10.igraph.df.csv"))

sn2<-as.data.frame(table(temp2$cp))
sn2$source<-str_split(sn2$Var1,"_",simplify = T)[,1]
sn2$target<-str_split(sn2$Var1,"_",simplify = T)[,2]
rownames(sn2)<-sn2[,1]
sn2<-sn2[,-1]
colnames(sn2)<-c("value","from","to")
sn2<-sn2[,c(2,3,1)]
sn2$pval<- -log10(aggregate(temp2$pval,list(temp2$cp),FUN=geometric.mean)[,2])
sn2$n_spot<- round(aggregate(temp2$n_spot,list(temp2$cp),FUN=mean)[,2],2)
sn2$exp<-round(aggregate(temp2$cp_lr_mean_scaled,list(temp2$cp),FUN=mean)[,2],2)
write.csv(sn2,paste0(sf,"3h.spmc-04.pval1e-3_spot10.igraph.df.csv"))

sn3<-as.data.frame(table(temp3$cp))
sn3$source<-str_split(sn3$Var1,"_",simplify = T)[,1]
sn3$target<-str_split(sn3$Var1,"_",simplify = T)[,2]
rownames(sn3)<-sn3[,1]
sn3<-sn3[,-1]
colnames(sn3)<-c("value","from","to")
sn3<-sn3[,c(2,3,1)]
sn3$pval<- -log10(aggregate(temp3$pval,list(temp3$cp),FUN=geometric.mean)[,2])
sn3$n_spot<- round(aggregate(temp3$n_spot,list(temp3$cp),FUN=mean)[,2],2)
sn3$exp<-round(aggregate(temp3$cp_lr_mean_scaled,list(temp3$cp),FUN=mean)[,2],2)
write.csv(sn3,paste0(sf,"24h.pval1e-4_spot10.igraph.df.csv"))

sn4<-as.data.frame(table(temp4$cp))
sn4$source<-str_split(sn4$Var1,"_",simplify = T)[,1]
sn4$target<-str_split(sn4$Var1,"_",simplify = T)[,2]
rownames(sn4)<-sn4[,1]
sn4<-sn4[,-1]
colnames(sn4)<-c("value","from","to")
sn4<-sn4[,c(2,3,1)]
sn4$pval<- -log10(aggregate(temp4$pval,list(temp4$cp),FUN=geometric.mean)[,2])
sn4$n_spot<- round(aggregate(temp4$n_spot,list(temp4$cp),FUN=mean)[,2],2)
sn4$exp<-round(aggregate(temp4$cp_lr_mean_scaled,list(temp4$cp),FUN=mean)[,2],2)
write.csv(sn4,paste0(sf,"72h.pval1e-3_spot10.igraph.df.csv"))

nodes<-data.frame(name=c('R.Microglia','Macrophage','Homeostatic.Microglia','Astrocyte.Gfap',
                         'Vascular',#'Astro.Svep1',
                         'Ependymal','Oligodendrocyte','Fibroblast',
                         'Neuron','Astrocyte.Slc7a10','OPC','Monocyte',
                         "Neutrophil","Myeloid"
                         ))
test<-sn1
#colnames(test)<-c('Var1','n','from','to')
test$n_spot<- scale(as.numeric(test$n_spot), center = 0)
        links <- data.frame(
            source = test$from,
            target = test$to,
            importance = as.numeric(test$n_spot)
          )

        graph_ntw <- igraph::graph_from_data_frame(links,
                vertices = nodes,
                directed = FALSE)
        
        deg <- degree(graph_ntw, mode="all")
        # Get color palette for difusion
        edge_importance <- E(graph_ntw)$importance #
        edge_col<-test$exp
        # Select a continuous palette
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
        
        # Create a color palette
        getPalette <- colorRampPalette(#viridisLite::viridis(n = 10,option = "G"))
                                       rev(c("#0C3321","#155637","#1a6a44","#1e7b4e","#2A6A63","#50817A","#6AA39A","#98BFB9","#C5DBD8")))

        
        #getPalette <- colorRampPalette(rev(c("#1D487C","#2864AD","#3d80d2","#6399db","#99bde7")))
        # Get how many values we need
        grad_edge <- seq(0, 1, 0.01)
        # Generate extended gradient palette dataframe
        graph_col_df <- data.frame(value = as.character(grad_edge),
                                   color = getPalette(length(grad_edge)),
                                   stringsAsFactors = FALSE)
        # Assign color to each edge
        color_edge <- data.frame(value = as.character(round(edge_col, 2)), stringsAsFactors = FALSE) %>%
          dplyr::left_join(graph_col_df, by = "value") %>%
          dplyr::pull(color)
        pdf(file=paste0(sf,"sham.spmc-04.pval1e-3_spot10.cp.lr.nspot_width.number_color.igraph.circle.plot_220725.pdf"))
        plot(graph_ntw,
             # Size of the edge
             edge.width = edge_importance*10, #*2 enhance the visualization
             edge.color = color_edge,
             # Size of the buble
             vertex.size = deg*4,
             vertex.color = "#E4791B",
             vertex.frame.color = "grey",
             vertex.label.color = "black",
             #vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
             layout = layout.circle)
       dev.off()

test<-sn2
#colnames(test)<-c('Var1','n','from','to')
test$n_spot<- scale(as.numeric(test$n_spot), center = 0)
        links <- data.frame(
            source = test$from,
            target = test$to,
            importance = as.numeric(test$n_spot)
          )
        #head(links)
        #nodes <- data.frame(name=unique(c(test$from,test$to)))
        #nodes
        graph_ntw <- igraph::graph_from_data_frame(links,
                vertices = nodes,
                directed = FALSE)
        
        deg <- degree(graph_ntw, mode="all")
        #pdf(paste0(result_dir,"spatial_celltype_withinspot_interaction_probability.circle.plot.png"))
        # Get color palette for difusion
        edge_importance <- E(graph_ntw)$importance #
        edge_col<-test$exp
        # Select a continuous palette
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
        
        # Create a color palette
        getPalette <- colorRampPalette(rev(c("#0C3321","#155637","#1a6a44","#1e7b4e","#2A6A63","#50817A","#6AA39A","#98BFB9","#C5DBD8")))
        #getPalette <- colorRampPalette(rev(c("#1D487C","#2864AD","#3d80d2","#6399db","#99bde7")))
        # Get how many values we need
        #grad_edge <- seq(3, 16
        #                 , 0.1)
        grad_edge <- seq(0, 1, 0.01)
        # Generate extended gradient palette dataframe
        graph_col_df <- data.frame(value = as.character(grad_edge),
                                   color = getPalette(length(grad_edge)),
                                   stringsAsFactors = FALSE)
        # Assign color to each edge
        color_edge <- data.frame(value = as.character(round(edge_col, 2)), stringsAsFactors = FALSE) %>%
          dplyr::left_join(graph_col_df, by = "value") %>%
          dplyr::pull(color)
        pdf(file=paste0(sf,"3h.spmc-04.pval1e-3_spot10.cp.lr.nspot_width.number_color.igraph.circle.plot.pdf"))
        plot(graph_ntw,
             # Size of the edge
             edge.width = edge_importance*10, #*2 enhance the visualization
             edge.color = color_edge,
             # Size of the buble
             vertex.size = deg*4,
             vertex.color = "#E4791B",
             vertex.frame.color = "grey",
             vertex.label.color = "black",
             #vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
             layout = layout.circle)
       dev.off()

test<-sn3
#colnames(test)<-c('Var1','n','from','to')
test$n_spot<- scale(as.numeric(test$n_spot), center = 0)
        links <- data.frame(
            source = test$from,
            target = test$to,
            importance = as.numeric(test$n_spot)
          )
        #head(links)
        #nodes <- data.frame(name=unique(c(test$from,test$to)))
        #nodes
        graph_ntw <- igraph::graph_from_data_frame(links,
                vertices = nodes,
                directed = FALSE)
        
        deg <- degree(graph_ntw, mode="all")
        #pdf(paste0(result_dir,"spatial_celltype_withinspot_interaction_probability.circle.plot.png"))
        # Get color palette for difusion
        edge_importance <- E(graph_ntw)$importance #
        edge_col<-test$exp
        # Select a continuous palette
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
        
        # Create a color palette
        getPalette <- colorRampPalette(rev(c("#0C3321","#155637","#1a6a44","#1e7b4e","#2A6A63","#50817A","#6AA39A","#98BFB9","#C5DBD8")))
        #getPalette <- colorRampPalette(rev(c("#1D487C","#2864AD","#3d80d2","#6399db","#99bde7")))
        # Get how many values we need
       # grad_edge <- seq(3, 16
       #                  , 0.1)
        
        grad_edge <- seq(0, 1, 0.01)

        # Generate extended gradient palette dataframe
        graph_col_df <- data.frame(value = as.character(grad_edge),
                                   color = getPalette(length(grad_edge)),
                                   stringsAsFactors = FALSE)
        # Assign color to each edge
        color_edge <- data.frame(value = as.character(round(edge_col, 2)), stringsAsFactors = FALSE) %>%
          dplyr::left_join(graph_col_df, by = "value") %>%
          dplyr::pull(color)
        pdf(file=paste0(sf,"24h.pval1e-3_spot10.cp.lr.nspot_width.number_color.igraph.circle.plot.pdf"))
        plot(graph_ntw,
             # Size of the edge
             edge.width = edge_importance*10, #*2 enhance the visualization
             edge.color = color_edge,
             # Size of the buble
             vertex.size = deg*2,
             vertex.color = "#E4791B",
             vertex.frame.color = "grey",
             vertex.label.color = "black",
             #vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
             layout = layout.circle)
        dev.off()

#remove Astro-Svep1
sn4<-sn4[sn4$from!="Astro.Svep1"&sn4$to!="Astro.Svep1",]
test<-sn4
#colnames(test)<-c('Var1','n','from','to')
test$n_spot<- scale(as.numeric(test$n_spot), center = 0)
        links <- data.frame(
            source = test$from,
            target = test$to,
            importance = as.numeric(test$n_spot)
          )
        #head(links)
        #nodes <- data.frame(name=unique(c(test$from,test$to)))
        #nodes
        graph_ntw <- igraph::graph_from_data_frame(links,
                vertices = nodes,
                directed = FALSE)
        
        deg <- degree(graph_ntw, mode="all")
        #pdf(paste0(result_dir,"spatial_celltype_withinspot_interaction_probability.circle.plot.png"))
        # Get color palette for difusion
        edge_importance <- E(graph_ntw)$importance #
        edge_col<-test$exp
        # Select a continuous palette
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
        
        # Create a color palette
        getPalette <- colorRampPalette(rev(c("#0C3321","#155637","#1a6a44","#1e7b4e","#2A6A63","#50817A","#6AA39A","#98BFB9","#C5DBD8")))
        #getPalette <- colorRampPalette(rev(c("#1D487C","#2864AD","#3d80d2","#6399db","#99bde7")))
        # Get how many values we need
       #grad_edge <- seq(3, 16
       #                 , 0.1)
        #grad_edge <- seq(11, 155, 0.1)
        grad_edge <- seq(0, 1, 0.01)
        
        # Generate extended gradient palette dataframe
        graph_col_df <- data.frame(value = as.character(grad_edge),
                                   color = getPalette(length(grad_edge)),
                                   stringsAsFactors = FALSE)
        # Assign color to each edge
        color_edge <- data.frame(value = as.character(round(edge_col, 2)), stringsAsFactors = FALSE) %>%
          dplyr::left_join(graph_col_df, by = "value") %>%
          dplyr::pull(color)
       pdf(file=paste0(sf,"72h.pval1e-3_spot10.cp.lr.nspot_width.number_color.igraph.circle.plot.pdf"))
        plot(graph_ntw,
             # Size of the edge
             edge.width = edge_importance*10, #*2 enhance the visualization
             edge.color = color_edge,
             # Size of the buble
             vertex.size = deg*4,
             vertex.color = "#E4791B",
             vertex.frame.color = "grey",
             vertex.label.color = "black",
             #vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
             layout = layout.circle)
       dev.off()


