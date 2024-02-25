library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)

#1.Load spatial position and deconvolution meta
meta<-read.csv("16celltypes.dec_conf.nor.meta.csv")
rownames(meta)<-meta[,1]
meta<-meta[,-1]

#2.celltype pair df
cp1<-read.csv(paste0(cpdir,"WT_sham.binaried.celltype_pair.csv"))
rownames(cp1)<-cp1[,1]
cp1<-cp1[,-1]
cp2<-read.csv(paste0(cpdir,"WT_3h.binaried.celltype_pair.csv"))
rownames(cp2)<-cp2[,1]
cp2<-cp2[,-1]
cp3<-read.csv(paste0(cpdir,"WT_24h.binaried.celltype_pair.csv"))
rownames(cp3)<-cp3[,1]
cp3<-cp3[,-1]
cp4<-read.csv(paste0(cpdir,"WT_72h.binaried.celltype_pair.csv"))
rownames(cp4)<-cp4[,1]
cp4<-cp4[,-1]

cn<-unique(c(colnames(cp1),colnames(cp2),colnames(cp3),colnames(cp4)))
length(cn)

test<-cp1
test[,setdiff(cn,colnames(cp1))]<-0
cp2[,setdiff(cn,colnames(cp2))]<-0
cp3[,setdiff(cn,colnames(cp3))]<-0
cp4[,setdiff(cn,colnames(cp4))]<-0
cp<-Reduce(rbind,list(test[,cn],cp2[,cn],cp3[,cn],cp4[,cn]))
cp<-apply(cp,2,as.factor)
me<-cbind(meta[rownames(cp),],cp)
write.csv(me,"EM.spr-04.ct.binaried.celltype_pair.meta.csv")

###CP aggregated proportion
cp_em_df<-data.frame(row.names = rownames(meta))
for(i in cn){
    ct1<-str_split(i,"_",simplify = T)[,1]
    ct2<-str_split(i,"_",simplify = T)[,2]
    tmp<-meta[,c(ct1,ct2)]
    cp_em_df[,i]<-apply(tmp,1,function(x) exp(mean(log(x))))
}
me_cp_em<-cbind(meta,cp_em_df)
write.csv(me_cp_em,"celltype_pair.EM_proportion.meta_220826.csv")

###Visualize

#pdf("Colocalize.spatial.Astrocyte.Gfap_Neuron_em_221027.pdf",width = 20,height = 16)
options(repr.plot.width=20,repr.plot.height=16)
    ggplot(me_cp_em,aes_string(x='x.ad',y='y.ad',color="Astrocyte.Gfap_Neuron"))+
        geom_point(size=1.2)+
        #scale_color_manual(values = rev(c("orange","#911eb4","#e6beff","blue","grey")))+
        theme(panel.background = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())+
        scale_color_gradientn(colours =c("#d6d4e0","#4B0082"))
        #scale_color_manual(values = c('#f2e394','#d5f4e6','#80ced6','#eea29a','#f9d5e5'#'#5b9aa0'#
         #                             ,'#618685'#,
                                      
         #                             ,'#d6cbd3','#a2836e'
         #                             ))
    ggsave("Colocalize.spatial.Astrocyte.Gfap_Neuron_em_221027.png",width = 20,height = 16,dpi = 400)
#dev.off()

p_df<-read.csv("4times.merged.spmc-04.pval005_jc005_spot15.lr.cp.long.df_sc.csv")
rownames(p_df)<-p_df[,1]
p_df<-p_df[,-1]

cc_lr<-read.csv("cellchat.all.lr.interaction.df.csv")
rownames(cc_lr)<-cc_lr[,1]
cc_lr<-cc_lr[,-1]
lrgenes<-unique(c(cc_lr$ligand,cc_lr$receptor1,cc_lr$receptor2))
lrgenes<-intersect(lrgenes,rownames(lr_exp))
## lr single express
lr_exp<-lr_exp[lrgenes,]
### lr-pair GM express
test<-as.data.frame(t(lr_exp))
for(i in 1:nrow(cc_lr)){
        lr<-as.character(cc_lr[i,][cc_lr[i,]!=""])
        nm<-rownames(cc_lr)[i]
        if(all(lr %in% lrgenes)){
            #expression
            lr1_exp<-test[,lr]
            head(lr1_exp)
            
            #the spot lr1 geometric mean exp
            sp_exp<-as.data.frame(apply(lr1_exp,1,function(x) exp(mean(log(x)))))
            colnames(sp_exp)<-nm
            test<-cbind(test,sp_exp)
                                        }
        }
saveRDS(test,"cellchat.lr.express.and.gm_express.rds")

me_exp<-cbind(test[rownames(meta),],meta)
saveRDS(me_exp,"cellchat.lr.express.and.gm_express.meta.rds")

#Line chart
temp<-cbind(me,me_exp[rownames(me),])
temp2<-temp[,c("time","Astrocyte.Gfap_Neuron","Astrocyte.Slc7a10_Neuron","NRXN1_NLGN2","NRXN1_NLGN3","NRXN3_NLGN2","NRXN3_NLGN3")]
temp2$cp<-""
temp2$cp[temp2$Astrocyte.Gfap_Neuron==1]<-"Astro.WM_Neuron"
temp2$cp[temp2$Astrocyte.Slc7a10_Neuron==1]<-"Astro.GM_Neuron"
temp3<-temp2[temp2$cp!="",]
temp3$cp<-factor(temp3$cp,levels = c("Astro.GM_Neuron","Astro.WM_Neuron"))
temp3$time<-factor(temp3$time,levels = c("WT_sham","WT_3h","WT_24h","WT_72h"))
temp3["row1",]<-c("WT_sham","1","0",0,0,0,0,"Astro.WM_Neuron")
temp3["row2",]<-c("WT_3h","1","0",0,0,0,0,"Astro.WM_Neuron")
temp3["row3",]<-c("WT_72h","0","1",0,0,0,0,"Astro.GM_Neuron")
temp3$group<-paste0(temp3$time,"-",temp3$cp)
temp3$NRXN1_NLGN2<-as.numeric(temp3$NRXN1_NLGN2)
temp3$NRXN1_NLGN3<-as.numeric(temp3$NRXN1_NLGN3)
temp3$NRXN3_NLGN2<-as.numeric(temp3$NRXN3_NLGN2)
temp3$NRXN3_NLGN3<-as.numeric(temp3$NRXN3_NLGN3)
temp_mean<-aggregate(temp3[,4:7],list(temp3$group),mean)
temp_mean<-as.data.frame(temp_mean)
temp_mean$time<-str_split(temp_mean$`Group.1`,"-",simplify = TRUE)[,1]
temp_mean$cp<-str_split(temp_mean$`Group.1`,"-",simplify = TRUE)[,2]
temp_mean$time<-factor(temp_mean$time,levels = rev(c("WT_sham","WT_3h","WT_24h","WT_72h")))
pdf("Astro-Neuron.NRXN3_NLGN3.mean.lineplot.pdf",width = 4,height = 10)
options(repr.plot.width=4,repr.plot.height=10)
ggplot(temp_mean, aes(x=time, y=NRXN3_NLGN3, group=cp)) + 
    geom_line(aes(col=cp))+ 
    scale_color_manual(values = c("#C71585","#4B0082"))+
    geom_point()+
    coord_flip()+
    theme_bw(base_line_size = 0.25)
#scale_linetype_manual(values=c("twodash", "dotted"))
dev.off()
