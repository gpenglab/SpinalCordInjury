#library required packages
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(igraph)
library(stringr)

#load 10x data
da<-readRDS("sci.rds")
meta<-read.csv("bs_sci.csv")
rownames(meta)<-meta[,1]
meta<-meta[,-1]
meta$celltype2<-da$celltype
n_cell<-rownames(meta)[meta$celltype2=="Neuron"]
l3_na<-rownames(meta)[is.na(meta$L3_taxon)]
meta[n_cell,"L3_taxon"]<-"Neuron"
meta[n_cell,"L2_taxon"]<-"Neuron"

gfap_state<-readRDS("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40_gfaphclust.all.rds")
me<-gfap_state@meta.data
me$gfap_hclust[is.na(me$gfap_hclust)]<-"Astrocyte.Slc7a10"
me$gfap_hclust[me$gfap_hclust==1]<-"WM_Gfap"
me$gfap_hclust[me$gfap_hclust==2]<-"GM_Gfap"

c_slc<-rownames(me)[me$subtype=="Astr_Slc7a10"]
c_gfap<-rownames(me)[me$subtype=="Astr_Gfap"]

c_wm<-rownames(me)[me$gfap_hclust=="WM_Gfap"]
c_gm_gfap<-rownames(me)[me$gfap_hclust=="GM_Gfap"]

meta$subtype_0624<-meta$L3_taxon
meta[c_gfap,"subtype_0624"]<-"Astrocyte.Gfap"
meta[c_slc,"subtype_0624"]<-"Astrocyte.Slc7a10"

#remove redundant Astrocyte
meta_sub<-meta[meta$subtype_0624!="Astrocyte",]

meta_sub$gfap_hclust_0624<-meta_sub$subtype_0624
meta_sub[c_wm,"gfap_hclust_0624"]<-"WM_Gfap"
meta_sub[c_gm_gfap,"gfap_hclust_0624"]<-"GM_Gfap"

### remove NA
mete_sub<-meta_sub[!is.na(meta_sub$gfap_hclust_0624),]

da_ad<-da[rownames(da@assays$RNA@counts),rownames(meta_sub)]
da_ad@meta.data<-meta_sub[rownames(da_ad@meta.data),]
saveRDS(da_ad,"sci.astr_hclust.rds")

#downsample to 1000
Idents(da_ad)<-da_ad$gfap_hclust_0624
set.seed(220624)
da_ad_s<-subset(da_ad,downsample=1000)
saveRDS(da_ad_s,"sci.astr_hclust.down1000.rds")
c_time<-rownames(da_ad_s@meta.data)[da_ad_s$time!="7dpi"]
### remove 7 dpi
f_da<-da_ad_s[rownames(da_ad_s@assays$RNA@counts),c_time]
ma<-as.matrix(f_da@assays$RNA@data)
saveRDS(ma,"10Xsc.hclust_astr.down1000.without7dpi.RNA.data.rds")

meta<-f_da@meta.data
meta$gfap_hclust_0624<-as.character(meta$gfap_hclust_0624)

meta$celltype<-meta$gfap_hclust_0624

### merged into 16 cell types
meta[meta$celltype %in% c("A-Endothelial","C-Endothelial",'Pericyte','Tip Cell','U-Vascular','V-Endothelial','VSMC'),"celltype"]<-"Vascular"
meta[meta$celltype %in% c('Div-OPC','OPC-A','OPC-B','Pre-Oligo'),"celltype"]<-"OPC"
meta[meta$celltype %in% c('Border-Associated Mac','Chemotaxis-Inducing Mac','Inflammatory Mac'),"celltype"]<-"Macrophage"
meta[meta$celltype %in% c('Dividing Microglia','Inflammatory Microglia','Migrating Microglia'),"celltype"]<-"R.Microglia"
meta[meta$celltype %in% c('Dividing Myeloid','Interferon Myeloid'),"celltype"]<-"Myeloid"
meta[meta$celltype %in% c('Astroependymal','Ependymal-A','Ependymal-B'),"celltype"]<-"Ependymal"

###change celltype names to fit our data
meta$celltype<-gsub("_",".",gsub(" ",".",meta$celltype))
meta$celltype<-gsub("_",".",gsub(" ",".",meta$celltype))

write.csv(meta,"10Xsc.hclust_astr.down1000.without7dpi.merged16celltype.meta.csv")

cell1<-rownames(meta[meta$time=="Uninjured",])
cell2<-rownames(meta[meta$time=="1dpi",])
cell3<-rownames(meta[meta$time=="3dpi",])

unj_s3_ma<-ma[,cell1]
unj_s3_me<-meta[cell1,]
d1_s3_ma<-ma[,cell2]
d1_s3_me<-meta[cell2,]
d3_s2_ma<-ma[,cell3]
d3_s2_me<-meta[cell3,]

me_list<-list(unj_s3_me,d1_s3_me,d3_s2_me)
names(me_list)<-c("uninj","d1","d3")

ma_list<-list(unj_s3_ma,d1_s3_ma,d3_s2_ma)
names(ma_list)<-c("uninj","d1","d3")
#head(ma_list)

###extract cellchat LR interaction to our data
CellChatDB<-CellChatDB.mouse
interaction_input<-CellChatDB$interaction
head(interaction_input)

#geneInfo<-CellChatDB$geneInfo
#head(geneInfo)

cc_lr<-as.data.frame(str_split(interaction_input$interaction_name_2," - ",simplify = T))
cc_r<-as.data.frame(str_split(gsub("\\)","",gsub("\\(","",cc_lr$V2)),"\\+",simplify = T))
head(cc_lr)
head(cc_r)
#colnames(cc_r)<-c("R1","R2")
#cc_r$R1<-gsub("[[:space:]]","",cc_r$R1)
#cc_r$R2<-gsub("[[:space:]]","",cc_r$R2,fixed = TRUE)
cc_lr<-cbind(cc_lr[,1],cc_r)
cc_lr1<-cc_lr[,1:2]
colnames(cc_lr1)<-c("ligand","receptor")
dim(cc_lr1)
head(cc_lr1)
cc_lr2<-cc_lr[!cc_lr[,3]=="",c(1,3)]
colnames(cc_lr2)<-c("ligand","receptor")
dim(cc_lr2)
head(cc_lr2)
cc_lr<-rbind(cc_lr1,cc_lr2)
#head(cc_lr)
#colnames(cc_lr)<-c("ligand","receptor")
cc_lr$LR<-paste0(cc_lr$ligand,"_",cc_lr$receptor)
#colnames(cc_lr)[1]<-"L"
#cc_lr<-cc_lr[,c("LR","ligand","receptor")]
cc_lr<-apply(cc_lr,2,function(x) gsub("[[:space:]]","",x))
#rownames(cc_lr)<-cc_lr$LR
cc_lr<-cc_lr[,c("LR","ligand","receptor")]
head(cc_lr)
write.csv(cc_lr,"cellchat.LR.df.csv")
#saveRDS(cc_lr_list,"cellchat.LR.list.rds")         

#calculate LR interaction at each timepoint
for(i in names(ma_list)){
    temp<-ma_list[[i]]
    temp2<-me_list[[i]]
    # subset cell type
    if(i=='uninj')
        cells<-rownames(temp2)[temp2$celltype %in% c('WM.Gfap','Homeostatic.Microglia',
                                                     'Astrocyte.Slc7a10','Fibroblast','Ependymal',
                                                     'Oligodendrocyte','Neuron','OPC','Vascular')]
    if(i=="d1")
        cells<-rownames(temp2)[temp2$celltype %in% c('WM.Gfap','GM.Gfap',
                                                     'Astrocyte.Slc7a10','Fibroblast','Neutrophil',
                                                     'Oligodendrocyte','Neuron','OPC','Vascular',
                                                      'Monocyte','R.Microglia','Myeloid','Macrophage'
                                                     )]
    if(i=="d3")
         cells<-rownames(temp2)[temp2$celltype %in% c('WM.Gfap','GM.Gfap','Homeostatic.Microglia',
                                                     'Astrocyte.Slc7a10','Fibroblast','Astro.Svep1',
                                                     'Oligodendrocyte','Neuron','OPC','Vascular',
                                                      'R.Microglia','Myeloid','Macrophage'
                                                     )]
    # subset gene
    fil <- apply(temp, 1, function(x) length(x[x>1])>= 10)
    temp<-as.matrix(temp[fil,cells])
    temp2<-temp2[cells,]
                 
    #create a cellchat object
    cellchat<-createCellChat(object = temp,meta = temp2,group.by = "celltype")
    
    #Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <-CellChatDB#subsetDB(CellChatDB, search = c("ligand","receptor")) # use Secreted Signaling
    
    cellchat@DB<-CellChatDB.use
    
    #Preprocessing the expression data for cell-cell communication analysis
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multiprocess", workers = 4) # do parallel
    
    cellchat <- identifyOverExpressedGenes(cellchat,thresh.fc = 0.5,thresh.pc = 0.5)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    computeAveExpr(cellchat, features = c("Spp1","Cd44"))
    
    #assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations
    #consider the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE.
    
    #Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat,population.size = F)
    
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10) 
    
    #Extract the inferred cellular communication network as a data frame
    
    df.net <- subsetCommunication(cellchat) #a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
    #df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) # gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
    #df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
    sf1<-"fc05.pc05.min10.interaction.net/"
    if(!dir.exists(sf1))
        dir.create(sf1)
    write.csv(df.net,paste0(i,".fc05.pc05.min10.interaction.net.csv"))
}

