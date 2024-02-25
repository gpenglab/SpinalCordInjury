library(WGCNA)
library(dplyr)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
library(Seurat)
library(pheatmap)
library(scales)
library(colorspace)
library(edgeR)

### preprocess: get tissue-domain average gene expression matrix and dattrait2
#load SCT data 
da<-readRDS("WT.merge.replace_v2.SCT.regress_CC.nC.mt.ident.pc20.k50.res02.rds")
meta<-da@meta.data
sctda<-da@assays$SCT@data
sctda<-as.matrix(sctda)
#remove low detection genes
keep<-rowSums(sctda)>=100
table(keep)
#remove low detection genes
y<-DGEList(counts=sctda)
#filter cpm>1>=2的基因
keep<-rowSums(cpm(y)>1)>=20
sctda_filt<-sctda[keep,]
dim(sctda_filt)
meta$group<-paste0(meta$orig.ident,"_",meta$domain_res02)
datexpr<-data.frame(row.names = rownames(sctda_filt))
for(i in unique(meta$group)){
    cells<-rownames(subset(meta,group==i))
    datexpr[,i]<-apply(sctda_filt[,cells],1,mean)
}
saveRDS(datexpr,"WT.SCT.domain.thr100.datexpr.rds")
# transverse expression matrix
datexpr<-as.data.frame(t(as.matrix(datexpr)))

#1.b checking data for excessive missing data and identification of outlier samples
gsg=goodSamplesGenes(datexpr,verbose = 3)
gsg$allOK  #last statement returns TRUE means all genes pass the cuts
#cluster samples to see if there are any obvious outliers
sampletree=hclust(dist(datexpr),method = "average")
# Remove the offending genes and samples from the data:
datexpr<-datexpr[gsg$goodSamples,gsg$goodGenes]
options(repr.plot.width=10, repr.plot.height=7#,font.size=1.5
       )
pdf("WT.domain.SCT.data.thr100.average_sample.hclust.pdf",width = 10,height = 7)
#
plot(sampletree,main = "WT.SCT.domain_sampletree.hclustering to detect outliers",
     sub = "",xlab = "",cex=0.7#,cex.axis=1,cex.main=3
    )
dev.off()
# create dattrait
dattrait=data.frame(row.names =rownames(datexpr))
for(i in 1:nrow(datexpr)){
    dattrait[,rownames(datexpr)[i]]=c(rep(0,i-1),1,rep(0,nrow(datexpr)-i))
}
dattrait
save(datexpr,dattrait,file = "WT.SCT.domain.thr100-01-datainput.RData")
#load data saved in the first part
#rm(list=ls())
lnames=load(file = "WT.SCT.domain.thr100-01-datainput.RData")
#the variable lnames contains the name of loaded variables
#lnames
#2.a.1 choosing the soft-thresholding power(beta):analysis of network topology
#co-expression similarity is raised to calculate adjacency
#choose a set of soft-thresholding powers
powers=c(1:30)
##call the network topology analysis function
sft=pickSoftThreshold(datexpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed")  

### pick a soft threshold power near the curve of the plot
cex1=0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
softpower=18
adjacency=adjacency(datexpr,power = softpower,type = "signed")
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "signed")
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average"
                 )
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#pdf("Gene clustering on TOM-based dissimilarity.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
#dev.off()

abline(h = 0.99,col="red")

minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, #distM = dissTOM,
                            method = "tree",
                deepSplit = 2,pamRespectsDendro = FALSE,
                #minSplitHeight=0.3,
                #cutHeight=0.99,
                minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(30,20)
options(repr.plot.width=10,repr.plot.height=7)
#pdf("Dynamic Tree Cut.pdf",width = 10,height = 7)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
abline(h = 0.99,col="red")
#dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datexpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
#sizeGrWindow(7, 6)
#pdf("Clustering of module eigengenes.pdf",width = 10,height = 7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.15
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
#abline(h=0.1,col="blue")
#dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datexpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "WT.SCT.domain-02-networkConstruction-stepByStep.RData")

#pdf("geneDendro-mergheight01.pdf", width = 10, height = 7)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# output modules
WGCNA.modules.full.gene.list = data.frame(module = moduleColors,
                                          gene = colnames(datexpr))
write.csv(WGCNA.modules.full.gene.list, file = "thr100.sp18.hclust_average.min30.deep2.h015.34modules.full.gene.list.csv",quote=FALSE)

WGCNA.modules.full.gene.list<-read.csv("thr100.sp18.hclust_average.min30.deep2.h015.34modules.full.gene.list.csv")
head(WGCNA.modules.full.gene.list)
table(WGCNA.modules.full.gene.list$module=="grey")

#### TOMplot
#transfer dissTOM with a power to make moderately strong, connections more visible in the heatmap
plotTOM<-dissTOM^10
#set diagonal to NA for a nicer plot
diag(plotTOM)<-NA

myheatcol<-rev(viridis_pal(option = "D")(10))#rev(c('#91BFDB','#FEE090','#FC8D59','#D73027'))
myheatcol

png("thr100.sp18.min30.deep2.h015.TOMplot_220720.col2.png")
TOMplot(plotTOM,geneTree,
        moduleColors,#dynamicColors#,
        col=myheatcol
       )
dev.off()

# Define numbers of genes and samples
nGenes = ncol(datexpr)
nSamples = nrow(datexpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datexpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
head(MEs)
## merge distance MEs together
mMEs<-MEs
mMEs$domain<-factor(str_sub(rownames(mMEs),-2,-1),levels=c("WM","MG","DH","VH"))
mMEs$time<-factor(str_split(rownames(MEs),"_[H*T]_",simplify = T)[,1],
                   levels = c('WT_sham','WT_3h','WT_24h','WT_72h'))

mMEs$group<-paste0(mMEs$time,"_",mMEs$domain)
#mMEs<-arrange(mMEs,group_by=domain)
#head(mMEs)
mMEs<-mMEs[,-which(colnames(mMEs)%in%c("domain","time"))]
head(mMEs)

test<-as.matrix(mMEs[,-which(colnames(mMEs)=="group")])
head(test)
group<-mMEs$group
names(group)<-rownames(mMEs)
head(group)

mergedMEs<-as.data.frame(rowsum(test,group = group))
head(mergedMEs)

mergedMEs<-mergedMEs[c('WT_sham_DH','WT_3h_DH','WT_24h_DH','WT_72h_DH',
                                  'WT_sham_MG','WT_3h_MG','WT_24h_MG','WT_72h_MG',
                                  'WT_sham_VH','WT_3h_VH','WT_24h_VH','WT_72h_VH',
                                  'WT_sham_WM','WT_3h_WM','WT_24h_WM','WT_72h_WM'),]

dattrait2<-data.frame(row.names = rownames(mergedMEs))
j=1
for(i in rownames(mergedMEs)){
    dattrait2[,i]<-c(rep(0,j-1),1,rep(0,nrow(dattrait2)-j))
    j=j+1
}
dattrait2

moduleTraitCor = cor(mergedMEs, dattrait2, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 16)

# Will display correlations and their p-values
pdf("module-time_domain relationships.pdf",width=30,height=30)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(13,12, 3, 3));
# Display the correlation values within a heatmap plot
options(repr.plot.width=30,repr.plot.height=30)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dattrait2),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1,1),
               plotLegend=T,
               main = paste("Module-trait relationships"))
dev.off()

row_anno<-data.frame(row.names = rownames(moduleTraitCor),"module"=gsub("ME","",rownames(moduleTraitCor)))
row_anno
row_color<-row_anno$module
names(row_color)<-row_anno$module
col_anno<-data.frame(row.names = colnames(moduleTraitCor),"domain"=c(rep("DH",4),rep("MG",4),rep("VH",4),rep("WM",4)),
                    "time"=c(rep(c("sham","3h","24h","72h"),4)))
col_anno
domain_col<-c('#BC3C29A8','#0072B5A8','#E18727A8','#20854EA8')
names(domain_col)<-c("DH", "MG", "VH", "WM")
time_col<-c('#374E55FF','#DF8F44FF','#00A1D5FF','#B24745FF')
names(time_col)<-c("sham","3h","24h","72h")
col_list<-list(module=row_color,domain=domain_col,time=time_col)
gap_col<-c(4,8,12)
p<-pheatmap(moduleTraitCor,#col=colorRampPalette(rev(brewer.pal(10,"Spectral")))(50),
         scale="none",
          cluster_cols=F,cluster_rows=F,#treeheight_row=10,cutree_rows=7,
         annotation_row=row_anno,
         annotation_col=col_anno,
         border_color="NA",
         fontsize=5,
         annotation_colors=col_list,
         legend=T,gaps_col=gap_col,
         #cellwidth=25,cellheight=25,
         show_rownames=T)
save_pheatmap_pdf(p,filename = "module-time_domain relationships.pheatmap.pdf")

### merge time into unique domain dattrait
mMEs<-MEs
mMEs$domain<-factor(str_sub(rownames(mMEs),-2,-1),levels=c("WM","MG","DH","VH"))
test<-as.matrix(mMEs[,-which(colnames(mMEs)=="domain")])
group<-mMEs$domain
names(group)<-rownames(mMEs)
mergedMEs<-as.data.frame(rowsum(test,group = group))
dattrait2<-data.frame(row.names = rownames(mergedMEs))
j=1
for(i in rownames(mergedMEs)){
    dattrait2[,i]<-c(rep(0,j-1),1,rep(0,nrow(dattrait2)-j))
    j=j+1
}
dattrait2
moduleTraitCor = cor(mergedMEs, dattrait2, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 4)
# Will display correlations and their p-values
pdf("module-domain relationships.withoutText.pdf",width=8,height=30)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(13,12, 3, 3));
# Display the correlation values within a heatmap plot
options(repr.plot.width=8,repr.plot.height=30)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dattrait2),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1,1),
               plotLegend=T,
               main = paste("Module-trait relationships"))
dev.off()

###export network
df<-WGCNA.modules.full.gene.list
df<-df[!df$module=="grey",]
probes=colnames(datexpr)
#select modules
modules=unique(df$module)
## export module one by one
sf<-"cytoscape/"
if(!dir.exists(sf))
    dir.create(sf)

for(i in 1:length(modules)){
    module=modules[i]
    inModule=(mergedColors==module)#is.finite(match(dynamicColors,i))
    head(inModule)
    modGenes=probes[inModule]
    length(modGenes)
    # select the corresponding Topological Overlap
    modTOM=TOM[inModule,inModule]
    dim(modTOM)
    dimnames(modTOM)=list(modGenes,modGenes)
    modTOM[1:5,1:5]
    # Export the nerwork into edge and node list files Cytoscape can read
    cyt=exportNetworkToCytoscape(modTOM,
                            edgeFile=paste(sf,"CytoscapeInput-edges-",paste(module,collapse = "-"),".txt",sep = ""),
                            nodeFile=paste(sf,"CytoscapeInput-nodes-",paste(module,collapse = "-"),".txt",sep = ""),
                            weighted=TRUE,
                            threshold=0.02,
                            nodeNames=modGenes,
                            altNodeNames=modGenes,
                            nodeAttr=dynamicColors[inModule])
    top100<-cyt$edgeData %>% top_n(.,wt=weight,n=100)
    write.csv(top100,paste(sf,"CytoscapeInput-edges-",paste(module,collapse = "-"),"wt.top100.csv",sep = ""),row.names=F,quote=F)
    top50<-cyt$edgeData %>% top_n(.,wt=weight,n=50)
    write.csv(top50,paste(sf,"CytoscapeInput-edges-",paste(module,collapse = "-"),"wt.top50.csv",sep = ""),row.names=F,quote=F)
}

### adjust cyan module(remove Rpl*, Rps* related genes)
module="cyan"
    inModule=which(WGCNA.modules.full.gene.list$module==module)#is.finite(match(dynamicColors,i))
    head(inModule)
    modGenes=probes[inModule]
    length(modGenes)
    # select the corresponding Topological Overlap
    modTOM=TOM[inModule,inModule]
    dim(modTOM)
    dimnames(modTOM)=list(modGenes,modGenes)
    modTOM[1:5,1:5]    

modTOM<-modTOM[!grepl("Rpl",rownames(modTOM)),!grepl("Rpl",rownames(modTOM))]
modTOM<-modTOM[!grepl("Rps",rownames(modTOM)),!grepl("Rps",rownames(modTOM))]
modGenes=rownames(modTOM)
# Export the nerwork into edge and node list files Cytoscape can read
cyt=exportNetworkToCytoscape(modTOM,
                            edgeFile=paste("cytoscape/","CytoscapeInput-edges-",paste(module,collapse = "-"),".txt",sep = ""),
                            nodeFile=paste("cytoscape/","CytoscapeInput-nodes-",paste(module,collapse = "-"),".txt",sep = ""),
                            weighted=TRUE,
                            threshold=0.02,
                            nodeNames=modGenes,
                            altNodeNames=modGenes#,
                            #nodeAttr=WGCNA.modules.full.gene.list$module[inModule]
                            )

    top50<-cyt$edgeData %>% top_n(.,wt=weight,n=50)
write.csv(top50,paste("cytoscape/","CytoscapeInput-edges-",paste("cyan",collapse = "-"),"wt.top50_rmRp.220812.csv",sep = ""),row.names=F,quote=F)












