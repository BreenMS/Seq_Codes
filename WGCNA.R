#THIS CODE IS A TEMPLATE FOR BUILDING WGCNA. THERE ARE SEVERAL WAYS ONE WOULD BUILD FROM THIS TO FIT A RANGE OF DIFFERENT EXPERIMENTS.
library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads(n=16)

########################################
########################################
##CREATE COMBINED CASE/CONTROL NETWORK##
########################################
#######################################

#INPUT NORMALIZED MATRIX
PRE <- read.delim("40_ASD_VOOM.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(PRE) =names(PRE)
dataExpr0<-as.data.frame(t(PRE))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

#INPUT META-DATA
traitData <- read.delim("Target.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH META-DATA TO EXPRESSION
rowsExpr <- rownames(dataExpr0)
test <- rownames(traitData$Subject.1)
traitRows <- match(rowsExpr,traitData$Subject.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

#CALCULATE SOFT-THRESHOLDS FOR EXPONENTIAL WEIGHTS
powers=c(1:20) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed") #SIGNED NETWORKS ARE PREFERRED OVER UNSIGNED NETWORKS
sft1=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="unsigned")


pdf("Soft_Threshold.pdf")
par(mfrow=c(2,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="steelblue4")
abline(h=0.8,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="steelblue4")
dev.off()

#plot(sft1$fitIndices[,1],-sign(sft1$fitIndices[,3])*sft1$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
#text(sft1$fitIndices[,1],-sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],labels=powers,col="steelblue4")
#abline(h=0.8,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
#plot(sft1$fitIndices[,1],sft1$fitIndices[,5],type="n",
#xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
#text(sft1$fitIndices[,1],sft1$fitIndices[,5],labels=powers,col="steelblue4")

#USING THE ABOVE OUTPUT, DESIGNATE THE IDEAL POWER FOR A WEIGHTED NETWORK AND CREATE A WEIGHTED GENE TREE FOR CUTTING (MODULE IDENTIFICATION)
adjacencyPre = adjacency((dataExpr0),power=5 ,type="signed") 
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = flashClust(as.dist(dissTOMPre), method="average")

#VISUALIZE TREE
pdf("GeneTree.pdf", w=11)
plot(geneTreePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity", labels=FALSE,hang=0.04, cex.axis=0.7);
dev.off()


#MODULE ASSIGNMENTS USING THE CUTREEHYBRID ALGORTIM (SET DEEP-SPLIT AND MODULE SIZE SPECS)
mColorh=NULL
for (ds in 0:4){
 tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
   minClusterSize = 15, cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMPre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

#MODULE ASSIGNMENTS USING THE CUTREEDYNAMIC ALGORTIM (SET DEEP-SPLIT AND MODULE SIZE SPECS)
DetectedColors = NULL;
DetectedColors = cbind(DetectedColors,labels2colors(cutreeDynamic(dendro = geneTreePre,
cutHeight = NULL, minClusterSize = 5,
method = "tree", deepSplit = TRUE)));

#COMPARE THE TREE CUTTING OUTPUTS
pdf("DeepSplit_Choices_hybrid.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Gene Co-Expression Network",dendroLabels=FALSE);
dev.off()

pdf("DeepSplit_Choices_dynamic.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, DetectedColors, paste("dpSplt =",0:4), main = "Gene Co-Expression Network",dendroLabels=FALSE);
dev.off()


#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,4]
table(modulesPRE)

modulesPRE =  DetectedColors[,1]
table(modulesPRE)

#Check to see if network modules can be cut and merged...
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
MEDiss2 = cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")
METree2 = hclust(as.dist(MEDiss2), method ="average")

pdf("Module_Relationships_hybrid.pdf")
plot(METree, main ="Clustering of Module Eigengenes",xlab ="",sub="")
abline(h=0.05,col="blue")
abline(h=0.1,col="red")
abline(h=0.15,col="blue")
abline(h=0.2,col="red")
plot(METree2, main ="Clustering of Module Eigengenes",xlab ="",sub="")
abline(h=0.05,col="blue")
abline(h=0.1,col="red")
abline(h=0.15,col="blue")
abline(h=0.2,col="red")
dev.off()

####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
MEDissThres = 0.2
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreePre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method ="average")
##########################################################################################################################


#CALCULATE MEs, PCA and TREES FOR VISUALIZATION OF MODULE RELATIONSHIPS
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
MDS_PD_n   = cmdscale(as.dist(distPCPD),3)
colorsPD = names(table(modulesPRE))
names = row.names((dataExpr0))

#Visualize module relationships
pdf("Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="", hang=0.3)
plot(MDS_PD_n, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesPRE)))
{
 par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
 plotMat(t(scale(dataExpr0[,modulesPRE==which.module])),
,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))

  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.5, font=2)
};
dev.off();




#HERE WE INPUT A VECTOR OF GENES WITH P-VALUES (I.E. DIFFERENTIAL EXPRESSION) TO DETERMINE MODULE ENRICHMENT WITH DEGS.
datSummary <- read.delim("Diff_Gene_P.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
GS_CM=-log10(datSummary$CM)
GS_CP=-log10(datSummary$CP)
GS_MP=-log10(datSummary$MP)
GS_Ave=-log10(datSummary$Average)
GS_Min=-log10(datSummary$Min)
whiteColor=datSummary$Color

#HERE WE VISUAL THE TREE WITH, COLORCODING THE DEGS
pdf("GeneTree_color.pdf",height=8,width=14)
par(mar=c(3,5,3,2))
plotDendroAndColors(geneTreePre, colors=modulesPRE, main="", dendroLabels=FALSE, hang=0.03, addGuide=FALSE, guideHang=0.05) 
dev.off()


#A MORE QUANTIATIVE WAY OF DOING THIS IS BY BARPLOTS AND KW TESTING (AS BELOW)
pdf("Module_MS_hybrid.pdf")
par(mfrow=c(3,2))
verboseBarplot(GS_CM,modulesPRE,color=colorsPD,main="CvM MS\n" ,xlab="",ylab="MS CvM (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    
dev.off()

#Multi-deminsional scaling plot of all data
cmd1=cmdscale(as.dist(dissTOMPre),2) #this may take a while
pdf("MDS_All_Modules.pdf")
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorsPD), main="MDS Plot", xlab="Scaling Dim. 1", ylab="Scaling Dim. 2")
dev.off()



colnames(ME_PD) =names(ME_PD)
dim(ME_PD)

#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs))
pdf("Modules_Heatmap.pdf", h=10, w=10)
plotEigengeneNetworks(MET,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

#INPUT META-DATA FOR CORRELATION ANALYSIS
traitData <- read.delim("Target.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH META-DATA TO EXPRESSION
rowsExpr <- rownames(dataExpr0)
test <- rownames(traitData$Subject.1)
traitRows <- match(rowsExpr,traitData$Subject.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

#CALCULATE ASSOCIATIONS BETWEEN EXPRESSION AND META-DATA
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf("ME_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8, 3, 1));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels=names(MEs),
ySymbols=names(MEs),
#yLabels = names(MEs0),
#ySymbols = paste (names(MEs), ": ", table (modulesPRE), sep=""),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-1,1),
cex.lab.x = 0.7,
cex.lab.y = 0.7,
main = paste("ME-Trait Relationships"))
dev.off()

#  CALCULATE MODULE MEMBERSHIP, P-VALUE SIGNIFICANT FOR EACH GENE IN EACH MODULE, AND HUB GENES AND OUTPUT ALL DATA INTO ONE MATRIX
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) 
nSamples = nrow(dataExpr0) 
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="Network_WGCNA_OUTPUT.csv")


#HERE WE EXPORT 20 TOP INTRAMODULAR HUB GENES PER MODULE OF INTEREST
modules = c("green") #"red", "yellow", "lightgreen") #,"yellow", "lightgreen", "black"
probes=names(dataExpr0)
inModule=is.finite(match(modulesPRE, modules));
modProbes=probes[inModule] ##
modTom=dissTOMPre[inModule, inModule]
dimnames(modTom)=list(modProbes, modProbes)
nTopHubs = 50
kIN = softConnectivity(dataExpr0[, modProbes])
selectHubs = (rank (-kIN) <= nTopHubs)

#EXPORT ENTIRE NETWORK
cyt=exportNetworkToCytoscape(modTom,
edgeFile=paste("CytoscapeInput_Pre_green_edges.txt", sep=""),
nodeFile=paste("CytoscapeInput_Pre_green_nodes.txt", sep=""),
weighted=TRUE,
threshold=0.50,
nodeNames=modProbes,
nodeAttr=modulesPRE[inModule])

##HERE WE PERFORM GENE-SET AND CELL-TYPE ENRICHMEN
Gene  = colnames(dataExpr0)
enrichments = userListEnrichment(Gene, modulesPRE,
fnIn = NULL,
catNmIn = NULL,
#fnIn = c("GeneList","ModuleColors"),
#catNmIn =  c("Genes","Modules"),
nameOut = "Network_Enrichment.csv", 
useBrainLists = TRUE,
useBloodAtlases = TRUE,
omitCategories = "grey", 
outputCorrectedPvalues = TRUE,
useStemCellLists = TRUE, 
outputGenes = TRUE, 
minGenesInCategory = 2, 
useBrainRegionMarkers = TRUE,
useImmunePathwayLists = TRUE,
usePalazzoloWang = TRUE)


