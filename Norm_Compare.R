#################################################################
################# Load the Libraries ############################
#################################################################
library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(rgl)
library(edgeR)
library(RUVSeq)
library(EDASeq)
#################################################################
################# UPPER QUARTILE NORM. ##########################
#################################################################

rawdata <- read.delim("Test.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)
dim(rawdata)

filter<-apply(rawdata,1, function(x) length(x[x>1])>=10) #require more than 5 reads in at least 10 samples
filtered <- rawdata[filter,]
dim(filtered)

targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor((targets$Group))

# create raw expression eset object
set<-newSeqExpressionSet(as.matrix(filtered),
phenoData=data.frame(Group,row.names=colnames(filtered)))

# create upper quantile normalization (between lane) eset object
set_uq <- betweenLaneNormalization(set,which="upper")

# extract and transpose expression matrix
upper_quartile_matrix <- normCounts(set_uq)
upper_quartile_matrix_t =t(upper_quartile_matrix)
write.table(upper_quartile_matrix, "Blood_Matrix_UQ.txt")

# Color code groups
colors <- c("firebrick1","brown","steelblue4", "steelblue2")
colnames(pData(set_uq))
table(pData(set)[,"Group"]);
table(pData(set_uq)[,"Group"]);

# create diagnostic plots for UQ normalization
par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Expression")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
#boxplot(normCounts(set_uq),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (UQ)", outline=FALSE,ylab="UQ Expression", ylim=c(0,200))
#legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_uq,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (UQ)", outline=FALSE, ylab="log UQ Expression",ylim=c(0,9))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
hist(log(counts(set_uq)), main="Histogram of log Norm UQ Counts", xlab="log UQ Expression")
plotRLE(set_uq,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Norm Counts UQ")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_uq,cex=0.6,col=colors[Group], main="PCA UQ", ylim=c(-0.3, 0.5), xlim=c(-0.25, 0.25), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plot(hclust(dist(upper_quartile_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set_uq,log=TRUE,ylim=c(0,20), main="UQ Over-Dispersion")





#################################################################
################# RUVg normalize to ERCC ########################
#################################################################
rawdata <- read.delim("40_ASD_Counts.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)
dim(rawdata)

filter<-apply(rawdata,1, function(x) length(x[x>5])>=10) #require more than 5 reads in at least 10 samples
filtered <- rawdata[filter,]
dim(filtered)

spikes <-rownames(filtered)[grep("ERCC",rownames(filtered))]

targets <- read.delim("Target.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor((targets$Group))

# create raw expression eset object
set<-newSeqExpressionSet(as.matrix(filtered),
phenoData=data.frame(Group,row.names=colnames(filtered)))

set_ercc <- RUVg (set, spikes, k=1)
pData (set_ercc) #factors of unwanted variation

# extract and transpose expression matrix
ercc_matrix <- normCounts(set_ercc)
ercc_matrix_t =t(ercc_matrix)
write.table(ercc_matrix, "40_ASD_RUVg_ERCC.txt")

# Color code groups
colors <- c("firebrick1","brown","steelblue4", "steelblue2")
colnames(pData(set_ercc))
table(pData(set)[,"Group"]);
table(pData(set_ercc)[,"Group"]);

par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Expression")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
#boxplot(normCounts(set_ercc),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (RUVg ERCC)", outline=FALSE,ylab="RUVg ERCC Expression", ylim=c(0,200))
#legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_ercc,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (RUVg ERCC)", outline=FALSE, ylab="log RUVg ERCC Expression",ylim=c(0,9))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
hist(log(counts(set_ercc)), main="Histogram of log Norm RUVg ERCC Counts", xlab="log RUVg ERCC Expression")
plotRLE(set_ercc,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Norm Counts RUVg ERCC")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_ercc,cex=0.6,col=colors[Group], main="PCA RUVg ERCC", ylim=c(-0.35, 0.55), xlim=c(-0.3, 0.3), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plot(hclust(dist(ercc_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set_ercc,log=TRUE,ylim=c(0,20), main="RUVg Over-Dispersion")


###############################################################
################ VOOM & TMM NORMALIZATION #####################
###############################################################
rawdata <- read.delim("HPC_Matrix.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(rawdata)
dim(rawdata)

y <- DGEList(counts=rawdata[,2:37], genes=rawdata[,1:1])
keep <- rowSums(cpm(y)>1) >= 10 # require 1 cpm in at least 10 samples
y <- y[keep,]
dim (y)

#apply scale normalization
y <- calcNormFactors(y)
y$samples

################
##### VOOM #####
################
VST <- voom(y, design = NULL, plot=TRUE)
voom_matrix <- cbind(VST$genes, VST$E)
write.table (voom_matrix, "HPC_Matrix_voom.txt", sep="\t")

exprsFile <- file.path("HPC_Matrix_voom.txt")
voom_matrix <- as.matrix(read.table(exprsFile, header=TRUE, sep="", row.names=1, as.is=TRUE))
dim(voom_matrix)
voom_matrix_t=t(voom_matrix)

#input meta-data
targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor((targets$Group))

#create e-set
set_voom<-newSeqExpressionSet(as.matrix(voom_matrix),
phenoData=data.frame(Group,row.names=colnames(voom_matrix)))

#draw QC plots
par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Expression")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
#boxplot(normCounts(set_voom),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (voom)", outline=FALSE,ylab="voom Expression", ylim=c(0,200))
#legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_voom,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (voom)", outline=FALSE, ylab="log voom Expression",ylim=c(-2,4))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
hist(log(counts(set_voom)), main="Histogram of log Norm voom Counts", xlab="log voom Expression")
plotRLE(set_voom,ylim=c(-1,1), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Norm Counts voom")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotMDS(voom_matrix,cex=0.8,col=colors[Group],main="PCA voom",ylim=c(-3, 3), xlim=c(-4, 4))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
#plotPCA(set_voom,cex=0.6,col=colors[Group], main="PCA voom", ylim=c(-0.35, 0.55), xlim=c(-0.3, 0.3), las=2)
plot(hclust(dist(voom_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set_voom,log=TRUE,ylim=c(0,20), main="VOOM Over-Dispersion")



library(limma)
exprsFile <- file.path("AMG_Matrix_voom.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

#create design matrix for DGE
design <- model.matrix(~0+Group)
colnames(design) <- c("Control_BNSS", "Control_F344", "Control_LEW", "Control_WKY", "CRS_BNSS", "CRS_F344", "CRS_LEW", "CRS_WKY", "PCRS_BNSS", "PCRS_F344", "PCRS_LEW", "PCRS_WKY")
colnames(design)
fit <- lmFit(exprs, design)

#specify contrasts of interest
contrast.matrix <- makeContrasts(
WKY=Control_WKY-PCRS_WKY,
F344=Control_F344-PCRS_F344,
LEW=Control_LEW-PCRS_LEW,
BNSS=Control_BNSS-PCRS_BNSS,
levels=design)

#apply linear model and output results
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef="WKY", adjust="BH")
WKY_DEG <- topTable(fit2, coef="WKY", adjust="BH",n=24000)
write.table(WKY_DEG, file="AMG_DEG_WKY.txt", quote=FALSE, row.names=TRUE)

topTable(fit2, coef="F344", adjust="BH")
F344_DEG <- topTable(fit2, coef="F344", adjust="BH",n=24000)
write.table(F344_DEG, file="AMG_DEG_F344.txt", quote=FALSE, row.names=TRUE)

topTable(fit2, coef="LEW", adjust="BH")
LEW_DEG <- topTable(fit2, coef="LEW", adjust="BH",n=24000)
write.table(LEW_DEG, file="AMG_DEG_LEW.txt", quote=FALSE, row.names=TRUE)

topTable(fit2, coef="BNSS", adjust="BH")
BNSS_DEG <- topTable(fit2, coef="BNSS", adjust="BH",n=24000)
write.table(BNSS_DEG, file="AMG_DEG_BNSS.txt", quote=FALSE, row.names=TRUE)



###############
##### TMM #####
###############
scale <- y$samples$lib.size*y$samples$norm.factors
X.tmm <- t(t(y$counts)/scale)*mean(scale)
head(X.tmm)
head(X.tmm)
TMM_Matrix <- cbind(y$genes, X.tmm)
write.table (TMM_Matrix, "40_ASD_TMM.txt", sep="\t")

exprsFile <- file.path("40_ASD_TMM.txt")
tmm_matrix <- as.matrix(read.table(exprsFile, header=TRUE, sep="", row.names=1, as.is=TRUE))
dim(tmm_matrix)
tmm_matrix_t=t(tmm_matrix)

set_tmm<-newSeqExpressionSet(as.matrix(tmm_matrix),
phenoData=data.frame(Group,row.names=colnames(tmm_matrix)))

par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Expression")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
#boxplot(normCounts(set_tmm),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (voom)", outline=FALSE,ylab="voom Expression", ylim=c(0,200))
#legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_tmm,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (TMM)", outline=FALSE, ylab="log TMM Expression",ylim=c(0,10))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
hist(log(counts(set_tmm)), main="Histogram of log Norm TMM Counts", xlab="log TMM Expression")
plotRLE(set_tmm,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Norm Counts TMM")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_tmm,cex=0.6,col=colors[Group], main="PCA TMM", ylim=c(-0.35, 0.4), xlim=c(-0.3, 0.3), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plot(hclust(dist(tmm_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set_tmm,log=TRUE,ylim=c(0,20), main="TMM Over-Dispersion")




#PLOT RAW COUNTS
par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Counts")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
hist(log(counts(set)), main="Histogram of log Counts", xlab="log Counts")
plotRLE(set,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Raw Counts")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set,cex=0.6,col=colors[Group], main="PCA Raw Counts", ylim=c(-0.3, 0.5), xlim=c(-0.25, 0.25), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plot(hclust(dist(t(filtered), method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set,log=TRUE,ylim=c(0,20), main="Raw Counts Over-Dispersion")

#BOXPLOTS
par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Counts")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_uq,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (UQ)", outline=FALSE, ylab="log UQ Expression",ylim=c(0,9))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_ercc,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (RUVg ERCC)", outline=FALSE, ylab="log RUVg ERCC Expression",ylim=c(0,9))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_voom,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (voom)", outline=FALSE, ylab="log voom Expression",ylim=c(-2,4))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
boxplot(set_tmm,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (TMM)", outline=FALSE, ylab="log TMM Expression",ylim=c(0,10))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)

#HISTOGRAMS
hist(log(counts(set)), main="Histogram of log Counts", xlab="log Raw Counts")
hist(log(counts(set_uq)), main="Histogram of log Norm UQ Counts", xlab="log UQ Expression")
hist(log(counts(set_ercc)), main="Histogram of log Norm RUVg ERCC Counts", xlab="log RUVg ERCC Expression")
hist(log(counts(set_voom)), main="Histogram of log Norm voom Counts", xlab="log voom Expression")
hist(log(counts(set_tmm)), main="Histogram of log Norm TMM Counts", xlab="log TMM Expression")

#RLE
plotRLE(set,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylab="RLE", main="Raw Counts")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotRLE(set_uq,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylab="RLE", main="UQ")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotRLE(set_ercc,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylab="RLE", main="RUVg ERCC")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotRLE(set_voom,ylim=c(-1,1), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylab="RLE", main="voom")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotRLE(set_tmm,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylab="RLE", main="TMM")
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)

#PCA
plotPCA(set,cex=0.6,col=colors[Group], main="PCA Raw Counts", ylim=c(-0.3, 0.5), xlim=c(-0.25, 0.25), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_uq,cex=0.6,col=colors[Group], main="PCA UQ", ylim=c(-0.3, 0.5), xlim=c(-0.25, 0.25), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_ercc,cex=0.6,col=colors[Group], main="PCA RUVg ERCC", ylim=c(-0.35, 0.55), xlim=c(-0.3, 0.3), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotMDS(voom_matrix,cex=0.8,col=colors[Group],main="PCA voom",ylim=c(-3, 3), xlim=c(-4, 4))
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)
plotPCA(set_tmm,cex=0.6,col=colors[Group], main="PCA TMM", ylim=c(-0.35, 0.4), xlim=c(-0.3, 0.3), las=2)
legend("topright", fill = colors, cex=0.6, legend = c("ASD BA9","ASD CBL", "HC BA9","HC CBL"), bty='n', horiz=T)

#HClust
plot(hclust(dist(t(filtered), method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7,main="Dendrogram Raw Counts")
plot(hclust(dist(upper_quartile_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7,main="Dendrogram UQ")
plot(hclust(dist(ercc_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7,main="Dendrogram RUVg ERCC")
plot(hclust(dist(voom_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7,main="Dendrogram voom")
plot(hclust(dist(tmm_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7,main="Dendrogram TMM")

#MeanVariance
meanVarPlot(set,log=TRUE,ylim=c(0,20), main="Raw Counts Over-Dispersion")
meanVarPlot(set_uq,log=TRUE,ylim=c(0,20), main="UQ Over-Dispersion")
meanVarPlot(set_ercc,log=TRUE,ylim=c(0,20), main="RUVg Over-Dispersion")
meanVarPlot(set_voom,log=TRUE,ylim=c(0,20), main="VOOM Over-Dispersion")
meanVarPlot(set_tmm,log=TRUE,ylim=c(0,20), main="TMM Over-Dispersion")