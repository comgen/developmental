
setwd("/Users/ibrahiaa/rnaseq_brain_data_sets/Maryland")

library(WGCNA)
library(DESeq2)
library(magrittr)
library(ggsci)
library(ggpubr)
library(plotrix)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(rutils)
library(dplyr)
library(rrvgo)
library(data.table)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(circlize)




# WGCNA CC ---------------------------------------------------------------------

# First run the PROCESS DATA module in ml_diagnostic.R

# Prepare data sets with outlier identification and removal --------------------

# normalize data using DESeq2
dds <- DESeqDataSetFromMatrix(countData=data.filt, colData=colData, design=~Region)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

# exclude outlier CC_1259 sample
cc <- subset(colData, Region=="CC")$SampleID
df <- assay(vsd)[,cc[-4]] 
# retain only protein coding genes
pc <- subset(genemap.filt, gene_biotype == "protein_coding") %>% rownames
df <- df[pc,]

# Prepare data sets
datExpr <- data.frame(df) %>% t
ctf <- read.delim("~/rnaseq_brain_data_sets/Maryland/results/ml_cell_fractions.txt",header=T)
rownames(ctf) <- ctf$SampleID
datTraits <- ctf[cc[-4],c(2:9,11,16,19,20)]
all(rownames(datTraits) == rownames(datExpr))


# Construct gene networks ------------------------------------------------------

enableWGCNAThreads() # Allow parallel execution
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function, note that we use "signed" networkType, which means we need higher threshold values
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType="signed")
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Construct gene network and identify modules using a power of 18 (large because selected signed network type)

cor <- WGCNA::cor # temporarily reassign the cor function to avoid mixing it across different packages

net.cc = blockwiseModules(datExpr, power = 18, maxBlockSize=15000,
                       TOMType = "signed", networkType="signed",minModuleSize = 40,
                       reassignThreshold = 0, mergeCutHeight = 0.4,
                       numericLabels = T, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ml_cc_TOM", 
                       verbose = 5)
cor <- stats::cor



# Save module assignment and module eigengene information for subsequent analyses
moduleLabels = net.cc$colors
me <- c("ME1","ME2","ME3","ME4","ME5","ME6","ME7","ME8","ME9","ME10","ME11","ME12","ME13")
moduleNames = labels2colors(net.cc$colors,colorSeq=me)
MEs0 = net.cc$MEs
geneTree = net.cc$dendrograms[[1]]
save(datExpr,datTraits,MEs0, moduleLabels, moduleNames, geneTree, net.cc,
     file = "ml_cc-networkConstruction-auto.RData")


# Relate modules to external variables -----------------------------------------

# Load network data saved in the previous part.
lnames = load(file = "ml_cc-networkConstruction-auto.RData")
l <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13")
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Gene relationship to trait (Age) and important modules -----------------

# Define variable 'Age' containing the Age column of datTraits
age = as.data.frame(datTraits$Age);
names(age) = "Age"
# names (colors) of the modules
modNames = l
MEs = orderMEs(MEs0[,-14]) # exclude module 0

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, age, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(age), sep="")
names(GSPvalue) = paste("p.GS.", names(age), sep="")



# Summary output of network analysis results -----------------------------------

# Create the starting data frame
a <- moduleNames
b <- c(names(MEs),"grey")
c <- c(l,"grey")
geneInfo0 = data.frame(ensemblGeneID = colnames(datExpr),
                       geneSymbol = genemap[colnames(datExpr),"external_gene_name"],
                       moduleName = moduleNames,
                       moduleID = a[a %in% b] <- c[match(a, b)],
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for Age
modOrder = order(-abs(cor(MEs, age, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module ID, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleID, -abs(geneInfo0$GS.Age))
geneInfo = geneInfo0[geneOrder, ]
geneInfo$moduleID <- factor(geneInfo$moduleID, levels=c(l,"grey"))
write.table(geneInfo,"~/rnaseq_brain_data_sets/Maryland/results/wgcna_geneInfo_cc.txt",sep="\t",row.names=F)


# Time series plots of all modules ---------------------------------------------

df <- MEs
colnames(df) <- l
df$SampleID <- rownames(df)
df <- merge(melt(df), datTraits[,c("Age","SampleID")],by="SampleID")

col <- get_palette("aaas", 13)


glist <- lapply(1:13, function(x){
  modSize <- subset(geneInfo, moduleID == levels(df$variable)[x]) %>% nrow
  df.sub <- subset(df, variable == levels(df$variable)[x])
  ggplot(data=df.sub, aes(x=Age, y=value, color=variable)) + 
    geom_smooth(se=T, method="loess", level=0.5,span= 0.76) +
    geom_jitter(width=0.2, height=0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(color="black",size=15),
          axis.text.y=element_text(color="black",size=15),
          plot.title=element_text(size=15,hjust=0.5),
          legend.position="",
          panel.background=element_rect(linetype="solid", color="black", size=1)) +
    scale_colour_manual(values=col[x]) +
    labs(title=paste0(levels(df$variable)[x]," (", modSize,")"), y="", x = "") +
    geom_vline(xintercept = c(5), linetype = "dotted", color = "black", size=0.5) +
    geom_vline(xintercept = c(10), linetype = "dotted", color = "black", size=0.5) +
    geom_vline(xintercept = c(15), linetype = "dotted", color = "black", size=0.5) +
    geom_vline(xintercept = c(20), linetype = "dotted", color = "black", size=0.5) +
    annotate("rect", xmin = -Inf, xmax = 5, ymin = -Inf, ymax = Inf, alpha = .1)
})

annotate_figure(ggarrange(plotlist=glist,ncol=5,nrow=3),
                left=text_grob("Module eigengene",size=16,rot=90),
                bottom=text_grob("Age (Y)",size=16,rot=0),
                top=text_grob("All modules, CC\n",size=16))


# Time series analysis to test significant deviation over time -----------------

mes <- MEs %>% t
rownames(mes) <- l
pheno <- datTraits
pheno$age_bins <- ifelse(pheno$Age>0 & pheno$Age<6, "age_0_5", 
                         ifelse(pheno$Age>5 & pheno$Age<11, "age_5_10",
                                ifelse(pheno$Age>10 & pheno$Age<16, "age_10_15",
                                       ifelse(pheno$Age>15 & pheno$Age<21, "age_15_20", "age_20_25"))))
all( colnames(mes) == rownames(pheno))

lev <- c("age_0_5", "age_5_10", "age_10_15", "age_15_20", "age_20_25")
f <- factor(pheno$age_bins, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(mes, design)
cm <- makeContrasts(
  A = age_5_10 - age_0_5,
  B = age_10_15 - age_0_5,
  C = age_15_20 - age_0_5,
  D = age_20_25 - age_0_5,
  levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2, trend=TRUE)
summary(decideTests(fit2, p.value=0.1))

resA <- topTable(fit2, coef="A", n=Inf, p.value=1) # no significant hits
resA$Period <- "5-10Y"
resA$Module <- rownames(resA)
resB <- topTable(fit2, coef="B", n=Inf, p.value=1)
resB$Period <- "10-15Y"
resB$Module <- rownames(resB)
resC <- topTable(fit2, coef="C", n=Inf, p.value=1)
resC$Period <- "15-20Y"
resC$Module <- rownames(resC)
resD <- topTable(fit2, coef="D", n=Inf, p.value=1)
resD$Period <- "20-25Y"
resD$Module <- rownames(resD)
res <- rbind(resA,resB,resC,resD)


## Heatmap of time series analysis results

df <- res
df$Module <- factor(df$Module, levels=l)
df$logFCsig <- ifelse(df$adj.P.Val<0.1, df$logFC, 0)
df$pAdjSig <- ifelse(df$adj.P.Val<0.1, "*", "")


fc <- reshape2::dcast(df[,c(1,7,8)], Module ~ Period, value.var="logFC")
rownames(fc) <- fc$Module
fc <- fc[,c(5,2:4)]

p <- reshape2::dcast(df[,c(7,8,10)], Module ~ Period, value.var="pAdjSig")
rownames(p) <- p$Module
p <- p[,c(5,2:4)]

col <- get_palette("aaas", 13)
colfunc <- colorRampPalette(c(col[6],"white", col[8]))
par(mar=c(7,4,4,2))
labeledHeatmap(Matrix = fc,
               xLabels = c("5-10Y", "10-15Y","15-20Y","20-25Y"),
               yLabels = rownames(fc),
               ySymbols = rownames(fc),
               colorLabels = F,
               colors = colfunc(50),
               textMatrix = p,
               setStdMargins = FALSE,
               textAdj = c(0.5, 0.75),
               cex.text = 1.9,
               cex.lab.x = 1.3,
               cex.lab.y = 1.3,
               cex.main=1.3,
               legendLabel = "logFC(eigengene)",
               cex.legendLabel = 1.3,
               xLabelsAngle=45,
               xLabelsAdj=1,
               main = "Change over time\nrelative to baseline (0-5Y)",
               font.main=1,
               adj=0)


# Time series plots of significant modules -------------------------------------

df <- MEs
colnames(df) <- l
df$SampleID <- rownames(df)
df <- merge(melt(df), datTraits[,c("Age","SampleID")],by="SampleID")

col <- get_palette("aaas", 13)

glist <- lapply(c(1:3,9,11,13), function(x){
  modSize <- subset(geneInfo, moduleID == levels(df$variable)[x]) %>% nrow
  df.sub <- subset(df, variable == levels(df$variable)[x])
  ggplot(data=df.sub, aes(x=Age, y=value, color=variable)) + 
    geom_smooth(se=T, method="loess", level=0.5,span= 0.76) +
    geom_jitter(width=0.2, height=0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(color="black",size=15),
          axis.text.y=element_text(color="black",size=15),
          plot.title=element_text(size=15,hjust=0.5),
          legend.position="",
          plot.margin=unit(c(0,0.2,0,0), "cm"),
          panel.background=element_rect(linetype="solid", color="black", size=1)) +
    scale_colour_manual(values=col[x]) +
    labs(title=paste0(levels(df$variable)[x]," (", modSize,")"), y="", x = "") 
})

annotate_figure(ggarrange(plotlist=glist,ncol=6,nrow=1),
                left=text_grob("Module eigengene",size=16,rot=90),
                bottom=text_grob("Age (Y)",size=15,rot=0),
                top=text_grob("Dynamic modules, CC\n",size=16,rot=0))




# GO analysis of top modules ---------------------------------------------------

egolist <- lapply(c(1:3,9,11,13), function(x){
  m <- levels(factor(geneInfo$moduleID))[x]
  dg <- subset(geneInfo,moduleID==m)
  dg <- dg[,c("geneSymbol","moduleID","GS.Age","p.GS.Age",paste0("MM.",m))]
  colnames(dg)[5] <- "MM"
  dg <- subset(dg, MM>0.7)  # only include genes with MM>0.7
  dg <- dg[order(dg$p.GS.Age,decreasing=F),]

  ego <- enrichGO(gene=dg$geneSymbol, ont="BP", OrgDb="org.Hs.eg.db", pvalueCutoff=0.1, qvalueCutoff=0.1, pAdjustMethod="fdr", keyType="SYMBOL") %>% simplify(cutoff=0.5) %>% data.frame 
  if( nrow(ego) > 1){
  simMatrix <- calculateSimMatrix(ego$ID, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
  scores <- setNames(-log10(ego$p.adjust), ego$ID)
  red <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Hs.eg.db")
  ego <- merge(ego, red, by=0)
  }
  ego$moduleID <- m
  assign(paste0("ego",m),ego)
})

egoFull <- rbindlist(egolist, fill=TRUE) %>% data.frame
egoFull <- egoFull[ with(egoFull, order(moduleID,pvalue)),]
write.table(egoFull,"~/rnaseq_brain_data_sets/Maryland/results/wgcna_egoFull_cc.txt",sep="\t",row.names=F)

egoFullRed <- distinct(egoFull, moduleID,parentTerm, .keep_all= TRUE) %>% data.frame
egoTop <- Reduce(rbind, 
                 by(egoFullRed, egoFullRed["moduleID"], head, n = 6)
                 )
write.table(egoTop[,c(6,7,10,13,17,21)],"~/rnaseq_brain_data_sets/Maryland/results/wgcna_egoTop_cc.txt",sep="\t",row.names=F) # shorten GO descriptions manually

# Plot GO results
egoTopShort <- read.delim("~/rnaseq_brain_data_sets/Maryland/results/wgcna_egoTop_short_cc.txt",sep="\t",header=T)
egoTopShort <- egoTopShort[ with(egoTopShort, order(moduleID,pvalue)),]
egoTopShort <- rbind( #need correct order for plotting
  subset(egoTopShort,moduleID=="M1"),
  subset(egoTopShort,moduleID=="M2"),
  subset(egoTopShort,moduleID=="M3"),
  subset(egoTopShort,moduleID=="M9"),
  subset(egoTopShort,moduleID=="M11"),
  subset(egoTopShort,moduleID=="M13")
)

egoTopShort$moduleID <- factor(egoTopShort$moduleID, levels= c("M13","M11","M9","M3","M2","M1"))

order <- egoTopShort$short_name
egoTopShort$short_name <- factor(egoTopShort$short_name, levels=order)
col <- get_palette("aaas", 13)[c(13,11,9,3,2,1)]

ggplot(egoTopShort, aes(x=short_name, y=moduleID, size = -log10(pvalue), fill=moduleID)) + 
  geom_point(alpha=0.8, shape=21, col="black") +
  scale_size(range = c(5, 11)) +
  scale_fill_manual(values=col) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=15,color="black", angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=15,color="black"),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  guides(fill = "none") +
  labs(x="",y="",title="")


