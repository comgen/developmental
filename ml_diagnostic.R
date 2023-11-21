
setwd("/Users/ibrahiaa/rnaseq_brain_data_sets/Maryland")

library(purrr)
library(tidyverse)
library(biomaRt)
library(edgeR)
library(DESeq2)
library(plyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gplots)
library(reshape2)
library(magrittr)
library(ggsci)
library(gridExtra)
library(ggfortify)
library(data.table)




# PROCESS DATA -----------------------------------------------------------------

# Read in data -----------------------------------------------------------------

f_files <- list.files("/Users/ibrahiaa/rnaseq_brain_data_sets/Maryland/count_data", pattern = ".counts", full.names = T)

read_in_feature_counts <- function(file){
  cnt <- read_tsv(file, col_names =T, comment = "#")
  cnt <- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  return(cnt)
}

count <- map(f_files, read_in_feature_counts)
count <- purrr::reduce(count, inner_join)
data <- as.data.frame(count)
rownames(data) <- data$Geneid
data <- data[,-1]
colnames(data) <- gsub("-BA11_sorted.bam", "_BA11", colnames(data))
colnames(data) <- gsub("sorted.bam", "CC", colnames(data))

colData <- read.csv("/Users/ibrahiaa/rnaseq_brain_data_sets/Maryland/maryland_colData.txt", sep="\t", header=TRUE)
rownames(colData) <- colData$SampleID
data <- data[, colData$SeqID]
colnames(data) <- colData$SampleID
all(colnames(data) == rownames(colData))

colData$Region <- factor(colData$Region)
colData$Period <- ifelse(colData$Age<12,"Childhood","Adolescence") # periods as defined by Kang et al. 2011

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "gene_biotype", "chromosome_name"), filters="ensembl_gene_id", values=rownames(data), mart=ensembl)
genemap$external_gene_name_2 <- ifelse(genemap$external_gene_name == "",genemap$ensembl_gene_id, genemap$external_gene_name)
genemap <- genemap[duplicated(genemap$external_gene_name_2) == F,]
rownames(genemap) <- genemap$ensembl_gene_id


# Pre-filter genes -------------------------------------------------------------

keep <- filterByExpr(data, design=colData$Region)
data.filt <- data[keep,]

genemap.filt <- genemap[rownames(data.filt),]
genemap.filt <- subset(genemap.filt,gene_biotype %in% c("protein_coding","lncRNA"))
data.filt <- data.filt[rownames(genemap.filt),] ## keep only proten coding and lncRNA 




# UNSUPERVISED CLUSTERING ------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=data.filt, colData=colData, design=~Region)
dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup=c("Region","Period","DonorID","Ethnicity","PMI", "Age", "AgeDays"), ntop=500, returnData=TRUE)
pca$AgeExact <- pca$Age + (pca$AgeDays)/365
pca$RegionOut <- as.character(pca$Region)
pca["CC_1259", "RegionOut"] <- "CC_out"

col <- get_palette("aaas", 10)

# normal PCA plot
ggplot(pca, aes(x = -PC1, y = -PC2, col=Region,fill=Period)) +
  scale_fill_manual(values=c(col[10],"grey")) +
  scale_color_manual(values=col[c(4,5)]) +
  geom_line(aes(group = DonorID), col="grey") +
  geom_point(size=6.5, alpha=0.7, shape=21, stroke=1) +
  guides(color = "none") +
  xlab("PC1: 73% variance") + 
  ylab("PC2: 6% variance") +
  theme_classic() +
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        legend.title=element_blank(),
        legend.direction="vertical",
        legend.position=c(0.001,0.999),
        legend.justification=c(0.001,0.999),
        legend.background=element_rect(linetype="solid", color="black"),
        panel.background=element_rect(linetype="solid", color="black", linewidth=1))

# PCA with Age
ggplot(pca, aes(x = AgeExact, y = PC1,fill=Region, col=RegionOut)) +
  scale_color_manual(values=c("black",col[2],"black")) +
  scale_fill_discrete(labels = c("CC, n=38", "DFC, n=41"), type=col[4:5]) +
  guides(color="none") +
  geom_line(aes(group = DonorID), col="grey") +
  geom_point(size=6.5, alpha=0.7, shape=21, stroke=1) +
  labs(x="Age",y="PC1") +
  lims(y=c(-80,110)) +
  theme_classic() +
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=17),
        legend.direction="horizontal",
        panel.background=element_rect(linetype="solid", color="black", linewidth=1)) +
  scale_x_continuous(limits = c(1, 24.2), breaks = c(seq(1,21,2),24)) +
  geom_vline(xintercept = c(12), linetype = "dotted", color = "black", size=0.7) +
  annotate("text", x=6, y=105, label="Childhood", size=6, hjust=0.5) +
  annotate("text", x=18, y=105, label="Adolescence", size=6, hjust=0.5) 





# COMPARE COVERAGE IN Maryland vs. BrainSpan -----------------------------------

ml <- subset(colData, Region == "DFC")[,c("DonorID","Age","Sex")]
ml$Sex <- "M"
ml$Stage <- ifelse(ml$Age<6,"1-5 Y",
                   ifelse(ml$Age>5 & ml$Age<12, "6-11 Y",
                          ifelse(ml$Age>11 & ml$Age<20, "12-19 Y","20-40 Y")))
ml$Study <- "Maryland"


bs <- read.csv2("/Users/ibrahiaa/organoids_scz/developmental_correlation/metaCol_ageStage.txt", header=T, sep="\t")
bs <- subset(bs, structure_acronym=="DFC")
bs <- bs[,c("donor_id","age_alt","gender","Stage")]
bs$Study <- "BrainSpan"
colnames(bs) <- c("DonorID","Age","Sex","Stage","Study")

df <- rbind(ml,bs)
#df$Stage <- factor(df$Stage, levels=c("8-9 PCW","10-12 PCW","13-15 PCW","16-18 PCW","19-24 PCW","25-38 PCW","Birth - 5 M", "6-12 M","1-5 Y","6-11 Y","12-19 Y", "20-40 Y"))
df$Age <- factor(df$Age, levels=c("8-9 PCW","10-12 PCW","13-15 PCW","16-18 PCW","19-24 PCW","25-38 PCW","0","1","2","3","4","5","6","8","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","30","36","37","40"))

col <- get_palette("aaas", 10)
ggplot(df, aes(x = Age, y = Study, fill=Study)) +
  scale_fill_manual(values=col[c(1,6)]) +
  geom_jitter(size=5.5, alpha=0.85, height=0.2, width=0.1, shape=21) +
  labs(x="",y="") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=17,color="black", angle=90, vjust=0.5, hjust=1), 
        axis.text.y=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        legend.position="n",
        panel.background=element_rect(linetype="solid", color="black", linewidth=1)) +
  geom_vline(xintercept = 7.5, linetype="dotted", color = "black", size=0.5) +
  geom_vline(xintercept = 29.5, linetype="dotted", color = "black", size=0.5)




## DECONVOLTION WITH InstaPrism ------------------------------------------------

library(InstaPrism)

# single-cell reference data source, alt 1: https://www.nature.com/articles/nbt.4038#Sec1, table S2
# single-cell reference data source, alt 2: https://www.sciencedirect.com/science/article/pii/S0092867422012582?via%3Dihub#da0010, http://brain.listerlab.org

# bulk mixture file
bk.dat <- data.filt  %>% data.frame
all(rownames(genemap.filt) == rownames(bk.dat))
rownames(bk.dat) <- genemap.filt$external_gene_name_2

# single-cell ref data, alt 1
#library(data.table)
#ref <- fread("/Users/ibrahiaa/completed_projects/cytokines_ueland/cognition/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt", sep="\t", header=F) 
#ref <- data.frame(ref)
#anno <- read.csv("/Users/ibrahiaa/completed_projects/cytokines_ueland/cognition/PMID29227469_snRef_frontalCortex_tableS2.txt", sep="\t",check.names=T) 
#colnames(ref) <- c("Gene",anno$Sample.Names..Library_Barcode.)
#rownames(ref) <- ref$Gene
#sc.dat <- ref[,-1] %>% t %>% data.frame

# single-cell ref data, alt 2
library(zellkonverter)
library(SingleCellExperiment)
library(BayesPrism)
library(InstaPrism)

dat <- readH5AD("/Users/ibrahiaa/rnaseq_brain_data_sets/rna_single_cell_refs/DFC_developmental/Processed_data-RNA-all_full-counts-and-downsampled-CPM.h5ad")
ref <- assays(dat)$X 

anno <- colData(dat) %>% data.frame
annoSubStage <- subset(anno, stage_id %in% c("Childhood","Adolescence","Adult"))
annoSubClust <- subset(annoSubStage, !(major_clust %in% c("PN_dev","CGE_dev","MGE_dev","Poor-Quality"))) # remove low-count cell types
annoSubClust$major_clust <- droplevels(annoSubClust$major_clust)
sub.sorted <- table(annoSubClust$sub_clust) %>% sort(decreasing=T) %>% data.frame
annoSubClust <- subset(annoSubClust, sub_clust %in% head(sub.sorted,67)$Var1) # remove sub clust with <10 cells
annoSubClust$sub_clust <- droplevels(annoSubClust$sub_clust)


# aggregate across cell states as BayesPrism cannot take sparse matrix as input. This procedure is recommended by BayesPrism authors
df <- ref[,rownames(annoSubClust)] %>% t
mlist <- lapply(levels(annoSubClust$sub_clust), function(x){
  s <- df[subset(annoSubClust, sub_clust==x) %>% rownames,]
  s.sum <- s %>% colSums %>% data.frame 
  colnames(s.sum) <- x
  t(s.sum) %>% data.frame
})
ref.coll <- rbindlist(mlist) %>% data.frame
rownames(ref.coll) <- levels(annoSubClust$sub_clust)
sc.dat <- ref.coll


annoSubClust$super_clust <- mapvalues(annoSubClust$major_clust, from = c("L2-3_CUX2","L4_RORB", "L5-6_THEMIS","L5-6_TLE4"), to = rep("Exc",4))
annoSubClust$super_clust <- mapvalues(annoSubClust$super_clust, from = c("VIP","ID2","LAMP5_NOS1","SST","PV","PV_SCUBE3"), to = rep("Inh",6))
t <- annoSubClust[ duplicated(annoSubClust$sub_clust) ==F, ] 

cell.state.labels <- rownames(ref.coll)
cell.type.labels <- t[ order(t$sub_clust), "super_clust"] %>% as.character

bk.dat <- data.filt %>% t %>% data.frame
all(rownames(genemap.filt) == colnames(bk.dat))
colnames(bk.dat) <- genemap.filt$external_gene_name_2

rm(dat,ref) # very large files

# Filter uninformative genes                     
sc.dat.filtered <- cleanup.genes(
  input=sc.dat,
  input.type="GEP",
  species="hs", 
  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1") ,
  exp.cells=5
  )

# Filter on protein coding genes as they show best correlation
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered, gene.type = "protein_coding")

# Construct Prism object
myPrism <- new.prism(
  reference= sc.dat.filtered.pc, # use protein coding sc.dat
  mixture=bk.dat,
  input.type="GEP", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1)

# Run deconvolution using InstaPrism to save time
instaPrismRes = InstaPrism_legacy(
  input_type = "prism",
  prismObj = myPrism,
  return.Z.cs = F,
  return.Z.ct  = T,
  update = T)
slotNames(instaPrismRes)


# Plot estimated cell type fractions
frac <- instaPrismRes@Post.updated.ct@theta %>% t
df <- merge(frac,colData, by=0)
write.table(df,"~/rnaseq_brain_data_sets/Maryland/results/ml_cell_fractions.txt", sep="\t", row.names=F, quote=F)
df <- reshape2::melt(df[,c(2:8,12)])
df$variable <- factor(df$variable, levels=c("Astro","Inh", "Exc", "OPC", "Oligo", "Micro","Vas"))

col <- get_palette("aaas", 10)
c <- c("CC"=col[4],"DFC"=col[5])

for( i in c("CC","DFC")){
g <- ggplot(subset(df,Region==i), aes(x=variable, y=value,col=Region)) + 
  geom_point(alpha=0.35, position=position_jitterdodge(0.2)) +
  geom_boxplot(alpha=0) + 
  scale_color_manual(values=c[i]) +
  theme_classic() +
  theme(axis.text.x=element_text(color="black",size=17, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(color="black",size=17),
        axis.title=element_text(size=17),
        legend.position="n",
        legend.background=element_rect(linetype="solid", color="black"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  labs(title="", y="Estimated fractions", x = "") +
  annotate("text", x=0.7, y=1, size=6, label = i, hjust=0, vjust=0.7) 
assign(paste0("g",i), g)
}  
ggarrange(gCC,gDFC,ncol=2,nrow=1)


# Plot estimated fractions highlighting outlier sample
frac <- instaPrismRes@Post.updated.ct@theta %>% t
df <- merge(frac,colData, by=0)
df$outlier <- ifelse(df$SampleID =="CC_1259", "Outlier", ifelse(df$Region=="CC", "CC", "DFC"))
df <- reshape2::melt(df[,c(2:8,10,12,21)])
df$variable <- factor(df$variable, levels=c("Astro","Inh", "Exc", "OPC", "Oligo", "Micro","Vas"))

col <- get_palette("aaas", 10)
ggplot(df, aes(fill=variable, y=value, x=outlier)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=col) +
  theme_classic() +
  theme(axis.text.x=element_text(color="black",size=15, angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=15),
        panel.background=element_rect(linetype="solid", color="black", size=1),
        legend.title=element_blank(),
        legend.text=element_text(size=15)) +
  labs(x="", y="Estimated fractions (avg.)")


# Plot estimated fractions by region across age
frac <- instaPrismRes@Post.updated.ct@theta %>% t
df <- reshape2::melt(frac)
df <- merge(df,colData[,c(2,4,7)], by.x="Var1", by.y="SampleID")
df$variable <- factor(df$variable, levels=c("Astro","Inh", "Exc", "OPC", "Oligo", "Micro","Vas"))

col <- get_palette("aaas", 10)
for( i in c("CC","DFC")){
g <- ggplot(subset(df,Region==i), aes(x=Age, y=value, col=Var2)) + 
  geom_smooth(se=T, method="loess", level=0.3,span= 0.9) +
  geom_point(alpha=0.35, position=position_jitter(widt=0.5)) +
  theme_classic() +
  theme(axis.text=element_text(color="black",size=17),
        axis.title=element_text(size=17),
        plot.title=element_text(size=17, hjust=0.5),
        legend.position="n", 
        legend.title=element_blank(),
        legend.text=element_text(size=17),
        legend.spacing=unit(0.5,"cm"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  scale_colour_manual(values=col) +
  scale_fill_manual(values=col) +
  labs(title="", y="Estimated fractions", x = "Age") +
  annotate("text", x=0, y=1, size=6, label = i, hjust=0, vjust=0.7) 
assign(paste0("g",i), g)
}  
ggarrange(gCC,gDFC,ncol=2,nrow=1,common.legend=T, legend="right")



# Adjust Oligo-specific gene expression ----------------------------------------

frac <- instaPrismRes@Post.updated.ct@theta %>% t %>% data.frame
expr <- instaPrismRes@Post.ini.ct@Z

# save cell type-specific expression values
write.table(expr[,,"Astro"],"~/rnaseq_brain_data_sets/Maryland/results/ml_astro_expression.txt", sep="\t", row.names=F, quote=F)
write.table(expr[,,"Oligo"],"~/rnaseq_brain_data_sets/Maryland/results/ml_oligo_expression.txt", sep="\t", row.names=F, quote=F)

y <- expr[,,"Oligo"]
y2 <- y * (1 - frac$Oligo) # adjust for oligo fractions

df <- edgeR::cpm(y) # normalize for library size
df <- df[ rev(rownames(y)),] # reverse row order for plotting
df <- df[,sample(ncol(df),1000)] # pick 1000 random genes 
df <- reshape2::melt(df)
df <- merge(df, colData[,c("SampleID","Region")], by.x="Var1", by.y="SampleID")

col <- get_palette("aaas", 10) 
g1 <- ggplot(df, aes(x=Var1, y=value, col=Region)) + 
  geom_boxplot(outlier.size=1) +
  scale_colour_manual(values=col[4:5]) +
  lims(y=c(0,2e5)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color="black"),
        axis.text.y=element_text(color="black"),
        axis.title=element_text(size=15),
        plot.title=element_text(size=15, hjust=0.5),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  labs(title="No adjustment",x="",y="") 

df <- edgeR::cpm(y2) # normalize for library size
df <- df[ rev(rownames(y)),] # reverse row order for plotting
df <- df[,sample(ncol(df),1000)] # pick 1000 random genes 
df <- reshape2::melt(df)
df <- merge(df, colData[,c("SampleID","Region")], by.x="Var1", by.y="SampleID")

col <- get_palette("aaas", 10) 
g2 <- ggplot(df, aes(x=Var1, y=value, col=Region)) + 
  geom_boxplot(outlier.size=1) +
  scale_colour_manual(values=col[4:5]) +
  lims(y=c(0,2e5)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color="black"),
        axis.text.y=element_text(color="black"),
        axis.title=element_text(size=15),
        plot.title=element_text(size=15, hjust=0.5),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  labs(title="Adjusted for Oligo fractions",x="",y="") 

annotate_figure(ggarrange(g1,g2,ncol=1,nrow=2,common.legend=T, legend="bottom"),
                left=text_grob("Oligo-specific expression",size=15,rot=90))




# Adjust Astro-specific gene expression ----------------------------------------

frac <- instaPrismRes@Post.updated.ct@theta %>% t %>% data.frame
expr <- instaPrismRes@Post.ini.ct@Z
y <- expr[,,"Astro"]
y2 <- y * (1 - frac$Astro) # adjust for astro fractions

df <- edgeR::cpm(y) # normalize for library size
df <- df[ rev(rownames(y)),] # reverse row order for plotting
df <- df[,sample(ncol(df),1000)] # pick 1000 random genes 
df <- reshape2::melt(df)
df <- merge(df, colData[,c("SampleID","Region")], by.x="Var1", by.y="SampleID")

col <- get_palette("aaas", 10) 
g1 <- ggplot(df, aes(x=Var1, y=value, col=Region)) + 
  geom_boxplot(outlier.size=1) +
  scale_colour_manual(values=col[4:5]) +
  lims(y=c(0,2e5)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color="black"),
        axis.text.y=element_text(color="black"),
        axis.title=element_text(size=15),
        plot.title=element_text(size=15, hjust=0.5),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  labs(title="No adjustment",x="",y="") 

df <- edgeR::cpm(y2) # normalize for library size
df <- df[ rev(rownames(y)),] # reverse row order for plotting
df <- df[,sample(ncol(df),1000)] # pick 1000 random genes 
df <- reshape2::melt(df)
df <- merge(df, colData[,c("SampleID","Region")], by.x="Var1", by.y="SampleID")

col <- get_palette("aaas", 10) 
g2 <- ggplot(df, aes(x=Var1, y=value, col=Region)) + 
  geom_boxplot(outlier.size=1) +
  scale_colour_manual(values=col[4:5]) +
  lims(y=c(0,2e5)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, color="black"),
        axis.text.y=element_text(color="black"),
        axis.title=element_text(size=15),
        plot.title=element_text(size=15, hjust=0.5),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        plot.margin = unit(c(0.1,0.1,0,0), "cm"),
        panel.background=element_rect(linetype="solid", color="black", size=1)) +
  labs(title="Adjusted for Astro fractions",x="",y="") 

annotate_figure(ggarrange(g1,g2,ncol=1,nrow=2,common.legend=T, legend="bottom"),
                left=text_grob("Astro-specific expression",size=15,rot=90))


# Celltype-specific gene expression --------------------------------------------


library(ggfortify)

frac <- instaPrismRes@Post.updated.ct@theta %>% t %>% data.frame
expr <- instaPrismRes@Post.ini.ct@Z

yAstro <- expr[,,"Astro"]
yAstro <- yAstro * (1 - frac$Astro) # adjust for Astro fractions
rownames(yAstro) <- paste0(rownames(yAstro), "_astro")

yOligo <- expr[,,"Oligo"]
yOligo <- yOligo * (1 - frac$Oligo) # adjust for Oligo fractions
rownames(yOligo) <- paste0(rownames(yOligo), "_oligo")

y <- rbind(yAstro, yOligo)
labels <- c( rep("DFC_Astro",41), rep("CC_Astro",38), rep("DFC_Oligo",41), rep("CC_Oligo",38) )
Region <- c( rep("DFC",41), rep("CC",38), rep("DFC",41), rep("CC",38) )
Celltype <- c( rep("Astro",41), rep("Astro",38), rep("Oligo",41), rep("Oligo",38) )
donorid <- gsub("BA11_", "ML",
                gsub("CC_", "ML",
                     gsub("_astro", "",
                          gsub("_oligo", "", rownames(y)))))
pheno <- data.frame(row.names=rownames(y), SampleID=rownames(y), Age=rep(colData$Age,2),DonorID=donorid, Region_Type=labels, Fractions=c(frac$Astro, frac$Oligo), Region=Region, Celltype=Celltype)

vsd <- vst(round(t(y)))
pca_res <- prcomp(t(vsd))
df <- merge(pheno,pca_res$x[,1:2], by=0)
col <- get_palette("futurama", 10)

# PCA plot
ggplot(df, aes(x = PC1, y = PC2, col=Celltype,fill=Celltype)) +
  scale_fill_manual(values=col[5:6]) +
  scale_color_manual(values=col[5:6]) +
  geom_point(size=6.5, alpha=0.9, shape=21, stroke=1) +
  labs(x="PC1: 73.9%", y="PC2: 11.7%") +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        legend.title=element_blank(),
        legend.position=c(0.5,0.999),
        legend.justification=c(0.5,0.999),
        legend.direction="horizontal",
        legend.background=element_rect(linetype="solid", color="black", linewidth=0.5),
        panel.background=element_rect(linetype="solid", color="black", linewidth=1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.8) +
  annotate("text", x=min(df$PC1), y=max(df$PC2), label="CC", size=6.5, hjust=0) +
  annotate("text", x=min(df$PC1), y=min(df$PC2), label="DFC", size=6.5, hjust=0)


# PCA plot colored by Age
ggplot(df, aes(x = PC1, y = PC2, col=Age,fill=Age)) +
  scale_fill_gradient(high = "#132B43", low = "#56B1F7", name="Age") +
  scale_color_gradient(high = "#132B43", low = "#56B1F7") +
  geom_point(size=6.5, alpha=0.9, shape=21, stroke=1) +
  labs(x="PC1: 73.9%", y="PC2: 11.7%") +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        legend.title=element_text(size=17),
        legend.position=c(0.5,0.999),
        legend.justification=c(0.5,0.999),
        legend.direction="horizontal",
        legend.key.width= unit(0.8, 'cm'),
        legend.background=element_rect(linetype="solid", color="black", linewidth=0.5),
        panel.background=element_rect(linetype="solid", color="black", linewidth=1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.8) +
  annotate("text", x=min(df$PC1), y=max(df$PC2), label="CC", size=6.5, hjust=0) +
  annotate("text", x=min(df$PC1), y=min(df$PC2), label="DFC", size=6.5, hjust=0)



# PCA plot colored by marker expression
markers <- c("AQP4","S100B","NDRG2","ALDH1L1", "MOG","CNP","SOX10","OLIG1")
markerExpr <- vsd[markers,] %>% t %>% data.frame
dfMark <- merge(df[,-1],markerExpr,by.x="SampleID", by.y=0)

labels <- c(paste0("Astro: ", markers[1:4]), paste0("Oligo: ", markers[5:8]))

glist <- lapply(10:17, function(x){
ggplot(dfMark, aes(x = PC1, y = PC2, col=dfMark[,x], fill=dfMark[,x])) +
  scale_fill_gradient(high = "#132B43", low = "#56B1F7", name="Expression") +
  scale_color_gradient(high = "#132B43", low = "#56B1F7") +
  geom_point(size=4, alpha=0.7, shape=21, stroke=1) +
  labs(x="", y="", title= labels[x-9]) +
  guides(color = "none") +
  theme_classic() +  
  theme(axis.text=element_text(size=17,color="black"), 
        axis.title=element_text(size=17),
        plot.title=element_text(size=17, hjust=0.5),
        legend.text=element_text(size=17),
        legend.title=element_text(size=17),
        legend.position="right",
        legend.key.width= unit(1.2, 'cm'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        panel.background=element_rect(linetype="solid", color="black", linewidth=1))
})
annotate_figure(ggarrange(plotlist=glist,ncol=4,nrow=2, common.legend=T, legend="bottom"))





# DE analysis of cell type-specific expression paired design -------------------

frac <- instaPrismRes@Post.updated.ct@theta %>% t %>% data.frame
expr <- instaPrismRes@Post.ini.ct@Z

astro <- expr[,,"Astro"]
rownames(astro) <- paste0(rownames(astro), "_astro")

oligo <- expr[,,"Oligo"]
rownames(oligo) <- paste0(rownames(yOligo), "_oligo")

df <- rbind(astro, oligo) %>% t
labels <- c( rep("DFC_Astro",41), rep("CC_Astro",38), rep("DFC_Oligo",41), rep("CC_Oligo",38) )
Region <- c( rep("DFC",41), rep("CC",38), rep("DFC",41), rep("CC",38) )
Celltype <- c( rep("Astro",41), rep("Astro",38), rep("Oligo",41), rep("Oligo",38) )
donorid <- gsub("BA11_", "ML",
                gsub("CC_", "ML",
                     gsub("_astro", "",
                          gsub("_oligo", "", colnames(df)))))
pheno <- data.frame(row.names=colnames(df), SampleID=colnames(df), DonorID=donorid, Region_Type=labels, Fractions=c(frac$Astro, frac$Oligo), Region=Region, Celltype=Celltype)

genes <- subset(genemap,external_gene_name %in% rownames(df))
rownames(genes) <- genes$external_gene_name
genes <- genes[rownames(df),]

y <- DGEList(counts = df, genes = genes, samples=pheno)
y <- calcNormFactors(y) # normalize 

design <- model.matrix(~ 0 + Region_Type + Fractions, data = pheno)
colnames(design)[1:4] <- levels(factor(pheno$Region_Type))

v <- voomWithQualityWeights(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = pheno$DonorID)
cor$consensus

fit <- lmFit(v, design, block = pheno$DonorID, correlation = cor$consensus)

cm <- makeContrasts(
  Oligo_DFC_vs_CC = (CC_Oligo - DFC_Oligo),
  Astro_DFC_vs_CC = (CC_Astro - DFC_Astro),
  DFC_Oligo_vs_Astro = (DFC_Oligo - DFC_Astro),
  CC_Oligo_vs_Astro = (CC_Oligo - CC_Astro),
  levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
summary(decideTests(fit2, p.value=0.05))

resOligo <- topTable(fit2, coef="Oligo_DFC_vs_CC", n=Inf, p.value=0.05)
resOligo$Comparison <- "Oligo_DFC_vs_CC"
resAstro <- topTable(fit2, coef="Astro_DFC_vs_CC", n=Inf, p.value=0.05)
resAstro$Comparison <- "Astro_DFC_vs_CC"
res<- rbind(resOligo,resAstro)



# plot top protein coding significant genes
col <- get_palette("futurama", 10)
samp <- subset(pheno, Celltype=="Oligo")$SampleID
g <- edgeR::cpm(y[,samp], log=TRUE, prior.count=3) %>% t %>% data.frame
g <- merge(g,pheno[samp,],by=0) %>% melt

glist <- lapply(1:10, function (x) {
  gene <- (resOligo$external_gene_name)[x]
  n <- subset(g, variable == gene)
  ggplot(n,
         aes(x = Region, y = value, color = Region, group = Region)) + 
    geom_jitter() + 
    geom_boxplot() +
    labs(title=gene,x="",y="") +
    scale_colour_manual(values=col[c(1,2)]) +
    theme_classic() +
    theme(axis.text.x=element_text(size=16,color="black", angle=50, hjust=1),
          axis.text.y=element_text(size=16,color="black"),
          axis.title=element_text(size=16),
          plot.title=element_text(size=16,hjust=0.5),
          legend.text=element_text(size=16),
          legend.title=element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm")) 
})

annotate_figure(ggarrange(plotlist=glist,ncol=5,nrow=2,common.legend=T,legend="bottom"),
                         left=text_grob("Normalized expression",size=16,rot=90),
                         top=text_grob("Oligo\n",size=16,rot=0))





