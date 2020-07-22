#### piphillin input prep ####
#piphillin requires a .csv file of the counts + a .fasta of the sequences:
seq.rare.trim <- read.csv("~/moorea_holobiont/mr_16S/seq.rare12k.trim_rd2.csv",row.names=1)
sqs.rare.trim <- colnames(seq.rare.trim)
taxa.rare.trim <- data.frame(taxa2[(rownames(taxa2) %in% sqs.rare.trim),])
colnames(seq.rare.trim) == rownames(taxa.rare.trim) #double check
colnames(seq.rare.trim) <- taxa.rare.trim$V8
ids.rare.trim <- rownames(taxa.rare.trim)

seq.rare.trim.matrix <- as.matrix(seq.rare.trim)

#fasta file
library(dada2)
path='~/moorea_holobiont/mr16s_fxn/mr16s.rare.trim.fasta'
uniquesToFasta(seq.rare.trim.matrix, path, ids = ids.rare.trim, mode = "w", width = 20000)

#csv file
seq.rare.trim <- read.csv("~/moorea_holobiont/mr_16S/seq.rare12k.trim_rd2.csv",row.names=1)
seq.rare.trim.t <- t(seq.rare.trim)
write.csv(seq.rare.trim.t,file="seq.rare.trim.t_forpiphillin.csv")
#manually renamed the 'X' as 'sqs'

#### array quality metrics - pathways ####
#BiocManager::install("arrayQualityMetrics")
library("arrayQualityMetrics")
#BiocManager::install("DESeq")
library("DESeq")
#library("DESeq2")
#install.packages('systemfonts')

setwd("~/moorea_holobiont/mr16s_fxn")
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id
#count data
countData <- read.table("ko_pathway_abund_table_unnorm.txt",row.names=1,header=TRUE)
countData2 = countData[rowSums(countData)>1, ]
countData_rounded <- floor(countData2)

real=newCountDataSet(countData_rounded,colData) 
real=estimateSizeFactors(real)
plot(sort(sizeFactors(real))) 

cds=estimateDispersions(real,method="blind",fitType='local')
vsdBlind=varianceStabilizingTransformation(cds)
arrayQualityMetrics(vsdBlind,intgroup=c("zone"), force=TRUE, outdir = "~/arrayQualityMetrics_path") # this makes a directory "arrayQualityMetrics" and outputs various reports on outliers

#identified outliers:
outliers_paths <- c("A11", "B1", "B9", "C10", "D6", "D7", "E3", "H11", "H12")
#same ones as for ko, but without B6 this time.. 

colData2 <- colData[!(row.names(colData) %in% outliers_paths),]
countData_rounded2 <- countData_rounded[,!(colnames(countData_rounded) %in% outliers_paths)]

#re-naming the no outliers ones so I can use them below:
countData_rounded <- countData_rounded2
colData <- colData2

### deseq + pathways ####
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
#resources to go back to:
#deseq stuffs
#https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
#KEGGREST
#https://bioconductor.org/packages/release/bioc/html/KEGGREST.html

#with outliers - skip if you removed outliers above
setwd("~/moorea_holobiont/mr16s_fxn")
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id
#count data
countData <- read.table("ko_pathway_abund_table_unnorm.txt",row.names=1,header=TRUE)

# Filter data where you only have 0 or 1 read count across all samples.
countData2 = countData[rowSums(countData)>1, ]
countData_rounded <- floor(countData2)

# Set up the DESeqDataSet Object and run the DESeq pipeline
#add '2' for no outliers, remove '2' for with outliers
dds = DESeqDataSetFromMatrix(countData=countData_rounded,
                             colData=colData,
                             design=~site*zone)
dds = DESeq(dds)
dds

res = results(dds, contrast=c("zone","in", "out"))
res = res[order(res$pvalue),]
summary(res) #nothing lol

#trying by site now
#MOOREA NW
col.mnw <- subset(colData,site=="MNW")
id.mnw <- row.names(col.mnw)
id.mnw

count.mnw <- countData_rounded[,colnames(countData_rounded) %in% id.mnw]
count.mnw2 = count.mnw[rowSums(count.mnw)>1, ] #same

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.mnw = DESeqDataSetFromMatrix(countData=count.mnw2,
                             colData=col.mnw,
                             design=~zone)

dds.mnw = DESeq(dds.mnw)
dds.mnw

res.mnw = results(dds.mnw, contrast=c("zone","in", "out"))
res.mnw = res.mnw[order(res.mnw$pvalue),]
summary(res.mnw) #nothing + 8 low counts
#removing outlier samples didn't change anything :(

res.mnw.sig <- subset(res.mnw,padj<=0.1)
#nothing

#MOOREA SE
col.mse <- subset(colData,site=="MSE")
id.mse <- row.names(col.mse)

count.mse <- countData_rounded[,colnames(countData_rounded) %in% id.mse]
count.mse2 = count.mse[rowSums(count.mse)>1, ] #some

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.mse = DESeqDataSetFromMatrix(countData=count.mse2,
                                 colData=col.mse,
                                 design=~zone)

dds.mse = DESeq(dds.mse)
dds.mse

res.mse = results(dds.mse, contrast=c("zone","in", "out"))
res.mse = res.mse[order(res.mse$pvalue),]
summary(res.mse) #26 up, 21 down, 146 low counts
#goes down to 0 of everything with removing outlier samples 

res.mse.sig <- subset(res.mse,padj<=0.1)
saveRDS(res.mse.sig,file="paths.res.mse.sig.RDS")

#Tahiti NW
col.tnw <- subset(colData,site=="TNW")
id.tnw <- row.names(col.tnw)

count.tnw <- countData_rounded[,colnames(countData_rounded) %in% id.tnw]
count.tnw2 = count.tnw[rowSums(count.tnw)>1, ] #some

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.tnw = DESeqDataSetFromMatrix(countData=count.tnw2,
                                 colData=col.tnw,
                                 design=~zone)

dds.tnw = DESeq(dds.tnw)
dds.tnw

res.tnw = results(dds.tnw, contrast=c("zone","in", "out"))
res.tnw = res.tnw[order(res.tnw$pvalue),]
summary(res.tnw) #10 up, 9 down, 155 low counts

res.tnw.sig <- subset(res.tnw,padj<=0.1)
saveRDS(res.tnw.sig,file="paths.res.tnw.sig.RDS")

#### kegg name functions ####
#now what are these guys?
#BiocManager::install("KEGGREST")
library("KEGGREST")

kegg_name <- function(database=database,query=query){
  kegg_result <- keggFind(database=database,query=query)
  return(kegg_result)
}

#pathways function
kegg_info_paths <- function(res=res,database="pathway"){
  df <- data.frame(res)
  df$paths <- rownames(df)
  paths <- rownames(df)
  df.paths <- data.frame(paths)
  df.paths$paths.noko <- gsub("ko","",df.paths$paths) #kegg fxn needs just numbers
  df.paths$paths.noko <- as.character(df.paths$paths.noko)
  n <- length(df.paths[,1])
  for (i in 1:n)
  {
    call <- df.paths[i,2] 
    df.paths[i,3] <- kegg_name(database=database,query=call)
  }
  df.out <- merge(df,df.paths,by="paths")
  return(df.out)
}

#features function
kegg_info_feats <- function(res=res,database="ko"){
  df <- data.frame(res)
  df$paths <- rownames(df)
  paths <- rownames(df)
  df.paths <- data.frame(paths)
  n <- length(df.paths[,1])
  for (i in 1:n)
  {
    call <- df.paths[i,1] #first column must be kegg call name
    df.paths[i,2] <- kegg_name(database=database,query=call)
  }
  df.out <- merge(df,df.paths,by="paths")
  return(df.out)
}

#### running new functions ####

mnw.paths.id <- kegg_info_paths(res=res.mnw.sig)
#empty - doesn't work

mse.paths.id <- kegg_info_paths(res=res.mse.sig)
write.csv(mse.paths.id,file="mse.paths.id.csv")

tnw.paths.id <- kegg_info_paths(res=res.tnw.sig)
write.csv(tnw.paths.id,file="tnw.paths.id.csv")

#### read back in deseq results & paths ####
setwd("~/moorea_holobiont/mr16s_fxn")
res.mse.sig <- readRDS("paths.res.mse.sig.RDS")
res.tnw.sig <- readRDS("paths.res.tnw.sig.RDS")
mse.paths.id <- read.csv(file="mse.paths.id.csv",row.names = 1,header=TRUE)
tnw.paths.id <- read.csv(file="tnw.paths.id.csv",row.names = 1,header=TRUE)

#are there any that are the same between them?
incommon <- merge(mse.paths.id,tnw.paths.id,by="V3") #yesssssss, 15 
#oops, need to sort by fore reef & back reef

#sorting by pos & neg [back & fore]
paths.msei <- mse.paths.id[mse.paths.id$log2FoldChange>0,]
paths.mseo <- mse.paths.id[mse.paths.id$log2FoldChange<0,]

paths.tnwi <- tnw.paths.id[tnw.paths.id$log2FoldChange>0,]
paths.tnwo <- tnw.paths.id[tnw.paths.id$log2FoldChange<0,]

paths.msei$V3 <- as.character(paths.msei$V3)
paths.mseo$V3 <- as.character(paths.mseo$V3)
paths.tnwi$V3 <- as.character(paths.tnwi$V3)
paths.tnwo$V3 <- as.character(paths.tnwo$V3)

msei_path_id <- paths.msei$V3
mseo_path_id <- paths.mseo$V3
tnwi_path_id <- paths.tnwi$V3
tnwo_path_id <- paths.tnwo$V3

shared.path.msei.tnwi <- merge(paths.msei,paths.tnwi,by="paths") #none in common
shared.path.mseo.tnwo <- merge(paths.mseo,paths.tnwo,by="paths") #none in common
#nothing :(

#### heat map - pathways ####
library(pheatmap)

#Moorea SE
mse.ids <- as.character(mse.paths.id$paths)
#mse.counts <- counts(dds.mse,normalized=T)
rlog.mse.path <- rlogTransformation(dds.mse, blind=TRUE,fitType='local')
mse.counts <- assay(rlog.mse.path)
exp.mse <- mse.counts[rownames(mse.counts) %in% mse.ids,]

colnames(exp.mse) == col.mse$id
colnames(exp.mse) <- col.mse$fullname

rownames(exp.mse) == mse.paths.id$paths
rownames(exp.mse) <- mse.paths.id$V3

means=apply(exp.mse,1,mean) # means of rows
explc.mse=exp.mse-means # subtracting them

#plot
colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)

exp2.mse <- exp.mse[,order(colnames(exp.mse))]
explc2.mse <- explc.mse[,order(colnames(explc.mse))]

col_order <- c("MSEI_B9","MSEI_C10","MSEI_C9","MSEI_A10","MSEI_E9","MSEI_F10","MSEI_G10","MSEI_H9",
               "MSEI_D10", "MSEI_H10","MSEI_A9","MSEI_G9","MSEI_D9","MSEI_B10","MSEI_E10",
               "MSEO_F12","MSEO_C11","MSEO_E12","MSEO_D12","MSEO_B12","MSEO_H11","MSEO_G11",
               "MSEO_G12","MSEO_A12","MSEO_B11","MSEO_E11","MSEO_D11","MSEO_F11","MSEO_A11","MSEO_H12")

exp2.mse <- exp.mse[, col_order]

quartz()
pheatmap(exp2.mse,scale="row",color=colz,show_rownames=T, border_color = NA,
         cluster_cols=F,main="Mo'orea SE")
dev.off()

#Tahiti NW
tnw.ids <- as.character(tnw.paths.id$paths)
tnw.counts <- counts(dds.tnw,normalized=T)
exp.tnw <- tnw.counts[rownames(tnw.counts) %in% tnw.ids,]

colnames(exp.tnw) == col.tnw$id
colnames(exp.tnw) <- col.tnw$fullname

rownames(exp.tnw) == tnw.paths.id$paths
rownames(exp.tnw) <- tnw.paths.id$V3

colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)

col_order <- c("TI_G5","TI_A6","TI_G6","TI_C6","TI_B6","TI_C5","TI_F5","TI_D5",
               "TI_H6","TI_E6","TI_F6","TI_A5","TI_E5","TI_D6","TO_D7","TO_F7","TO_H7","TO_H8",
               "TO_B7","TO_F8","TO_C7","TO_E7","TO_G8","TO_D8","TO_C8","TO_G7","TO_A7")
  
exp2.tnw <- exp.tnw[, col_order]

quartz()
pheatmap(exp2.tnw,scale="row",color=colz,show_rownames=T, border_color = "grey",
         cluster_cols=F,main="Tahiti NW")
dev.off()

#### array quality metrics - ko ####
#BiocManager::install("arrayQualityMetrics")
library("arrayQualityMetrics")
#BiocManager::install("DESeq")
library("DESeq")
#library("DESeq2")
#install.packages('systemfonts')

setwd("~/moorea_holobiont/mr16s_fxn")
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id
#count data
countData <- read.table("ko_abund_table_unnorm.txt",row.names=1,header=TRUE)
countData2 = countData[rowSums(countData)>1, ]
countData_rounded <- floor(countData2)

real=newCountDataSet(countData_rounded,colData) 
real=estimateSizeFactors(real)
plot(sort(sizeFactors(real))) 

cds=estimateDispersions(real,method="blind",fitType='local')
vsdBlind=varianceStabilizingTransformation(cds)
arrayQualityMetrics(vsdBlind,intgroup=c("zone"), force=TRUE, outdir = "~/arrayQualityMetrics_ko") # this makes a directory "arrayQualityMetrics" and outputs various reports on outliers

#identified outliers:
outliers_ko <- c("A11", "B1","B6", "B9", "C10", "D6", "D7", "E3", "H11", "H12")

colData2 <- colData[!(row.names(colData) %in% outliers_ko),]
countData_rounded2 <- countData_rounded[,!(colnames(countData_rounded) %in% outliers_ko)]

#re-naming the no outliers ones so I can use them below:
countData_ko_rounded <- countData_rounded2
colData <- colData2

#### deseq + ko ####
#count data
#BiocManager::install("DESeq2")
library(DESeq2)

setwd("~/moorea_holobiont/mr16s_fxn")
#skip these if you're doing without outliers - made just above
countData_ko <- read.table("ko_abund_table_unnorm.txt",row.names=1,header=TRUE)
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id

# Filter data where you only have 0 or 1 read count across all samples.
countData2_ko = countData_ko[rowSums(countData_ko)>1, ] #no difference
countData_ko_rounded <- floor(countData2_ko)

#MOOREA NW
col.mnw <- subset(colData,site=="MNW")
id.mnw <- row.names(col.mnw)

count.ko.mnw <- countData_ko_rounded[,colnames(countData_ko_rounded) %in% id.mnw]
count.ko.mnw2 = count.ko.mnw[rowSums(count.ko.mnw)>1, ] #some fewer

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.ko.mnw = DESeqDataSetFromMatrix(countData=count.ko.mnw2,
                                 colData=col.mnw,
                                 design=~zone)

dds.ko.mnw = DESeq(dds.ko.mnw)
dds.ko.mnw

res.ko.mnw = results(dds.ko.mnw, contrast=c("zone","in", "out"))
res.ko.mnw = res.ko.mnw[order(res.ko.mnw$pvalue),]
summary(res.ko.mnw) #281 up, 88 down, 17% low counts

res.ko.mnw.sig <- subset(res.ko.mnw,padj<=0.1)

#different normalization
rlog.mnw <- rlogTransformation(dds.ko.mnw, blind=TRUE,fitType='local')
rlog.mnw.out <- assay(rlog.mnw)
saveRDS(rlog.mnw.out,file="rlog.mnw.ko.RDS")

#MOOREA SE
col.mse <- subset(colData,site=="MSE")
id.mse <- row.names(col.mse)

count.ko.mse <- countData_ko_rounded[,colnames(countData_ko_rounded) %in% id.mse]
count.ko.mse2 = count.ko.mse[rowSums(count.ko.mse)>1, ] #some

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.ko.mse = DESeqDataSetFromMatrix(countData=count.ko.mse2,
                                 colData=col.mse,
                                 design=~zone)

dds.ko.mse = DESeq(dds.ko.mse)

res.ko.mse = results(dds.ko.mse, contrast=c("zone","in", "out"))
res.ko.mse = res.ko.mse[order(res.ko.mse$pvalue),]
summary(res.ko.mse) #17 up, 207 down, 1168 low counts

res.ko.mse.sig <- subset(res.ko.mse,padj<=0.1)

#different normalization
rlog.mse <- rlogTransformation(dds.ko.mse, blind=TRUE,fitType='local')
rlog.mse.out <- assay(rlog.mse)
saveRDS(rlog.mse.out,file="rlog.mse.ko.RDS")

#Tahiti NW
col.tnw <- subset(colData,site=="TNW")
id.tnw <- row.names(col.tnw)

count.ko.tnw <- countData_ko_rounded[,colnames(countData_ko_rounded) %in% id.tnw]
count.ko.tnw2 = count.ko.tnw[rowSums(count.ko.tnw)>1, ] #some

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.ko.tnw = DESeqDataSetFromMatrix(countData=count.ko.tnw2,
                                 colData=col.tnw,
                                 design=~zone)

dds.ko.tnw = DESeq(dds.ko.tnw)
dds.ko.tnw

res.ko.tnw = results(dds.ko.tnw, contrast=c("zone","in", "out"))
res.ko.tnw = res.ko.tnw[order(res.ko.tnw$pvalue),]
summary(res.ko.tnw) #95 up, 9 down, 386 low counts

res.ko.tnw.sig <- subset(res.ko.tnw,padj<=0.1)

rlog.tnw <- rlogTransformation(dds.ko.tnw, blind=TRUE,fitType='local')
rlog.tnw.out <- assay(rlog.tnw)
saveRDS(rlog.tnw.out,file="rlog.tnw.ko.RDS")

#got an error at 01572 - removing for now
res.ko.mnw.sig2 <- res.ko.mnw.sig[!(row.names(res.ko.mnw.sig) %in% "K01572"), ]
mnw.ko.id <- get_kegg_info(database="ko",res=res.ko.mnw.sig2)

res.ko.mse.sig2 <- res.ko.mse.sig[!(row.names(res.ko.mse.sig) %in% "K02756"), ]
mse.ko.id <- get_kegg_info(database="ko",res=res.ko.mse.sig2)

tnw.ko.id <- get_kegg_info(database="ko",res=res.ko.tnw.sig)

merge(tnw.ko.id,mnw.ko.id,by="V2")

#need to save these dfs right now!
write.csv(mnw.ko.id,file="mnw.ko.id.csv")
write.csv(mse.ko.id,file="mse.ko.id.csv")
write.csv(tnw.ko.id,file="tnw.ko.id.csv")

#### read back in results for features ####
setwd("~/moorea_holobiont/mr16s_fxn")
mnw.ko.id <- read.csv(file="mnw.ko.id.csv",row.names=1,header=TRUE)
mse.ko.id <- read.csv(file="mse.ko.id.csv",row.names=1,header=TRUE)
tnw.ko.id <- read.csv(file="tnw.ko.id.csv",row.names=1,header=TRUE)

rlog.mnw.out <- readRDS("rlog.mnw.ko.RDS")
rlog.mse.out <- readRDS("rlog.mse.ko.RDS")
rlog.tnw.out <- readRDS("rlog.tnw.ko.RDS")

mnw.mse <- merge(mnw.ko.id,mse.ko.id,by="paths") #6 in common
mnw.tnw <- merge(mnw.ko.id,tnw.ko.id,by="paths") #7 in common
mse.tnw <- merge(mse.ko.id,tnw.ko.id,by="paths") #28 in common
#just realized I need to split by reef zone (pos or neg log change) first

#sorting by pos & neg [back & fore]
mnwi <- mnw.ko.id[mnw.ko.id$log2FoldChange>0,]
mnwo <- mnw.ko.id[mnw.ko.id$log2FoldChange<0,]

msei <- mse.ko.id[mse.ko.id$log2FoldChange>0,]
mseo <- mse.ko.id[mse.ko.id$log2FoldChange<0,]

tnwi <- tnw.ko.id[tnw.ko.id$log2FoldChange>0,]
tnwo <- tnw.ko.id[tnw.ko.id$log2FoldChange<0,]

mnwi$V4 <- as.character(mnwi$V4)
mnwo$V4 <- as.character(mnwo$V4)
msei$V4 <- as.character(msei$V4)
mseo$V4 <- as.character(mseo$V4)
tnwi$V4 <- as.character(tnwi$V4)
tnwo$V4 <- as.character(tnwo$V4)

mnwi_id <- mnwi$V4
mnwo_id <- mnwo$V4
msei_id <- msei$V4
mseo_id <- mseo$V4
tnwi_id <- tnwi$V4
tnwo_id <- tnwo$V4

shared.mnwi.msei <- merge(mnwi,msei,by="paths") #4 in common
shared.mnwo.mseo <- merge(mnwo,mseo,by="paths") #1 in common
shared.mnwi.tnwi <- merge(mnwi,tnwi,by="paths") #6 in common
mseo.tnwo <- merge(mseo,tnwo,by="paths") #just testing, supposed to be 0 & there is 0

#just by site
mnw <- mnw.ko.id$V4
mse <- mse.ko.id$V4
tnw <- tnw.ko.id$V4

#install.packages("VennDiagram")
library(VennDiagram)

#site
all_shared_site = list("mnw" = mnw, "msei" = mse, "tnwi"=tnw)
prettyvenn=venn.diagram(
  x = all_shared_site,
  filename=NULL,
  #main="Back reef",
  #col = "transparent",
  fill = c("darkslategray3","darkslategray4","#000004"),
  alpha = 0.5,
  #label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  # cex = 2.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  #cat.pos=1
);
quartz()
grid.draw(prettyvenn)
dev.off()

all_shared = list("mnwi" = mnwi_id, "msei" = msei_id, "tnwi"=tnwi_id)
prettyvenn=venn.diagram(
  x = all_shared,
  filename=NULL,
  #main="Back reef",
  #col = "transparent",
  fill = c("darkslategray3","darkslategray4","#000004"),
  alpha = 0.5,
  #label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  # cex = 2.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos=1
);
grid.draw(prettyvenn)
dev.off()

all_shared_off = list("mnwo" = mnwo_id, "mseo" = mseo_id, "tnwo"=tnwo_id)
prettyvenn_off=venn.diagram(
  x = all_shared_off,
  filename=NULL,
  #main="Fore reef",
  #col = "transparent",
  fill = c("darkslategray3","darkslategray4","#000004"),
  alpha = 0.5,
  #label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  # cex = 2.5,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos=1
);

grid.draw(prettyvenn_off)
dev.off()

library(gridExtra)
quartz()
grid.arrange(gTree(children=prettyvenn), gTree(children=prettyvenn_off),nrow=2)

#### heat map - ko ####
library(pheatmap)
library(gplots)

#MOOREA NW
#I want fewer, or else it's too much:
mnw.order <- mnw.ko.id[order(mnw.ko.id$padj),]
mnw.top50 <- mnw.order[1:20,]
#mnw.p <- filter(mnw.ko.id, padj <= 0.01)

mnw.ids <- as.character(mnw.top50$paths)
#mnw.counts <- counts(dds.ko.mnw,normalized=T)
mnw.counts <- rlog.mnw.out #different method of normalizing
exp.mnw <- mnw.counts[rownames(mnw.counts) %in% mnw.ids,]

mnw.top50 <- mnw.top50[order(mnw.top50$paths),]
mnw.top50$paths == rownames(exp.mnw)
rownames(exp.mnw) <- mnw.top50$V4
colnames(exp.mnw) == col.mnw$id
colnames(exp.mnw) <- col.mnw$fullname

colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)

exp2.mnw <- exp.mnw[,order(colnames(exp.mnw))]
quartz()
pheatmap(exp2.mnw,scale="row",color=colz,
         border_color="grey",show_rownames=T, 
         cluster_cols=F,main="Mo'orea NW")
dev.off()

#means method
means=apply(exp2.mnw,1,mean) # means of rows
explc=exp2.mnw-means # subtracting them

col_order <- c("MNWI_D1","MNWI_B2","MNWI_E1","MNWI_C1","MNWI_G2","MNWI_B1","MNWI_A1",
               "MNWI_G1","MNWI_F2","MNWI_E2","MNWI_F1","MNWI_H2","MNWI_A2","MNWI_H1",
               "MNWO_E4","MNWO_E3","MNWO_F4","MNWO_C4","MNWO_A3","MNWO_G4","MNWO_C3",
               "MNWO_D3","MNWO_F3","MNWO_B4","MNWO_B3","MNWO_G3","MNWO_H3")

exp2.mnw <- exp2.mnw[, col_order]

zone_col <- data.frame(c(rep("BR",14),rep("FR",13)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp2.mnw)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

pheatmap(exp2.mnw,color=colz,
         border_color="grey",show_rownames=T, 
         cluster_cols=F,main="Mo'orea NW",scale="row",
         annotation_col=zone_col,
         annotation_colors=my_colour,
         show_colnames=FALSE)

#other heat map
# heatmap.2(as.matrix(explc), col = colz, Rowv = TRUE, Colv = FALSE, scale = "row",
#           dendrogram = "both",
#           trace = "none",
#           main = "GO:0004298;GO:0070003 threonine-type endopeptidase activity",
#           margin = c(5,15))

#MOOREA SE
#I want fewer, or else it's too much:
mse.order <- mse.ko.id[order(mse.ko.id$padj),]
mse.top50 <- mse.order[1:20,]
#mse.p <- filter(mse.ko.id, padj <= 0.01)

mse.ids <- as.character(mse.top50$paths)
#mse.counts <- counts(dds.ko.mse,normalized=T)
mse.counts <- rlog.mse.out #different method of normalizing
exp.mse <- mse.counts[rownames(mse.counts) %in% mse.ids,]

mse.top50 <- mse.top50[order(mse.top50$paths),]
mse.top50$paths == rownames(exp.mse)
rownames(exp.mse) <- mse.top50$V4
colnames(exp.mse) == col.mse$id
colnames(exp.mse) <- col.mse$fullname

colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)

exp2.mse <- exp.mse[,order(colnames(exp.mse))]
quartz()
pheatmap(exp2.mse,scale="row",color=colz,
         border_color="grey",show_rownames=T, 
         cluster_cols=T,main="Mo'orea SE")
dev.off()

#means method
means=apply(exp2.mse,1,mean) # means of rows
explc=exp2.mse-means # subtracting them

col_order <- c("MSEI_G10","MSEI_H10","MSEI_H9","MSEI_C10","MSEI_D10","MSEI_G9","MSEI_A9",
"MSEI_C9","MSEI_D9","MSEI_A10","MSEI_F10","MSEI_B9","MSEI_E9","MSEI_B10","MSEI_E10",
"MSEO_A11","MSEO_H12","MSEO_H11","MSEO_D12","MSEO_B12","MSEO_D11","MSEO_A12","MSEO_B11",
"MSEO_E11","MSEO_E12","MSEO_F12","MSEO_C11","MSEO_F11","MSEO_G11","MSEO_G12")
  
exp2.mse <- exp2.mse[, col_order]

zone_col <- data.frame(c(rep("BR",15),rep("FR",15)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp2.mse)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
heat.mse <- pheatmap(exp2.mse,color=colz,
                     border_color="grey",show_rownames=T, 
                     cluster_cols=F,main="Mo'orea SE",scale="row",
                     annotation_col=zone_col,
                     annotation_colors=my_colour,
                     show_colnames=FALSE)

#TAHITI NW
#I want fewer, or else it's too much:
tnw.order <- tnw.ko.id[order(tnw.ko.id$padj),]
tnw.top50 <- tnw.order[1:20,]
#tnw.p <- filter(tnw.ko.id, padj <= 0.01)

tnw.ids <- as.character(tnw.top50$paths)
#tnw.counts <- counts(dds.ko.tnw,normalized=T)
tnw.counts <- rlog.tnw.out #different method of normalizing
exp.tnw <- tnw.counts[rownames(tnw.counts) %in% tnw.ids,]

tnw.top50 <- tnw.top50[order(tnw.top50$paths),]
tnw.top50$paths == rownames(exp.tnw)
rownames(exp.tnw) <- tnw.top50$V4
colnames(exp.tnw) == col.tnw$id
colnames(exp.tnw) <- col.tnw$fullname

colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)

exp2.tnw <- exp.tnw[,order(colnames(exp.tnw))]
quartz()
pheatmap(exp2.tnw,scale="row",color=colz,
         border_color="grey",show_rownames=T, 
         cluster_cols=T,main="Mo'orea SE")
dev.off()

#means method
means=apply(exp2.tnw,1,mean) # means of rows
explc=exp2.tnw-means # subtracting them

col_order <- c("tnwI_G10","tnwI_H10","tnwI_H9","tnwI_C10","tnwI_D10","tnwI_G9","tnwI_A9",
               "tnwI_C9","tnwI_D9","tnwI_A10","tnwI_F10","tnwI_B9","tnwI_E9","tnwI_B10","tnwI_E10",
               "tnwO_A11","tnwO_H12","tnwO_H11","tnwO_D12","tnwO_B12","tnwO_D11","tnwO_A12","tnwO_B11",
               "tnwO_E11","tnwO_E12","tnwO_F12","tnwO_C11","tnwO_F11","tnwO_G11","tnwO_G12")

exp2.tnw <- exp2.tnw[, col_order]

zone_col <- data.frame(c(rep("BR",15),rep("FR",15)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp2.tnw)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp2.tnw,color=colz,
         border_color="grey",show_rownames=T, 
         cluster_cols=F,main="Mo'orea SE",scale="row",
         annotation_col=zone_col,
         annotation_colors=my_colour,
         show_colnames=FALSE)

#### heat maps - shared ones ####
library(pheatmap)

#backreef first
sh.mnwi.tnwi.ids <- as.character(shared.mnwi.tnwi$paths)
counts.tnw <- rlog.tnw.out #different method of normalizing
exp.tnwi.sh <- counts.tnw[rownames(counts.tnw) %in% sh.mnwi.tnwi.ids,]

shared.mnwi.tnwi$paths == rownames(exp.tnwi.sh)
rownames(exp.tnwi.sh) <- shared.mnwi.tnwi$V4.x

colnames(exp.tnwi.sh) == col.tnw$id
colnames(exp.tnwi.sh) <- col.tnw$fullname

# col_order <- c("TI_F5","TI_A5","TI_C5","TI_C6","TI_G6","TI_D6","TI_E5","TI_B6","TI_D5",
#                "TI_H6","TI_G5","TI_A6","TI_E6","TI_F6",
#                "TO_D7","TO_A7","TO_G7","TO_H7","TO_C8","TO_G8","TO_C7","TO_E7","TO_B7",
#                "TO_F8","TO_D8","TO_F7","TO_H8")

#exp.tnwi.sh2 <- exp.tnwi.sh[, col_order]

exp.tnwi.sh2 <- exp.tnwi.sh[,order(colnames(exp.tnwi.sh))]

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#8405A7FF","white","#ED7953FF"))(100)

zone_col <- data.frame(c(rep("BR",14),rep("FR",13)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.tnwi.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.tnwi.sh2,scale="row",color=colz_rz,
         border_color="grey",show_rownames=T,show_colnames=F, 
         cluster_cols=F,main="TNW-B",
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()

counts.mnw <- rlog.mnw.out #different method of normalizing
exp.mnwi.sh <- counts.mnw[rownames(counts.mnw) %in% sh.mnwi.tnwi.ids,]

shared.mnwi.tnwi$paths == rownames(exp.mnwi.sh)
rownames(exp.mnwi.sh) <- shared.mnwi.tnwi$V4.x

colnames(exp.mnwi.sh) == col.mnw$id
colnames(exp.mnwi.sh) <- col.mnw$fullname

# col_order <- c("MNWI_B2","MNWI_D1","MNWI_B1","MNWI_E2","MNWI_G2","MNWI_C1","MNWI_G1",
#                "MNWI_H1","MNWI_H2","MNWI_A2","MNWI_A1","MNWI_E1","MNWI_F1","MNWI_F2",
#                "MNWO_E4","MNWO_G3","MNWO_B4","MNWO_C3","MNWO_H3","MNWO_F3","MNWO_E3",
#                "MNWO_A3","MNWO_F4","MNWO_G4","MNWO_B3","MNWO_C4","MNWO_D3")
# 
# exp.mnwi.sh2 <- exp.mnwi.sh[, col_order]

exp.mnwi.sh2 <- exp.mnwi.sh[,order(colnames(exp.mnwi.sh))]

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#8405A7FF","white","#ED7953FF"))(100)

zone_col <- data.frame(c(rep("BR",14),rep("FR",13)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.mnwi.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.mnwi.sh2,scale="row",color=colz_rz,
         border_color=NA,show_rownames=T,show_colnames = F, 
         cluster_cols=F,main="MNW-B",
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()

#fore reef next
sh.mnwo.mseo.ids <- as.character(shared.mnwo.mseo$paths)
counts.mse <- rlog.mse.out #different method of normalizing
exp.mseo.sh <- data.frame(counts.mse[rownames(counts.mse) %in% sh.mnwo.mseo.ids,])
exp.mseo.sh <- unname(exp.mseo.sh)
colnames(exp.mseo.sh) <- "K01270"
exp.mseo.sh <- t(exp.mseo.sh)

shared.mnwo.mseo$paths == rownames(exp.mseo.sh)
rownames(exp.mseo.sh) <- shared.mnwo.mseo$V4.x

colnames(exp.mseo.sh) == col.mse$id
colnames(exp.mseo.sh) <- col.mse$fullname

# col_order <- c("MSEI_C9","MSEI_D9","MSEI_F10","MSEI_H9","MSEI_C10","MSEI_E9","MSEI_D10",
#                "MSEI_G10","MSEI_E10","MSEI_H10","MSEI_A10","MSEI_B10","MSEI_G9","MSEI_B9",
#                "MSEI_A9","MSEO_H12","MSEO_G11","MSEO_E11","MSEO_G12","MSEO_E12","MSEO_D12",
#                "MSEO_B12","MSEO_F12","MSEO_C11","MSEO_D11","MSEO_H11","MSEO_A12","MSEO_B11",
#                "MSEO_A11","MSEO_F11")

#exp.mseo.sh2 <- t(data.frame(exp.mseo.sh[, col_order]))
exp.mseo.sh2 <- t(data.frame(exp.mseo.sh[,order(colnames(exp.mseo.sh))]))

rownames(exp.mseo.sh2) <- c("pepD; dipeptidase D [EC:3.4.13.-]")

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#ED7953FF","white","#8405A7FF"))(100)

zone_col <- data.frame(c(rep("BR",15),rep("FR",15)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.mseo.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.mseo.sh2,scale="row",color=colz_rz,
         border_color=NA,show_rownames=T,show_colnames=F, 
         cluster_cols=F,main="MSE-F",cluster_rows=F,
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()

counts.mnw <- rlog.mnw.out #different method of normalizing
exp.mnwo.sh <- data.frame(counts.mnw[rownames(counts.mnw) %in% sh.mnwo.mseo.ids,])
exp.mnwo.sh <- unname(exp.mnwo.sh)
colnames(exp.mnwo.sh) <- "K01270"
exp.mnwo.sh <- t(exp.mnwo.sh)

shared.mnwo.mseo$paths == rownames(exp.mnwo.sh)
rownames(exp.mnwo.sh) <- shared.mnwo.mseo$V4.x

colnames(exp.mnwo.sh) == col.mnw$id
colnames(exp.mnwo.sh) <- col.mnw$fullname

#exp.mnwo.sh2 <- t(data.frame(exp.mnwo.sh[,order(colnames(exp.mnwo.sh))]))
exp.mnwo.sh2 <- t(data.frame(exp.mnwo.sh[,order(colnames(exp.mnwo.sh))]))

rownames(exp.mnwo.sh2) <- c("pepD; dipeptidase D [EC:3.4.13.-]")

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#ED7953FF","white","#8405A7FF"))(100)

zone_col <- data.frame(c(rep("BR",14),rep("FR",13)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.mnwo.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.mnwo.sh2,scale="row",color=colz_rz,
         border_color=NA,show_rownames=T,show_colnames=F, 
         cluster_cols=F,main="MNW-F",cluster_rows=F,
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()

#last one! 4 shared between MNW-B & MSE-B
sh.mnwi.msei.ids <- as.character(shared.mnwi.msei$paths)
counts.mse <- rlog.mse.out #different method of normalizing
exp.msei.sh <- counts.mse[rownames(counts.mse) %in% sh.mnwi.msei.ids,]

shared.mnwi.msei$paths == rownames(exp.msei.sh)
rownames(exp.msei.sh) <- shared.mnwi.msei$V4.x

colnames(exp.msei.sh) == col.mse$id
colnames(exp.msei.sh) <- col.mse$fullname

# col_order <- c("TI_F5","TI_A5","TI_C5","TI_C6","TI_G6","TI_D6","TI_E5","TI_B6","TI_D5",
#                "TI_H6","TI_G5","TI_A6","TI_E6","TI_F6",
#                "TO_D7","TO_A7","TO_G7","TO_H7","TO_C8","TO_G8","TO_C7","TO_E7","TO_B7",
#                "TO_F8","TO_D8","TO_F7","TO_H8")

#exp.tnwi.sh2 <- exp.tnwi.sh[, col_order]

exp.msei.sh2 <- exp.msei.sh[,order(colnames(exp.msei.sh))]

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#8405A7FF","white","#ED7953FF"))(100)

zone_col <- data.frame(c(rep("BR",15),rep("FR",15)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.msei.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.msei.sh2,scale="row",color=colz_rz,
         border_color="grey",show_rownames=T,show_colnames=F, 
         cluster_cols=F,main="MSE-B",
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()

counts.mnw <- rlog.mnw.out #different method of normalizing
exp.mnwi.sh <- counts.mnw[rownames(counts.mnw) %in% sh.mnwi.msei.ids,]

shared.mnwi.msei$paths == rownames(exp.mnwi.sh)
rownames(exp.mnwi.sh) <- shared.mnwi.msei$V4.x

colnames(exp.mnwi.sh) == col.mnw$id
colnames(exp.mnwi.sh) <- col.mnw$fullname

# col_order <- c("MNWI_B2","MNWI_D1","MNWI_B1","MNWI_E2","MNWI_G2","MNWI_C1","MNWI_G1",
#                "MNWI_H1","MNWI_H2","MNWI_A2","MNWI_A1","MNWI_E1","MNWI_F1","MNWI_F2",
#                "MNWO_E4","MNWO_G3","MNWO_B4","MNWO_C3","MNWO_H3","MNWO_F3","MNWO_E3",
#                "MNWO_A3","MNWO_F4","MNWO_G4","MNWO_B3","MNWO_C4","MNWO_D3")
# 
# exp.mnwi.sh2 <- exp.mnwi.sh[, col_order]

exp.mnwi.sh2 <- exp.mnwi.sh[,order(colnames(exp.mnwi.sh))]

#colz = colorRampPalette(c("#440154","#365D8D","white","#58C765","#FDE725"))(100)
colz_rz = colorRampPalette(c("#8405A7FF","white","#ED7953FF"))(100)

zone_col <- data.frame(c(rep("BR",14),rep("FR",13)))
colnames(zone_col) <- "Reef zone"
row.names(zone_col) <- colnames(exp.mnwi.sh2)

my_colour = list(
  "Reef zone" = c(BR = "#ED7953FF", FR = "#8405A7FF")
)

quartz()
pheatmap(exp.mnwi.sh2,scale="row",color=colz_rz,
         border_color=NA,show_rownames=T,show_colnames = F, 
         cluster_cols=F,main="MNW-B",
         annotation_col=zone_col,
         annotation_colors=my_colour)
dev.off()



#### archive ####
#### deseq without outliers ####
library("DESeq2")
library("dplyr")
library("ggplot2")

setwd("~/moorea_holobiont/mr16s_fxn")
countData_ko <- read.table("ko_abund_table_unnorm.txt",row.names=1,header=TRUE)
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
mnw.ko.id <- read.csv(file="mnw.ko.id.csv",row.names=1,header=TRUE)
row.names(colData) <- colData$id

# Filter data where you only have 0 or 1 read count across all samples.
countData2_ko = countData_ko[rowSums(countData_ko)>1, ] #no difference

countData_ko_rounded <- floor(countData2_ko)

#outliers identified in part 1:
outliers <- c("A11", "B1", "B6", "B9", "C10", "D6", "D7", "E3", "H11", "H12")
colData2 <- colData[!(row.names(colData) %in% outliers),]
countData_ko_rounded2 <- countData_ko_rounded[,!(colnames(countData_ko_rounded) %in% outliers)]

#MOOREA NW
col.mnw <- subset(colData2,site=="MNW")
id.mnw <- row.names(col.mnw)

count.ko.mnw <- countData_ko_rounded2[,colnames(countData_ko_rounded2) %in% id.mnw]
count.ko.mnw2 = count.ko.mnw[rowSums(count.ko.mnw)>1, ] #some fewer

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.ko.mnw = DESeqDataSetFromMatrix(countData=count.ko.mnw2,
                                    colData=col.mnw,
                                    design=~zone)

dds.ko.mnw = DESeq(dds.ko.mnw)
dds.ko.mnw

res.ko.mnw = results(dds.ko.mnw, contrast=c("zone","in", "out"))
res.ko.mnw = res.ko.mnw[order(res.ko.mnw$pvalue),]
summary(res.ko.mnw) #without outliers: 272 br, 8 fr

res.ko.mnw.sig <- subset(res.ko.mnw,padj<=0.1)

#different normalization
rlog.mnw <- rlogTransformation(dds.ko.mnw, blind=TRUE,fitType='local')
rlogged.mnw <- assay(rlog.mnw)

DESeq2::plotPCA(rlog.mnw, returnData = TRUE, intgroup = c("zone") ) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = zone)) #+
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #labs(title = "~ genotype, p<0.001")+
  #theme_cowplot()

#heat map
library("gplots")
library("pheatmap")

res.mnw.df <- data.frame(res.ko.mnw.sig)
mnw.order <- res.mnw.df[order(res.mnw.df$padj),]
mnw.top50 <- mnw.order[1:50,]
mnw.top50.ids <- rownames(mnw.top50)

rlogged.mnw50 <- rlogged.mnw[(row.names(rlogged.mnw) %in% mnw.top50.ids),]

means=apply(rlogged.mnw50,1,mean) # means of rows
explc.mnw=rlogged.mnw50-means # subtracting them

colnames(explc.mnw) == col.mnw$id
colnames(explc.mnw) <- col.mnw$fullname

##need to update (if doing this analysis)
#rownames(explc.mnw) == mnw.ko.id$paths
#rownames(explc.mnw) <- mnw.ko.id$V3

heatmap.2(as.matrix(explc.mnw), Rowv = TRUE, Colv = TRUE, scale = "row",
          dendrogram = "both",
          trace = "none",
          main = "site",
          margin = c(5,15))
dev.off()

colnames(rlogged.mnw50) == col.mnw$id
colnames(rlogged.mnw50) <- col.mnw$fullname

pheatmap(rlogged.mnw50,scale="row",
         border_color="grey",show_rownames=T, 
         cluster_cols=T,main="Mo'orea NW")


