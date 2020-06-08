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

### deseq + kegg ####
#resources to go back to:
#deseq stuffs
#https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
#KEGGREST
#https://bioconductor.org/packages/release/bioc/html/KEGGREST.html

setwd("~/moorea_holobiont/mr16s_fxn")
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id

#### pathways ####
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)

#count data
countData <- read.table("ko_pathway_abund_table_unnorm.txt",row.names=1,header=TRUE)

# Filter data where you only have 0 or 1 read count across all samples.
countData2 = countData[rowSums(countData)>1, ]

countData_rounded <- floor(countData2)

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData_rounded,
                             colData=colData,
                             design=~site*zone)
dds = DESeq(dds)
dds

res = results(dds, contrast=c("zone","in", "out"))
res = res[order(res$pvalue),]
summary(res) #nothing lol

dds.all <- estimateSizeFactors(dds,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")
res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
#sigtab.all = cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.all), ], "matrix"))
sigtab.all
dim(sigtab.all) #nothing

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

res.mnw.sig <- subset(res.mnw,padj<=0.1)

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

res.mse.sig <- subset(res.mse,padj<=0.1)

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

#setting up data frames for calling the ko pathway names
mnw.paths <- rownames(res.mnw.sig)
df.mnw.paths <- data.frame(mnw.paths)
df.mnw.paths$mnw.paths <- gsub("ko","",df.mnw.paths$mnw.paths) #kegg fxn needs just numbers
df.mnw.paths$mnw.paths <- as.character(df.mnw.paths$mnw.paths)

mnw.paths.id <- get_kegg_info(df=df.mnw.paths,n=19,database="pathway")
#empty - doesn't work

mse.paths <- rownames(res.mse.sig)
df.mse.paths <- data.frame(mse.paths)
df.mse.paths$mse.paths <- gsub("ko","",df.mse.paths$mse.paths) #kegg fxn needs just numbers
df.mse.paths$mse.paths <- as.character(df.mse.paths$mse.paths)

mse.paths.id <- get_kegg_info(df=df.mse.paths,n=19,database="pathway")

tnw.paths <- rownames(res.tnw.sig)
df.tnw.paths <- data.frame(tnw.paths)
df.tnw.paths$tnw.paths <- gsub("ko","",df.tnw.paths$tnw.paths) #kegg fxn needs just numbers
df.tnw.paths$tnw.paths <- as.character(df.tnw.paths$tnw.paths)

tnw.paths.id <- get_kegg_info(df=df.tnw.paths,n=19,database="pathway")

#are there any that are the same between the two?
merge(mse.paths.id,tnw.paths.id) #yesssssss
# V2 mse.paths tnw.paths
# 1   Fatty acid degradation     00071     00071
# 2    Fatty acid metabolism     01212     01212
# 3 Homologous recombination     03440     03440
# 4        Purine metabolism     00230     00230
# 5    Pyrimidine metabolism     00240     00240
# 6                 Ribosome     03010     03010
# 7      Tyrosine metabolism     00350     00350

#### just 'features' ####
#count data
#BiocManager::install("DESeq2")
library(DESeq2)
setwd("~/moorea_holobiont/mr16s_fxn")
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

#venn diagram
#BiocManager::install()
#BiocManager::install("limma")
library(limma)
library(tidyverse)
library(ggforce)
set.seed((123))
mydata <- data.frame(A = rbinom(100, 1, 0.8),
                     B = rbinom(100, 1, 0.7),
                     C = rbinom(100, 1, 0.6)) %>%
  mutate_all(., as.logical)

df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = c('A', 'B', 'C'))

vdc <- vennCounts(mydata)
class(vdc) <- 'matrix'
df.vdc <- as.data.frame(vdc)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))

df.vdc$Counts2 <- c(1,2,3,4,5,6,7)

ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
  scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts2, size = 5)

#1 = MNW
#2 = TNW
#3 = MNW + TNW
#4 = MSE
#5 = MNW + MSE
#6 = MSE + TNW
#7 = all!

#### heat map ####
setwd("~/moorea_holobiont/mr16s_fxn")
countData_ko <- read.table("ko_abund_table_unnorm.txt",row.names=1,header=TRUE)
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_samdf.rare_12k.csv")
row.names(colData) <- colData$id

mnw.ids <- as.character(mnw.ko.id$paths)
mnw.counts <- counts(dds.ko.mnw)
exp <- mnw.counts[rownames(mnw.counts) %in% mnw.ids,]
#mnw.ko.id$paths == rownames(hm)
#rownames(hm) <- mnw.ko.id$V4
#colnames(hm) <- col.mnw$fullname

# iso2go = read_tsv('astrangia_iso2go.tab') %>%
#   rename(Iso = Gene_id)
# cold_results_df = as.data.frame(cold_results) %>%
#   rownames_to_column(var = 'Iso')
# hot_results_df = as.data.frame(hot_results) %>%
#   rownames_to_column(var = 'Iso')
# cold_rlog = as.data.frame(assay(cold_rlogged)) %>%
#   rownames_to_column(var = 'Iso') %>%
#   left_join(cold_results_df) %>%
#   filter(padj < 0.1) %>%
#   dplyr::select(-baseMean, -log2FoldChange, -lfcSE, -stat, -pvalue, -padj)
# hot_rlog = as.data.frame(assay(hot_rlogged)) %>%
#   rownames_to_column(var = 'Iso') %>%
#   left_join(hot_results_df) %>%
#   filter(padj < 0.1) %>%
#   dplyr::select(-baseMean, -log2FoldChange, -lfcSE, -stat, -pvalue, -padj)
# 
hot_colour = colorRampPalette(rev(c("#EA6227", "#F09167", "white", "grey40", "black")))(100)
                                     GO_hot = iso2go %>%
                                       filter(str_detect(GO_id, "GO:0004298") | str_detect(GO_id, "GO:0070003")) %>%
#   left_join(gene) %>%
#   left_join(hot_rlog) %>%
#   mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
#   column_to_rownames(var = "gene_symbol") %>%
#   dplyr::select(-GO_id, -Gene, -Iso) %>%
#   drop_na()%>%
#   dplyr::select(sort(current_vars()))
# exp = GO_hot

means=apply(exp,1,mean) # means of rows

explc=exp-means # subtracting them
library(gplots)

hot_colour = colorRampPalette(rev(c("#EA6227", "#F09167", "white", "grey40", "black")))(100)

heatmap.2(as.matrix(explc), col = hot_colour, Rowv = TRUE, Colv = FALSE, scale = "row",
          dendrogram = "both",
          trace = "none")
          #main = "GO:0004298;GO:0070003 threonine-type endopeptidase activity",
          #margin = c(5,15))
dev.off()

#alternative
pheatmap(trans,annotation_col=ann,annotation_colors=anno_colors,cex=1.2,border_color="grey",show_rownames=T, cluster_cols=F, color = colorRampPalette(c("white", "coral2"))(25))
#install.packages("pheatmap")
library('pheatmap')
pheatmap(explc,scale="column",
         cex=1.2,border_color="grey",show_rownames=T, 
         cluster_cols=T)

#### outdated heat map ####

## Turning into z-score table
hm2 <- counts(dds.ko.mnw)
dev.off()

hm.z = data.matrix(hm2)
hm.z = sweep(hm.z, 1L, rowMeans(hm.z), check.margin = FALSE)
hm.z.sx = apply(hm.z, 1L, sd)
hm.z = sweep(hm.z, 1L, hm.z.sx, "/", check.margin = FALSE)
hm.z = data.matrix(hm.z)

hm.z <- hm.z[rownames(hm.z) %in% mnw.ids,]
rownames(hm.z) <- mnw.ko.id$V4
colnames(hm.z) <- col.mnw$zone

top50 <- hm.z[1:50,]
top50 <- data.matrix(order(top50))

colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
heatmap.2(top50, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
          dendrogram = "both",
          trace = "none")
dev.off()





