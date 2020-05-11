### deseq + kegg? ####
#resources to go back to:
#deseq stuffs
#https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
#KEGGREST
#https://bioconductor.org/packages/release/bioc/html/KEGGREST.html

setwd("~/moorea_holobiont/mr16s_fxn")
colData <- read.csv("~/moorea_holobiont/mr_16S/mr16s_sampledata.csv")
row.names(colData) <- colData$id

# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)

#count data
countData <- read.table("ko_pathway_abund_table_unnorm copy.txt",row.names=1,header=TRUE)

# Filter data where you only have 0 or 1 read count across all samples.
countData2 = countData[rowSums(countData)>1, ]
head(countData2)

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
count.mnw2 = count.mnw[rowSums(count.mnw)>1, ] #some

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds.mnw = DESeqDataSetFromMatrix(countData=count.mnw2,
                             colData=col.mnw,
                             design=~zone)

dds.mnw = DESeq(dds.mnw)
dds.mnw

res.mnw = results(dds.mnw, contrast=c("zone","in", "out"))
res.mnw = res.mnw[order(res.mnw$pvalue),]
summary(res.mnw) #nothing lol

#MOOREA SE
col.mse <- subset(colData,site=="MSE")
id.mse <- row.names(col.mse)
id.mse

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
summary(res.mse) #things!

#Tahiti NW
col.tnw <- subset(colData,site=="TNW")
id.tnw <- row.names(col.tnw)
id.tnw

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
summary(res.tnw) #lots

#now what are these guys?
#BiocManager::install("KEGGREST")
library("KEGGREST")

keggFind(database="pathway",query="00572")


#BiocManager::install("gage")
library(gage)
BiocManager::install(c("pathview", "gage", "gageData", "GenomicAlignments",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
3
library(GenomicAlignments)
fls <- list.files("tophat_all/", pattern="bam$", full.names =T)
bamfls <- BamFileList(fls)
flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
param <- ScanBamParam(flag=flag)
#to run multiple core option: library(parallel); options("mc.cores"=4)
  gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="Union",
                               + ignore.strand=TRUE, singleEnd=FALSE, param=param)
hnrnp.cnts=assay(gnCnt)

data(hnrnp.cnts)
cnts=hnrnp.cnts
dim(cnts)
[1] 22932 8
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
[1] 17192 8
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
[1] 0.0 999179.9
cnts.norm=log2(cnts.norm+8)

ref.idx=5:8
samp.idx=1:4
data(kegg.gs)
cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
                      samp = samp.idx, compare ="unpaired")
##step 4: pathview
  cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
  + !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &
  + !is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
library(pathview)
pv.out.list <- sapply(path.ids2, function(pid) pathview(
  + gene.data = cnts.d, pathway.id = pid,
  + species = "hsa"))