setwd("~/moorea_holobiont/mr_2brad/3.pop_structure")
##just bams names:
#bams=read.table("bams")[,1] # list of bam files, in the same order as you did your analysis
##renamed the long bam file names for plotting:
bams=read.table("bams_renamedforplot")[,1] # list of bam files, in the same order as you did your analysis
goods=c(1:length(bams))

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
ma = as.matrix(read.table("myresult.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
quartz()
plot(hc,cex=0.3,xlab="") # this shows how similar clones are
#my output is named 'part3_clones_dendro.pdf' in github folder
#all of my 'd' (technical replicates) named ones are clones, as expected
#also going to get rid of two which were unintentionally clones: MNW-F_125 & MNW-B_58