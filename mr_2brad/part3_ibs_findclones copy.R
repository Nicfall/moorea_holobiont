bams=read.table("~/Google Drive/Moorea/2brad_moorea/part3_bams")[,1] # list of bam files, in the same order as you did your analysis
goods=c(1:length(bams))

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
ma = as.matrix(read.table("~/Google Drive/Moorea/2brad_moorea/part3_myresult.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
quartz()
plot(hc,cex=0.3,xlab="") # this shows how similar clones are
#my output is named 'part3_clones_dendro.pdf' in github folder
#all of my 'd' (technical replicates) named ones are clones, as expected
#also going to get rid of two which were unintentionally clones: MNW-F_125 & MNW-B_58