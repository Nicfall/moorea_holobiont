#### angsd call in terminal ####
#FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 91 -snp_pval 1e-5 -minMaf 0.05 -minIndDepth 8"
#TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doVcf 1 -doPost 1 -doGlf 2 -dumpCounts 2"
#angsd -b bamscl_no7 -GL 1 $FILTERS $TODO -P 1 -out clresult_no7.05.nohwe

#### pcangsd in terminal ####
#python pcangsd.py -beagle clresult_no8.05.nohwe -o beagle_out -threads 4

# perform PCA on covarince matrix in R
## open R
setwd("/Users/nicolakriefall/moorea_holobiont/mr_2brad/3.pop_structure/tests")
C <- as.matrix(read.table("beagle_out.cov"))
e <- eigen(C)

i2p=read.table("/Users/nicolakriefall/moorea_holobiont/mr_2brad/3.pop_structure/bamscl_no7_pops.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
site_zone=i2p[,2]
site=i2p[,3]
zone=i2p[,4]

palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

conds=data.frame(cbind(site_zone,site,zone))

plot(e$vectors[,1:2],xlab="PC1",ylab="PC2",col=transp(colors))
ordiellipse(e$vectors[,1:2],groups=conds$site,draw="polygon",col=colpops,label=T)

adonis(C~site*zone,data=conds)