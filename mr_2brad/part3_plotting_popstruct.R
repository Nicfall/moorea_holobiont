setwd("~/moorea_holobiont/mr_2brad")

#### PCA from .vcf file ####

library(vcfR)
library(adegenet)
library(vegan) 
library(stringr)
#install.packages("poppr")
library("poppr")

gl=vcfR2genlight(read.vcfR("donresult.05.vcf.gz")) #output from angsd
#gl=vcfR2genlight(read.vcfR("~/Google Drive/Moorea/2brad_moorea/don_nolink.recode.vcf")) #output from angsd

##some visualizations
# glPlot(gl, posi="topleft")
##white blocks would indicate missing data

##allele frequency spectrum
# myFreq <- glMean(gl)
# hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
#      main="Distribution of (second) allele frequencies")
# temp <- density(myFreq)
# lines(temp$x, temp$y*1.8,lwd=3)

# assign populations
pops=read.table("part3_inds2pops_no7.txt",sep="\t") #bams file with a 2nd column describing variable of interest
#pops=read.table("part3_bams_no7_year.txt",sep="\t") #another variable I looked at
pop(gl)=pops$V2 #you can have other columns with other variables in other columns, select which one for analysis here
pca=glPca(gl,nf=3,parallel=F) #make pca

#nice pca
quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

##I didn't do this, but useful:
## manually identifying outliers: click on outlier points in the plot.
## adjust n=3 parameter to identify more than 3 points
#outliers=identify(pca$scores[,1:2],labels=gl@ind.names,n=3,cex=0.8)

## re-making PCA without outliers
#gl2=gl[-outliers]
#pca2=glPca(gl2,nf=3,parallel=F)
#colors=as.numeric(as.numeric(as.factor(levels(pop(gl))))) # or use your own colors for populations
#s.class(pca2$scores[],pop(gl2),col=transp(colors,0.5),cstar=1,cellipse=1,clabel=1,axesell=F,grid=F,cpoint=2)

#### PCA stats ####
site <- pops$V2
site <- sub("TO","TNWO",site) #just had to make all sites have 4 characters so I could look at some other variables
site <- sub("TI","TNWI",site)
zone <- str_sub(site,4,4)
realsite <- str_sub(site,1,3)
scores <- pca$scores
adonis(scores ~ zone, strata=realsite, perm=999,method="euclidean")
#this example helped me figure out how to do this: https://stackoverflow.com/questions/43736314/error-in-g-that-non-conformable-arrays
#also this: https://rdrr.io/rforge/vegan/man/adonis.html#heading-1

# Call:
#   adonis(formula = scores ~ site, method = "manhattan") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# site        5    165.86  33.173  1.2078 0.05295  0.268
# Residuals 108   2966.34  27.466         0.94705       
# Total     113   3132.21                 1.00000       

#### ggplot PCA ####
library(ggplot2)
#install.packages("ggrepel")

p1 <- ggordiplots::gg_ordiplot(ord=scores,groups=site,scaling=0.8,choices=c(1,2),conf=0.7,spiders=FALSE,pt.size=2)
#names(p1)
spid <- p1$df_spiders
#p1$df_spiders
gg <- p1$df_ord
gg$zone <- zone
gg$reef <- realsite
quartz()
ggplot(gg,aes(x=x,y=y,color=zone,shape=reef,fill=zone))+
  geom_point(size=2)+
 # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=zone),geom="polygon",alpha=0.1)+
  xlab('PC1 (4.11%)')+
  ylab('PC2 (4.06%)')+
  theme_classic()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(15,16,17),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_linetype_manual(values=c("solid","twodash"),labels=c("Backreef","Forereef"))+
  labs(shape="Site",color="Reef zone",linetype="Reef zone",fill="Reef zone")
  # geom_label(x=-0.2618,y=-0.1116,label="MNW-B",color="black",fill="white")+
  # geom_label(x=-0.3347,y=-0.1352,label="MNW-F",color="black",fill="white")+
  # geom_label(x=-0.1423,y=-0.9173,label="MSE-B",color="black",fill="white")+
  # geom_label(x=0.6595,y=0.8860,label="MSE-F",color="black",fill="white")+
  # geom_label(x=-0.0550,y=-0.1702,label="TNW-B",color="black",fill="white")+
  # geom_label(x=-0.0671,y=0.4204,label="TNW-F",color="black",fill="white")

#now for sequencer instead of sites
years <- read.table("~/Google Drive/Moorea/Host/bamscl_year_no7.txt")
year <- years[,2]
year <- gsub('2015', 'HiSeq 2500',year)
year <- gsub('2017', 'HiSeq 4000',year)

quartz()
ggplot(gg,aes(x=x,y=y,color=year,shape=year,fill=year))+
  geom_point(size=2)+
  # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=year),geom="polygon",alpha=0.1)+
  xlab('PC1 (4.11%)')+
  ylab('PC2 (4.06%)')+
  theme_classic()+
  labs(shape="Sequencer",color="Sequencer",linetype="Sequencer",fill="Sequencer")

##### MDS plot/CCA from IBS matrix #####
library(vegan)
library(adegenet) # for transp()

#[from Misha's angsd_ibs_pca.R script]
bams=read.table("part3_bams_no7")[,1] # list of bam files
goods=c(1:length(bams))

## reading table of pairs of replicates (tab-delimited) - skip if there are no clones
#clonepairs=read.table("clonepairs.tab",sep="\t")
#repsa= clonepairs[,1]
#repsb= clonepairs[,2]
## removing "b" replicates
#goods=which(!(bams %in% repsb))

#--
# loading individual to population correspondences
i2p=read.table("part3_inds2pops_no7.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
#i2p=read.table("part3_bams_no7_year.txt") #for sequencer instead
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]

# setting up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

ma = as.matrix(read.table("donresult.05.ibsMat"))

#---
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

##skipping because removed clones already
#hc=hclust(as.dist(ma),"ave")
#plot(hc,cex=0.5)  # this shows how similar clones are

#ma=ma[goods,goods]
#dimnames(ma)=list(bams[goods],bams[goods])
#hc=hclust(as.dist(ma),"ave")
#plot(hc,cex=0.7) # without clones

# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site) 
# Call:
#   adonis(formula = ma ~ site, data = conds) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# site        1    0.0615 0.061453   1.052 0.00931  0.171
# Residuals 112    6.5427 0.058417         0.99069       
# Total     113    6.6042                  1.00000    

## eigenvectors
#plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
cmd=pp0 #change to pp for CAP, pp0 for MDS
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

## unscaled, to identify outliers
#plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
#ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
#ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
#identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

#### ggplot MDS plot ####
library(ggplot2)
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
library("ggordiplots")
library("stringr")
library("cowplot")
library("extrafont")
loadfonts()

mds <- pp0[["CA"]][["u"]]

p1 <- ggordiplots::gg_ordiplot(ord=mds,groups=site,scaling=0.8,choices=c(1,2),conf=0.7,spiders=FALSE,pt.size=2)
#names(p1)
#spid <- p1$df_spiders
#p1$df_spiders
gg <- p1$df_ord

site <- sub("TO","TNW-F",site) #just had to make all sites have 4 characters so I could look at some other variables
site <- sub("TI","TNW-B",site)
site <- sub("I","-B",site)
site <- sub("O","-F",site)

zone <- str_sub(site,5,5)
realsite <- str_sub(site,1,3)

gg$zone <- zone
gg$reef <- realsite
quartz()
ggplot(gg,aes(x=x,y=y,color=zone,shape=reef,fill=zone))+
  geom_point(size=2)+
  # geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,aes(lty=zone),geom="polygon",alpha=0.1)+
  xlab('MDS1 (5.85%)')+
  ylab('MDS2 (5.76%)')+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(15,16,17),labels=c("MNW","MSE","TNW"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_linetype_manual(values=c("solid","twodash"),labels=c("Backreef","Forereef"))+
  labs(shape="Site",color="Reef zone",linetype="Reef zone",fill="Reef zone")+
  theme(text=element_text(family="Gill Sans MT"))
  
#### K plot from NGSadmix & beagle file ####

# assembling the input table
dir="~/Google Drive/Moorea/2brad_moorea/" # path to input files
inName="part3_don_k2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="part3_inds2pops_no7.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)

names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))

source("~/Google Drive/Moorea/2brad_moorea/plot_admixture_v5_function.R")
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.1)

## recording cluster affiliations
#cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
#save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))


#### K plot from .vcf ####

# primitive look at admixture data:
tbl=read.table("~/Google Drive/Moorea/2brad_moorea/part3_done.2.Q")
barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)

#---
# prettier:

# assembling the input table
dir="~/Google Drive/Moorea/2brad_moorea/" # path to input files
inName="part3_done.2.Q" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
npops=2
pops="part3_inds2pops_no7.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look

source("~/Google Drive/Moorea/2brad_moorea/plot_admixture_v4_function.R")

# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
#tbl$pop=factor(tbl$pop,levels=c("O","K"))

quartz()
ords=plotAdmixture(data=tbl,npops=npops,angle=0,vshift=0,hshift=0)

#### FST from .vcf file ####
#good webiste on calculating this stuff (by hand though):
#http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

library(vcfR)
#install.packages("hierfstat")
library(hierfstat)
library(adegenet)

vcf <- read.vcfR("~/Google Drive/Moorea/2brad_moorea/part3_donresult.vcf")
genind <- vcfR2genind(vcf)
pop(genind) <- i2p$pop

basic.stats(genind)
# $overall
# Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.3114 0.3383 0.3386 0.0003 0.3386 0.0003 0.0007 0.0009 0.0795 0.0005 

fstat(genind, fstonly = FALSE, pop=NULL) #pop=null means you inherit the pops already there 
pairwise.fst(genind, pop = NULL, res.type = c("dist", "matrix"))

test.between(genind,test.lev="Locality",rand.unit="Patch")

# 1          2          3          4          5
# 2 0.01955097                                            
# 3 0.01797238 0.01623331                                 
# 4 0.01732276 0.01493712 0.01528038                      
# 5 0.01559150 0.01361887 0.01387200 0.01315708           
# 6 0.01580877 0.01421282 0.01412182 0.01329455 0.01175710