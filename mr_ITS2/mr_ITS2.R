#Nicola's Moorea ITS2 analysis
#Almost entirely based on DADA2 ITS Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/ITS_workflow.html
#with edits by Carly D. Kenkel and modifications for my data by Nicola Kriefall
#7/23/19

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

##in Terminal home directory:
##following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
##1. download BBMap package, sftp to installation directory
##2. untar: 
#tar -xvzf BBMap_(version).tar.gz
##3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

## my adaptors, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

##Note: Illumina should have cut these out already, normal if you don't get any

##primers for ITS:
# >forward
# GTGAATTGCAGAACTCCGTG
# >reverse
# CCTCCGCTTACTTATATGCTT

##Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list

##cuts off the extra words in the .fastq files
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##only keeping reads that start with the primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq k=15 restrictleft=21 literal=GTGAATTGCAGAACTCCGTG,CCTCCGCTTACTTATATGCTT outm1=${file}_R1_NoIll_ITS.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_ITS.fastq outu2=${file}_R2_check.fastq; done &>bbduk_ITS.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove adapters & reads with Ns in them
#module load cutadapt 
# for file in $(cat samples.list)
# do
# cutadapt -g GTGAATTGCAGAACTCCGTG -G CCTCCGCTTACTTATATGCTT -o ${file}_R1.fastq -p ${file}_R2.fastq --max-n 0 ${file}_R1_NoIll_ITS.fastq ${file}_R2_NoIll_ITS.fastq
# done &> clip.log
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

##did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2)
#packageVersion("dada2")
#I have version 1.10.1 - tutorial says 1.8 but I think that's OK, can't find a version 1.10 walkthrough
library(ShortRead)
#packageVersion("ShortRead")
library(Biostrings)
#packageVersion("Biostrings")

path <- "/Users/nicolakriefall/Desktop/its2/files_rd2" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

#### check for primers ####
FWD <- "GTGAATTGCAGAACTCCGTG"  ## CHANGE ME to your forward primer sequence
REV <- "CCTCCGCTTACTTATATGCTT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
#no primers - amazing

###### Visualizing raw data

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,2,3,4)])
plotQualityProfile(fnFs[c(90,91,92,93)])

#Then look at quality profile of R2 reads

plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(90,91,92,93)])

#starts to drop off around 220 bp for forwards and 200 for reverse, will make approximately that the cutoff values for filter&trim below

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse, added "trimleft" to cut off primers [20 for forward, 21 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(220,200), #leaves overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     minLen = 50,
                     #trimLeft=c(20,21), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce mtching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #if necessary, increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

#~############################~#
##### Dereplicate reads ########
#~############################~#
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaRs[[1]]

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#mostly at 300 bp, which is just what was expected [271-303 bp]

plot(table(nchar(getSequences(seqtab))))

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

#if I wanted to remove some lengths - not recommended by dada2 for its2 data
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(268,327)] 

table(nchar(getSequences(seqtab)))
dim(seqtab)
plot(table(nchar(getSequences(seqtab))))

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 38 bimeras out of 156 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#0.9989
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="its2_reads.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta",tryRC=TRUE,minBoot=70,verbose=TRUE)
unname(head(taxa))

#to come back to later
saveRDS(seqtab.nochim, file="~/Desktop/its2/mrits2_seqtab.nochim.rds")
saveRDS(taxa, file="~/Desktop/its2/mrits2_taxa.rds")

write.csv(seqtab.nochim, file="mrits2_seqtab.nochim.csv")
write.csv(taxa, file="~/Desktop/its2/mrits2_taxa.csv")

#### Reading in prior data files ####
setwd("~/moorea_holobiont/mr_ITS2")
seqtab.nochim <- readRDS("mrits2_seqtab.nochim.rds")
taxa <- readRDS("mrits2_taxa.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')

#import dataframe holding sample information
samdf<-read.csv("mrits_sampledata.csv")
head(samdf)
rownames(samdf) <- samdf$Sample

# Construct phyloseq object (straightforward from dada2 outputs)
# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                sample_data(samdf), 
#                tax_table(taxa))
# 
# ps

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences
path='~/moorea_holobiont/mr_ITS2/mrits2.fasta'
uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)
colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa, rownames(taxa)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))

ps

#removing sample 87 from ps object, no data
ps_no87 <- subset_samples(ps, Sample != "87")
ps_no87

#removing 87 from other things
samdf.no87 <- samdf[-92,]
seqtab.no87 <- seqtab.nochim[-92,]

#### Alpha diversity #####

#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

plot_richness(ps_no87, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

df.div <- data.frame(estimate_richness(ps_no87, split=TRUE, measures =c("Shannon","InvSimpson")))
df.div

df.div$site <- samdf.no87$site
df.div$zone <- samdf.no87$zone
df.div$site_zone <- samdf.no87$site_zone

#write.csv(df.div,file="~/moorea_holobiont/mr_ITS2/mrits_diversity.csv") #saving
#df.div <- read.csv("~/Desktop/mrits/mrits_diversity.csv") #reading back in 

##checking out how diversity works with variables - nothing interesting
#df.size <- read.csv("mr_size.csv",header=TRUE)
#df.div$coral_id <- row.names(df.div)
# 
# mergeddf <- merge(df.div, df.size, by="coral_id", sort=FALSE)
# plot(Shannon~size,data=mergeddf)
# shapiro.test(mergeddf$Shannon)
# shapiro.test(mergeddf$size)
# kruskal.test(Shannon~size,data=mergeddf)
## no significant relationship between coral size & diversity! interesting

## now by sample host heterozygosity 
# df.het <- read.table("host_het.txt") #read in data
# df.het$het <- df.het$V3/(df.het$V2+df.het$V3) #heterozygosity calculations
# df.het$V1 <- sub("TO","TNWO",df.het$V1) #just renaming some sites
# df.het$V1 <- sub("TI","TNWI",df.het$V1) #just renaming some sites
# df.het$site <- substr(df.het$V1, 0, 4)
# df.het$site <- as.factor(df.het$site)
# df.het$site <- sub("I","-B", df.het$site)
# df.het$site <- sub("O","-F", df.het$site)
# str(df.het)

## just getting sample names on the same page
# df.het$V1 <- sub(".trim.bt2.bam.out.saf.idx.ml","",df.het$V1)
# df.het$coral_id <- substr(df.het$V1, 6, 8)
# df.div$coral_id <- row.names(df.div)
# 
# mergeddf2 <- merge(df.div, df.het, by="coral_id", sort=TRUE)
# plot(Shannon~het,data=mergeddf2)
# kruskal.test(Shannon~het,data=mergeddf2)
## not significant!

library(ggplot2)
#install.packages("extrafontdb")
library(extrafontdb)
library(extrafont)
library(Rmisc)
library(cowplot)
library(ggpubr)
#font_import(paths = "/Library/Fonts/")
loadfonts()

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("zone","site"))
diver.si <- summarySE(data=df.div,measurevar=c("InvSimpson"),groupvars=c("zone","site"))

quartz()
gg.sh <- ggplot(diver.sh, aes(x=site, y=Shannon,color=zone,fill=zone,shape=zone))+
  geom_errorbar(aes(ymin=Shannon-se,ymax=Shannon+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15))+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"))+
  guides(fill=guide_legend(title="Reef zone"),color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.sh

gg.si <- ggplot(diver.si, aes(x=site, y=InvSimpson,color=zone,fill=zone,shape=zone))+
  geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Inverse Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15))+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"))+
  guides(fill=guide_legend(title="Reef zone"),color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.si

quartz()
ggarrange(gg.sh,gg.si,common.legend=TRUE,legend="right")

mnw <- subset(df.div,site=="MNW")
mse <- subset(df.div,site=="MSE")
tah <- subset(df.div,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#not different
wilcox.test(InvSimpson~zone,data=mnw)
#p = 0.01

wilcox.test(Shannon~zone,data=mse)
#p = 0.001
wilcox.test(InvSimpson~zone,data=mse)
#p = 0.01

wilcox.test(Shannon~zone,data=tah)
#p = 0.04
wilcox.test(InvSimpson~zone,data=tah)
#p = 0.019

#### *move below ####

# #renaming samples so i can see which pop they're from in bar plot
# sam$newname <- paste(sam$site_zone,sam$Sample,sep="_")
# row.names(sam) <- sam$newname
# row.names(seq2) <- sam$newname
# #row.names(taxa) <- ids
# 
# ps.newname <- phyloseq(otu_table(seq2, taxa_are_rows=FALSE), 
#                     sample_data(sam), 
#                     tax_table(taxa))
# ps.newname

#### Bar plot - raw table ####
top <- names(sort(taxa_sums(ps_no87), decreasing=TRUE))[1:30]
ps.top <- transform_sample_counts(ps_no87, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
plot_bar(ps.top, x="Sample",fill="Class") + facet_wrap(~zone, scales="free_x")

bot <- names(sort(taxa_sums(ps_no87), decreasing=TRUE))[5:118]
ps.bot <- transform_sample_counts(ps_no87, function(OTU) OTU/sum(OTU))
ps.bot <- prune_taxa(bot, ps.bot)
plot_bar(ps.bot, x="Sample",fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#bar plot but without the lines between OTUs & no "NAs"
ps_glom <- tax_glom(ps_no87, "Class")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Class")

#~##########################################~#
###### Apply LULU to cluster ASVs ############
#~##########################################~#

##Necessary pre-steps after producing ASV table and associated fasta file:
##Produce a match list using BLASTn

##IN TERMINAL#

##First produce a blastdatabase with the OTUs
#module load blast+
#makeblastdb -in mrits2.fasta -parse_seqids -dbtype nucl

##Then blast the OTUs against the database to produce the match list 
#blastn -db mrits2.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query mrits2.fasta
##HSP = high scoring pair 
##perc_identity = percent of nucleotides in the highly similar pairings that match

##transfer match_list.txt to R working directory
##BACK IN R#

#first, read in ASV table
#install.packages("remotes")
library("remotes")
#install_github("https://github.com/tobiasgf/lulu.git")
library("lulu")

#must setting up table for lulu, saving it as a csv then reading back in 
write.csv(seqtab.no87,"seqtab.no87.csv")
seq.lulu <- read.csv("seqtab.no87.csv")
rownames(seq.lulu)<-seq.lulu$X
ASVs<-data.frame(t(seq.lulu[,2:119])) 
head(ASVs)

#just made this in terminal
matchList<-read.table("match_list.txt")

#Now, run the LULU curation
curated_result <- lulu(ASVs, matchList,minimum_match=99)
#default: minimum_relative_cooccurence = 0.95 default, changed to 0.70 to see what happens, nothing
#default: minimum_match = 84 default, only 1 OTU different between 97 & 84
summary(curated_result)
#me: leaves 7 OTUs with 84% match - take home story is that they're all very similar
#19 otus with match of 99$

#Pull out the curated OTU list, re-transpose
lulu.out <- data.frame(t(curated_result$curated_table))

#Continue on to your favorite analysis
write.csv(lulu.out,"~/moorea_holobiont/mr_ITS2/lulu_output.csv")

#removing X from row names to match up with previous data frames
rownames(lulu.out) <- sub("X","",rownames(lulu.out))

#phyloseq object with lulu results
ps.lulu <- phyloseq(otu_table(lulu.out, taxa_are_rows=FALSE), 
                    sample_data(samdf.no87), 
                    tax_table(taxa2))

#### Bar plot & alpha diversity - lulu results ####
ps.top <- transform_sample_counts(ps.lulu, function(OTU) OTU/sum(OTU))
plot_bar(ps.top, x="Sample",fill="Class") + facet_wrap(~zone, scales="free_x")

plot_richness(ps.lulu, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

df.div <- data.frame(estimate_richness(ps.lulu, split=TRUE, measures =c("Shannon","InvSimpson")))
df.div

df.div$site <- samdf.no87$site
df.div$zone <- samdf.no87$zone
df.div$site_zone <- samdf.no87$site_zone

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("zone","site"))

ggplot(diver.sh, aes(x=site,y=Shannon,color=zone,shape=zone))+
  geom_errorbar(aes(ymin=Shannon-se,ymax=Shannon+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))

kruskal.test(Shannon~zone,data=df.div)

mnw <- subset(df.div,site=="MNW")
mse <- subset(df.div,site=="MSE")
tah <- subset(df.div,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#not different
wilcox.test(InvSimpson~zone,data=mnw)
#p = 0.01

wilcox.test(Shannon~zone,data=mse)
#p = 0.001
wilcox.test(InvSimpson~zone,data=mse)
#p = 0.01

wilcox.test(Shannon~zone,data=tah)
#p = 0.04
wilcox.test(InvSimpson~zone,data=tah)
#p = 0.019

#### MCMC.OTU to remove underrepresented ASVs ####
library(MCMC.OTU)

#added a column with a blank name in the beginning, with 1-95 in the column, mcmc.otu likes this
#also removed the X from the beginning of sample names
lulu.mcmc <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_formcmc.csv")

goods <- purgeOutliers(lulu.mcmc,count.columns=3:21,otu.cut=0.001,zero.cut=0.02)
#otu.cut = 0.1% of reads represented by ASV 
colnames(goods)
#sq 1, 3, 6, 7 retained with min 84% matching in lulu
#sq 1, 2, 3, 6, 7, 12, 18, 24, 32, 33 with min 99% matching in lulu

#### Bar plot - clustered & trimmed ####
rownames(goods) <- goods$sample
counts <- goods[,3:11]

ps.lulu <- phyloseq(otu_table(counts, taxa_are_rows=FALSE), 
                    sample_data(samdf.no87), 
                    tax_table(taxa2))

ps_glom <- tax_glom(ps.lulu, "Class")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Class")

#### *move to norm section ####
#just checking out normalized stuff 
plot_richness(ps.norm, x="site", measures=c("Shannon", "Simpson"), color="in_off") + theme_bw()

df.norm <- data.frame(estimate_richness(ps.norm, split=TRUE, measures =c("Shannon","InvSimpson")))
df.norm$site <- sam$site
df.norm$zone <- sam$in_off

ggplot(df.norm,aes(x=site,y=Shannon,color=zone))+
  geom_boxplot()

#checking if my variables have sig different sample reads before normalizing
sam$new <- sample_sums(ps_no87)
plot(new~site_zone,data=sam) #MSEO has way more than everyone else
shapiro.test(sam$new)
a1 <- aov(new~site_zone,data=sam)
summary(a1)
plot(new~island,data=sam) #even
plot(new~site,data=sam) #MSE still more
plot(new~in_off,data=sam) #outer more
a1 <- aov(new~in_off,data=sam)
summary(a1) # p = 0.0587 lol

#### DESEQ to normalize counts ####
#BiocManager::install("DESeq2")
library("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("DEFormats")
library(edgeR)
library(DEFormats)

#count data
seq <- read.csv("mrits2_seqtab.nochim_no87.csv")
#getting rid of crap ones
#mcmc.otu requires a blank column with sample names first, did manually in excel

seq <- purgeOutliers(seq,count.columns=4:120,sampleZcut=-2.5,otu.cut=0.001)
#bad samples: 513, 530, 76, 505
#49 OTUs pass frequency cutoff

row.names(seq) <- seq$sample #setting my sample names as row names so I can remove the sample name column
seq2 <- seq[,3:52] #now removing sample name columns
seqt <- t(seq2) #deseq requires count data in columns

#sample data
sam <- read.csv("mrits_sampledata_no87.csv")
row.names.remove <- c("513","530","76","505")
sam <- sam[!(row.names(sam) %in% row.names.remove), ]
row.names(sam) <- sam$Sample

#method from here:
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex#dereplication
deseq_counts <- DESeqDataSetFromMatrix(seqt, colData = sam, design = ~site*zone) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. 

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
#vst_trans_count_tab[vst_trans_count_tab < 0.0] <- 0.0 #getting rid of 0s

# back into phyloseq
vstt <- t(vst_trans_count_tab)
ps.norm <- phyloseq(otu_table(vstt, taxa_are_rows=FALSE), 
                    sample_data(sam), 
                    tax_table(taxa))
ps.norm

sam$new.norm <- sample_sums(ps.norm)
plot(new.norm~site_zone,data=sam) #more normalized now
a1 <- aov(new.norm~site_zone,data=sam)
summary(a1) #p = 0.0646 lol
plot(new.norm~island,data=sam) #even
plot(new.norm~site,data=sam) #MSE still a little more
plot(new.norm~in_off,data=sam) #outer more
a1 <- aov(new.norm~in_off,data=sam)
summary(a1) # p = significant

#### Pcoa #####
write.csv(vstt,"~/moorea_holobiont/mr_ITS2/vstt.csv")

#before normalizing
ps.ord <- ordinate(ps_no87, method="MDS", distance="euclidean")
ps.ord.eigs <- ps.ord$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
plot_ordination(ps_no87, ps.ord, color="in_off") + 
  geom_point(size=1) #+ 

#after normalizing
vst_pcoa <- ordinate(ps.norm, method="MDS", distance="euclidean")
#vst_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#^this looks horrible
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
quartz()
plot_ordination(ps.norm, vst_pcoa, color="zone",shape="zone") + 
  geom_point(size=1) +
  facet_wrap(~site)+
  stat_ellipse(level=0.95)+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(color="Reef zone",shape="Reef zone")
  #geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  #coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")
  #scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
  #theme(legend.position="none")

#stats
library(dendextend)
# and calculating our Euclidean distance matrix
euc_dist <- dist(vstt)
euc_clust <- hclust(euc_dist, method="ward.D2")

anova(betadisper(euc_dist, sam$site_zone)) #not sig
#significant means I can't do straight adonis
#wasn't significnt for inshore/offshore
adonis(euc_dist~sam$site*sam$zone) # all sig

#by site
adonis( ~ zone, data = , method="manhattan")


#install.packages('dendextend')
#### hierarchical clustering of samples ####
library(dendextend)
# and calculating our Euclidean distance matrix
euc_dist <- dist(vstt)
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)

sam$color <- sam$site_zone
sam$color <- gsub("MNWO","darkred",sam$color)
sam$color <- gsub("MNWI","red",sam$color)
sam$color <- gsub("MSEI","yellow",sam$color)
sam$color <- gsub("MSEO","darkgoldenrod1",sam$color)
sam$color <- gsub("TI","green",sam$color)
sam$color <- gsub("TO","darkgreen",sam$color)
dend_cols <- as.character(sam$color[order.dendrogram(euc_dend)])

labels_colors(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")
#no clustering by sites

#now by in/off
sam$color2 <- sam$in_off
sam$color2 <- gsub("Inner","orange",sam$color2)
sam$color2 <- gsub("Outer","purple",sam$color2)
dend_cols2 <- as.character(sam$color2[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols2
plot(euc_dend, ylab="VST Euc. dist.")

#trying to select only the variables I want
#resource for analysis here:
#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

#ps.mnw = subset_samples(ps_no87, site=="MNW")
#ps.mnw = subset_samples(ps.norm, site=="MNW")
ps.mnw = subset_samples(ps.rare, site=="MNW")
ds.mnw = phyloseq_to_deseq2(ps.mnw, ~ zone)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
#nothing for MNW, even at 0.1 prior to normalized & after

#ps.mse = subset_samples(ps.norm, site=="MSE")
ps.mse = subset_samples(ps.rare, site=="MSE")
ds.mse = phyloseq_to_deseq2(ps.mse, ~ zone)
dds.mse <- estimateSizeFactors(ds.mse,type="poscounts")
stat.mse = DESeq(dds.mse, test="Wald", fitType="parametric")

res = results(stat.mse, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.mse)[rownames(sigtab), ], "matrix"))
dim(sigtab)
#5 for MSE - not normalized
#2 for MSE normalized at 0.1, none for 0.05
#1 for MSE after 

ps.t = subset_samples(ps.rare, site=="TNW")
#ps.t = subset_samples(ps.norm, site=="TNW")
ds.t = phyloseq_to_deseq2(ps.t, ~ in_off)
dds.t <- estimateSizeFactors(ds.t,type="poscounts")
stat.t = DESeq(dds.t, test="Wald", fitType="parametric")

res = results(stat.t, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.t)[rownames(sigtab), ], "matrix"))
dim(sigtab)
#1 for tahiti at 0.05 & 0.1 lol
#none for tahiti after normalized

# #visualizing deseq results - can't do with no results
# theme_set(theme_bw())
# scale_fill_discrete <- function(palname = "Set1", ...) {
#   scale_fill_brewer(palette = palname, ...)
# }
# # Phylum order
# x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# # Class order
# x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#### stats ####
library(vegan)
#assumption of homogeneity of dispersion within groups must be tested
anova(betadisper(euc_dist, sam$in_off)) #not sig
#significant means I can't do straight adonis
#wasn't significnt for inshore/offshore
adonis(euc_dist~sam$site*sam$in_off) # all sig

#### Rarefy ####
library(vegan)

rare <- seq2 #creating a copy so I can mess with it
rarecurve(rare, step = 100, label=TRUE)
#going to remove samples < 2500 reads

rare$total <- rowSums(rare)
seq.less <- subset(rare,total >= 2500)
#lose a lot of samples here :(

rarecurve(seq.less, step = 100, label=FALSE)
rowSums(seq.less)

seq.rare <- rrarefy(seq.less,sample=2500)
rarecurve(seq.rare,step=100,label=FALSE)

write.csv(seq.rare,"~/moorea_holobiont/mr_ITS2/seq.rare_1.csv")

ps.rare <- phyloseq(otu_table(seq.rare, taxa_are_rows=FALSE), 
               sample_data(sam), 
               tax_table(taxa))

ps.rare

plot_richness(ps.rare, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

df.rare <- data.frame(estimate_richness(ps.rare, split=TRUE, measures =c("Shannon","InvSimpson")))

sam$total <- rare$total
sam.rare <- subset(sam,total >= 2500)

df.rare$site <- sam.rare$site
df.rare$zone <- sam.rare$zone

ggplot(df.rare,aes(x=site,y=InvSimpson,color=zone))+
  geom_boxplot()

#stats:

wilcox.test(Shannon~zone,data=df.rare)

mnw <- subset(df.rare,site=="MNW")
mse <- subset(df.rare,site=="MSE")
tah <- subset(df.rare,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#not different
wilcox.test(InvSimpson~zone,data=mnw)
#p = 0.06
wilcox.test(Shannon~zone,data=mse)
#very different p-value = 0.005113
wilcox.test(InvSimpson~zone,data=mse)
#p = 0.04
wilcox.test(Shannon~zone,data=tah)
#p = 0.03
wilcox.test(InvSimpson~zone,data=tah)
#p = 0.03

####Bar plot for rarefied results#####
#renaming samples so i can see which pop they are on the x axis
sam.rare$newname <- paste(sam.rare$site_zone,sam.rare$Sample,sep="_")
row.names(sam.rare) <- sam.rare$newname
row.names(seq.rare) <- sam.rare$newname
#row.names(taxa) <- ids

ps.rare2 <- phyloseq(otu_table(seq.rare, taxa_are_rows=FALSE), 
                    sample_data(sam.rare), 
                    tax_table(taxa))

top3 <- names(sort(taxa_sums(ps.rare2), decreasing=TRUE))[1:118]
ps.top3 <- transform_sample_counts(ps.rare2, function(OTU) OTU/sum(OTU))
ps.top3 <- prune_taxa(top3, ps.top3)
quartz()
plot_bar(ps.top3, x="newname",fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

bot30 <- names(sort(taxa_sums(ps.rare2), decreasing=FALSE))[1:100]
ps.bot30 <- transform_sample_counts(ps.rare2, function(OTU) OTU/sum(OTU))
ps.bot30 <- prune_taxa(bot30, ps.bot30)
plot_bar(ps.bot30, x="newname", fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")


#~############################~#
##### output 'OTU' table #######
#~############################~#

#seqtab.nochim is the 'OTU' table...but is a little unwieldy
#For Symbiodinium, sequence classification is not so great...
#want fasta file of 'OTUs' and table designated by 'OTU'

#First, output fasta file for 'OTUs'
path='~/Desktop/mrits/mrits.fasta'
uniquesToFasta(seqtab.nochim, path, ids = NULL, mode = "w", width = 20000)
#the headings don't match up:
#install.packages("phylotools")
library("phylotools")
name_table <- read.table("~/Google Drive/Moorea/ITS2/rename_fasta.txt")
rename.fasta(infile="~/Google Drive/Moorea/ITS2/mrits.fasta",ref_table =name_table,outfile="~/Google Drive/Moorea/ITS2/mrits_renamed.fasta")

#then, rename output table and write it out
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim)<-ids

write.csv(seqtab.nochim,file="mrits_seqtab_newids.csv",quote=F)
##replace sequences with shorter names (correspondence table output below)
taxa_names(ps)<-ids
str(seqtab.nochim)


#trimming low abundance ones using mcmc.otu
library("MCMC.OTU")

#had to add a first column, blank title with just 1-96 listed
lulu.mc <- read.csv("~/moorea_holobiont/mr_ITS2/lulu.mcmc.csv")

goods2=purgeOutliers(lulu.mc,count.columns=8:15,otu.cut=0.001,zero.cut=0.03)
goods2
#after lulu & mcmc.otu, only keeps 1, 3, 6, 7, 12
less.taxa <- taxa[c(1,3,6,7,12),]
row.names(less.taxa) <- c("sq1","sq3","sq6","sq7","sq12")

counts.less <- goods2[,c(1,8:12)]
data.less <- goods2[,1:7]
row.names(counts.less) <- goods2$sample
row.names(data.less) <- goods2$sample

ps.wayless <- phyloseq(otu_table(counts.less, taxa_are_rows=FALSE), 
                             sample_data(data.less), 
                             tax_table(less.taxa))

plot_bar(ps.all,x="site_zone",y="Abundance",fill="Class")

library(dplyr)
ps.all <- transform_sample_counts(ps.wayless, function(OTU) OTU/sum(OTU))
tb <- psmelt(ps.all)%>%
  filter(!is.na(Abundance))%>%
  group_by(site_zone,OTU)%>%
  summarize_at("Abundance",mean)

tb$sym <- gsub("sq12","Clade C",tb$OTU)
tb$sym <- gsub("sq1","C3k.1",tb$sym)
tb$sym <- gsub("sq3","C3k.2",tb$sym)
tb$sym <- gsub("sq6","C3k.3",tb$sym)
tb$sym <- gsub("sq7","A1",tb$sym)

tb$site_zone <- gsub("MNWI","MNW-B",tb$site_zone)
tb$site_zone <- gsub("MNWO","MNW-F",tb$site_zone)
tb$site_zone <- gsub("MSEI","MSE-B",tb$site_zone)
tb$site_zone <- gsub("MSEO","MSE-F",tb$site_zone)
tb$site_zone <- gsub("TI","TNW-B",tb$site_zone)
tb$site_zone <- gsub("TO","TNW-F",tb$site_zone)

#library(cowplot)
tb$sym <- factor(tb$sym, levels=c("C3k.1","C3k.2","C3k.3","Clade C","A1"))
quartz()
ggplot(tb,aes(x=site_zone,y=Abundance,fill=sym))+
  geom_bar(stat="identity", colour="black")+
  theme_cowplot()+
  theme(text=element_text(family="Gill Sans MT"))+
  xlab('Site')+
  scale_fill_manual(name="Sym.",values=c("seagreen2","seagreen3","seagreen4","green","blue"))
# Kingdom        Phylum     Class 
# sq7  "Symbiodinium" " Clade A" " A1" 
# sq1  "Symbiodinium" " Clade C" " C3k"
# sq12 "Symbiodinium" " Clade C" NA    
# sq3  "Symbiodinium" " Clade C" " C3k"
# sq6  "Symbiodinium" " Clade C" " C3k"

#now redoing deseq results with lulu output
row.names(taxa) <- ids
counts <- alldat[,7:15]
#removing sample 87 - no data
counts.no87 <- counts[c(1:91,93:96),]

sams <- alldat[,1:6]
sams.no87 <- sams[c(1:91,93:96),]

ps.lulu <- phyloseq(otu_table(counts.no87, taxa_are_rows=FALSE), 
                    sample_data(sams.no87), 
                    tax_table(taxa))
ps.lulu

ps.mnw = subset_samples(ps.lulu, site=="MNW")
ds.mnw = phyloseq_to_deseq2(ps.mnw, ~ in_off)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab), ], "matrix"))
head(sigtab)
#nothing for MNW, even at 0.1

ps.mse = subset_samples(ps.lulu, site=="MSE")
ds.mse = phyloseq_to_deseq2(ps.mse, ~ in_off)
dds.mse <- estimateSizeFactors(ds.mse,type="poscounts")
stat.mse = DESeq(dds.mse, test="Wald", fitType="parametric")

res = results(stat.mse, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.mse)[rownames(sigtab), ], "matrix"))
dim(sigtab)
#5 for MSE - all c3k
#1 for MSE with lulu data

ps.t = subset_samples(ps.lulu, site=="T")
ds.t = phyloseq_to_deseq2(ps.t, ~ in_off)
dds.t <- estimateSizeFactors(ds.t,type="poscounts")
stat.t = DESeq(dds.t, test="Wald", fitType="parametric")

res = results(stat.t, cooksCutoff = FALSE)
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.t)[rownames(sigtab), ], "matrix"))
dim(sigtab)
#1 for tahiti at 0.05 & 0.1 lol
#also 1 for tahiti with lulu - A1

#visualizing deseq results
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

####Bar plot for lulu results#####
#renaming samples so i can see which pop they are on the x axis
sams.no87$newname <- paste(sams.no87$site_zone,sams.no87$X,sep="_")
row.names(sams.no87) <- sams.no87$newname
row.names(counts.no87) <- sams.no87$newname
row.names(taxa) <- ids

ps.lulu <- phyloseq(otu_table(counts.no87, taxa_are_rows=FALSE), 
                    sample_data(sams.no87), 
                    tax_table(taxa))
tax_table(ps.lulu)

top3 <- names(sort(taxa_sums(ps.lulu), decreasing=TRUE))[1:9]
ps.top3 <- transform_sample_counts(ps.lulu, function(OTU) OTU/sum(OTU))
ps.top3 <- prune_taxa(top3, ps.top3)
plot_bar(ps.top3, x="newname",fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

bot30 <- names(sort(taxa_sums(ps.lulu), decreasing=FALSE))[1:30]
ps.bot30 <- transform_sample_counts(ps.lulu, function(OTU) OTU/sum(OTU))
ps.bot30 <- prune_taxa(bot30, ps.bot30)
plot_bar(ps.bot30, x="newname", fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

library("vegan")
div <- as.data.frame(diversity(counts.no87, index = "shannon", MARGIN = 1, base = exp(1)))
div$name <- row.names(div)
div$dive <- div$`diversity(counts.no87, index = "shannon", MARGIN = 1, base = exp(1))`
div$site_zone <- sams.no87$site_zone

div.se <- summarySE(data=div,measurevar="dive",groupvars="site_zone")
ggplot(div.se,aes(x=site_zone,y=dive))+
  geom_point()+
  geom_errorbar(aes(ymin=dive-se,ymax=dive+se))

#### DESEQ normalization of lulu counts ####
ds.lulu = phyloseq_to_deseq2(ps.lulu, ~ site_zone)
dds.lulu <- estimateSizeFactors(ds.lulu,type="poscounts")
stat.lulu = DESeq(dds.lulu)
newcounts <- counts(stat.lulu,normalize=TRUE)
counts.luluseq <- t(newcounts)

ps.luluseq <- phyloseq(otu_table(counts.luluseq, taxa_are_rows=FALSE), 
                    sample_data(sams.no87), 
                    tax_table(taxa))

top3 <- names(sort(taxa_sums(ps.luluseq), decreasing=TRUE))[1:9]
ps.top3 <- transform_sample_counts(ps.luluseq, function(OTU) OTU/sum(OTU))
ps.top3 <- prune_taxa(top3, ps.top3)
plot_bar(ps.top3, x="newname",fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

ps.luluseq.melt <- psmelt(ps.luluseq)
ps.luluseq.melt$newclass <- paste(ps.luluseq.melt$Class,ps.luluseq.melt$OTU,sep="_")

ggplot(ps.luluseq.melt,aes(x=site_zone,y=Abundance,color=newclass))+
  geom_boxplot()

c3k <- subset(ps.luluseq.melt,Class==" C3k")
ggplot(c3k,aes(x=site_zone,y=Abundance,color=Class))+
  geom_boxplot()

c3k <- subset(ps.luluseq.melt,Class==" C3k")
ggplot(c3k,aes(x=site_zone,y=Abundance,color=Class))+
  geom_boxplot()

#Ordinate Samples
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="zone", title="Bray NMDS")
#doesn't converge, skipping

#Bar-plots

top30 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:30]
ps.top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top30 <- prune_taxa(top30, ps.top30)
plot_bar(ps.top30, x="site_zone", fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

btm30 <- names(sort(taxa_sums(ps), decreasing=FALSE))[1:30]
ps.btm30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.btm30 <- prune_taxa(btm30, ps.btm30)
plot_bar(ps.btm30, x="site_zone", fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")


#Mess around with other stuff in phyloseq here...

#### Principal coordinate analysis ####

library(vegan)
#install.packages("MCMC.OTU")
library(MCMC.OTU)
#install.packages("ggfortify")
library(ggfortify)
library(cluster)
#install.packages('labdsv')
library(labdsv)

#have to add a first column with no name numbered 1-93
alldat<-read.csv("~/Desktop/its2/seqtab.csv")
names(alldat)
str(alldat)
head(alldat) 
#dat<-alldat[c(1:33,35:91,93:96),] #to remove rows with missing data. mine were #34 [380 TI] & #92 [87 MNWI]

# purging under-sequenced samples; and OTUs represented in less than 3% of all samples (~must occur at least 3 separate times)
goods2=purgeOutliers(alldat,count.columns=3:159,sampleZcut=-2.5,otu.cut=0.001)
# [1] "samples with counts below z-score -2.5 :"
# [1] "513" "530" "76" 
# [1] "zscores:"
# [1] -3.952059 -2.521370 -4.035209
# [1] "OTUs passing frequency cutoff  0.001 : 38"
write.csv(goods2,"~/Desktop/its2/goods2.csv")

# creating a log-transfromed normalized dataset for PCoA:
goods.log=logLin(data=goods2,count.columns=8:length(names(goods2)))

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist=vegdist(goods.log,method="manhattan")
goods.pcoa=pcoa(goods.dist)

scores=goods.pcoa$vectors
conditions=goods2[,1:8]
margin=0.01
quartz()
plot(scores[,1], scores[,2],type="n",
     xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     #mgp=c(2.3,1,0),
     xlab=paste("Axis1 (",round(goods.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
     ylab=paste("Axis2 (",round(goods.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
     main="")
points(scores[conditions$in_off=="Inner",1],scores[conditions$in_off=="Inner",2],pch=1)
points(scores[conditions$in_off=="Outer",1],scores[conditions$in_off=="Outer",2],pch=19)
legend("bottomright", c("Inshore","Offshore"), pch=c(1, 19), cex=0.8)
#super weird separation 

unname(taxa)
#############################################
# exploring correlations between OTUs 

#goods=purgeOutliers(dat,count.columns=8,otu.cut=0.030,zero.cut=0.10) #just to reduce number represented in scatter plot; could also do subsets of interesting majority types

# creating a log-transformed normalized dataset ignoring zero counts:
nl=startedLog(data=goods2,count.columns=8:length(names(goods2)),logstart=0)

# displaying a matrix of scatterplots and p-values of OTU correlations
# (onlu p-values better than 0.1 are displayed)
pairs(nl,lower.panel=panel.smooth,upper.panel=panel.cor.pval) 
#some significantly correlated - will have to explore more


#~############################~#
######### Sig OTUs only ########
#~############################~#

goods2=purgeOutliers(dat,count.columns=8:82,otu.cut=0.001,zero.cut=0.0288)
goods2
# [1] "samples with counts below z-score -2.5 :"
# [1] "178" "513" "58"  "76" 
# [1] "zscores:"
# 15        56        68        91 
# -3.659871 -3.334992 -3.412355 -3.360231 
# [1] "OTUs passing frequency cutoff  0.001 : 34"
# [1] "OTUs with counts in 0.0288 of samples:"
# 
# TRUE 
# 34 

# must get rid of all sqs past 34, and samples 15, 34, 56, 68, 91, 92
seqtab.sub <- seqtab.nochim[c(1:14,16:33,35:55,57:67,69:90,93:96),c(1:34)]
str(seqtab.sub)

# must get rid of samples 15, 34, 56, 68, 91, 92
samdf.sub <- samdf[c(1:14,16:33,35:55,57:67,69:90,93:96),]
row.names(samdf.sub) <- samdf.sub$Sample
str(samdf.sub)

# the sqs again
taxa.sub <- taxa[c(1:34),]

ps2 <- phyloseq(otu_table(seqtab.sub, taxa_are_rows=FALSE), 
                sample_data(samdf.sub), 
                tax_table(taxa.sub))

all <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:34]
ps.all <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps.bar <- prune_taxa(all, ps.all)
plot_bar(ps.bar,x="site_zone", fill="Class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps.bar)
write.csv(psz, file="mrits_phyloseqmelt.csv")

psz2 <- read.csv("mrits_phyloseqmelt.csv")
psz2$Sample <- as.factor(psz2$Sample)
ggplot(psz2, aes(x=Sample, y=Abundance, fill=Class))+
  geom_bar(stat="identity", colour="black") #+ facet_wrap("tank")

library(plyr)
library(dplyr)

class <- psz2$Class
id <- psz2$Sample
abun <- psz2$Abundance
new <- cbind.data.frame(class,id,abun)

daily <- group_by(new, class, id)
per_day <- summarize(daily, sumabun = sum(abun))

ggplot(per_day, aes(x=id, y=sumabun, fill=class))+
  geom_bar(stat="identity", colour="black") #+ facet_wrap("tank")

sz <- psz2$site_zone
class <- psz2$Class
abun <- psz2$Abundance

new2 <- cbind.data.frame(sz,class,abun)
group <- group_by(new2, class, sz)
grouped <- summarize(group, sumabun = sum(abun))
groupednew <- group_by(grouped,sz)
mnwi <- subset(grouped,sz=="MNWI")
mnwo <- subset(grouped,sz=="MNWO")
msei <- subset(grouped,sz=="MSEI")
mseo <- subset(grouped,sz=="MSEO")
ti <- subset(grouped,sz=="TI")
to <- subset(grouped,sz=="TO")

sum(mnwi$sumabun)
13
sum(mnwo$sumabun)
15
sum(msei$sumabun)
16
sum(mseo$sumabun)
16
sum(ti$sumabun)
14
sum(to$sumabun)
16

lut <- c("MNWI"=13, "MNWO"=15, "MSEI"=16, "MSEO"=16,"TI"=14,"TO"=16)
grouped$relabun <- lut[grouped$sz]
grouped$rel <- grouped$sumabun/grouped$relabun

grouped$sz <- factor(grouped$sz,levels=c("MNWI","MSEI","TI","MNWO","MSEO","TO"))

#site & zone
quartz()
ggplot(grouped, aes(x=sz, y=rel, fill=class))+
  geom_bar(stat="identity", colour="black") +
  theme_classic()+
  xlab("Site & Reef Zone")+
  ylab("Relative Abundances")+
  guides(fill=guide_legend(title="Class"))+
  theme(text=element_text(family="Gill Sans MT"))

#just zone
zone <- psz2$in_off
class <- psz2$Class
abun <- psz2$Abundance

new2 <- cbind.data.frame(class,zone,abun)
group <- group_by(new2, class, zone)
grouped <- summarize(group, sumabun = sum(abun))

ins <- subset(grouped,zone=="Inner")
out <- subset(grouped,zone=="Outer")

sum(ins$sumabun)
43
sum(out$sumabun)
47

lut <- c("Inner"=43, "Outer"=47)
grouped$relabun <- lut[grouped$zone]
grouped$rel <- grouped$sumabun/grouped$relabun

quartz()
ggplot(grouped, aes(x=zone, y=rel, fill=class))+
  geom_bar(stat="identity", colour="black") +
  theme_classic()+
  xlab("Reef Zone")+
  ylab("Relative Abundances")+
  guides(fill=guide_legend(title="Class"))+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_x_discrete(labels=c("Inshore","Offshore"))

#just island
isl <- psz2$island
class <- psz2$Class
abun <- psz2$Abundance

isla <- cbind.data.frame(class,isl,abun)
group <- group_by(isla, class, isl)
grouped <- summarize(group, sumabun = sum(abun))

mor <- subset(grouped,isl=="Moorea")
tah <- subset(grouped,isl=="Tahiti")

sum(mor$sumabun)
57.86458
sum(tah$sumabun)
31.24531

lut <- c("Moorea"=57.86458, "Tahiti"=31.24531)
grouped$relabun <- lut[grouped$isl]
grouped$rel <- grouped$sumabun/grouped$relabun

quartz()
ggplot(grouped, aes(x=isl, y=rel, fill=class))+
  geom_bar(stat="identity", colour="black") +
  theme_classic()+
  xlab("Island")+
  ylab("Relative Abundances")+
  guides(fill=guide_legend(title="Class"))+
  theme(text=element_text(family="Gill Sans MT"))

#PCA
library(vegan)
#BiocManager::install("ggordiplots")
#install.packages('devtools')
library(devtools)
#install_github("jfq3/ggordiplots")
library('ggordiplots')

otu <- read.csv("mrits_mcmc.csv")
head(otu)
rownames(otu) <- otu$sample #make sample IDs the rownames
length(names(otu))
dat<-otu[c(1:33,35:91,93:96),] #to remove rows with missing data. mine were #34 [380 TI] & #92 [87 MNWI]

goods=purgeOutliers(dat,count.columns=8:82,otu.cut=0.001,zero.cut=0.0288) #just to reduce number represented in scatter plot; could also do subsets of interesting majority types

goods.log=logLin(data=goods,count.columns=8:length(names(goods)))

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist=vegdist(goods.log,method="manhattan")
goods.pcoa=pcoa(goods.dist)
goods.pcoa
3.116452e-01
2.079159e-01

scores=goods.pcoa$vectors
conditions=goods[,1:7]
conditions
zone <- conditions$in_off
zone
p1 <- ggordiplots::gg_ordiplot(ord=goods.pcoa$vectors,groups=zone,choices=c(1,2),conf=0.95,spiders=TRUE,pt.size=2)
names(p1)
spid <- p1$df_spiders
p1$df_spiders
gg <- p1$df_ord
quartz()
ggplot(gg,aes(x=x,y=y,group=Group,color=Group,shape=Group,fill=Group))+
  geom_point(size=2)+
  geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.95)+
  xlab('Axis 1 (31.2%)')+
  ylab('Axis 2 (20.8%)')+
  theme_classic()+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_shape_manual(values=c(21,24),labels=c("Inshore","Offshore"))+
  scale_fill_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  scale_color_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  labs(shape="Reef zone",color="Reef zone",fill="Reef zone")

adonis(goods.dist ~ zone, data = goods.log, method="manhattan")

# scorez <- as.data.frame(scores)
# scorez$zone <- conditions$in_off
# adonis(Axis.1~zone, data=scorez, permutations = 999, method = "manhattan", strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", parallel = getOption("mc.cores"))
# adonis(formula, data, permutations = 999, method = "bray", strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)