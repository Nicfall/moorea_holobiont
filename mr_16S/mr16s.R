#Nicola's Moorea 16S analysis
#Almost entirely based on DADA2 Pipeline 1.8 Walkthrough:
#https://benjjneb.github.io/dada2/tutorial.html
#with edits by Carly D. Kenkel and modifications for my data by Nicola Kriefall
#12/30/18

#~########################~#
##### PRE-PROCESSING #######
#~########################~#

#fastq files should have R1 & R2 designations for PE reads
#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#in Terminal home directory:
#following instructions of installing BBtools from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/
#1. download BBMap package, sftp to installation directory
#2. untar: 
#tar -xvzf BBMap_(version).tar.gz
#3. test package:
#cd bbmap
#~/bin/bbmap/stats.sh in=~/bin/bbmap/resources/phix174_ill.ref.fa.gz

# my adaptors for 16S, which I saved as "adaptors.fasta"
# >forward
# AATGATACGGCGACCAC
# >forwardrc
# GTGGTCGCCGTATCATT
# >reverse
# CAAGCAGAAGACGGCATAC
# >reverserc
# GTATGCCGTCTTCTGCTTG

#primers for 16S: 
# >forward
# GTGYCAGCMGCCGCGGTA
# >reverse
# GGACTACHVGGGTWTCTAAT

##Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list

##cuts off the extra words in the .fastq files
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 

##gets rid of reads that still have the adaptor sequence, shouldn't be there, I didn't have any
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log

##getting rid of first 4 bases (degenerate primers created them)
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log

##only keeping reads that start with the 16S primer
#for file in $(cat samples.list); do ~/bin/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=20 k=10 literal=GTGYCAGCMGCCGCGGTA,GGACTACHVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log
##higher k = more reads removed, but can't surpass k=20 or 21

##using cutadapt to remove primer
# for file in $(cat samples.list)
# do
# cutadapt -g GTGYCAGCMGCCGCGGTA -a ATTAGAWACCCVHGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1.fastq -p ${file}_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq
# done &> clip.log
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed
##this overwrote my original renamed files 

# did sftp of *_R1.fastq & *_R2.fastq files to the folder to be used in dada2

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2); packageVersion("dada2")
#I have version 1.10.0 - tutorial says 1.8 but I think that's OK, can't find a version 1.10 walkthrough
library(ShortRead)
#packageVersion("ShortRead")
library(Biostrings)
#packageVersion("Biostrings")
path <- "~/Desktop/mr16s/files" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))
#no primers - amazing

#### Visualizing raw data ####

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])
plotQualityProfile(fnFs.filtN[c(90,91,92,93)])
#looks mostly good up to 200

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])
plotQualityProfile(fnRs.filtN[c(90,91,92,93)])
#180

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
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

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(244,264)] #again, being fairly conservative wrt length

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 105 bimeras out of 2485 input sequences.

sum(seqtab.nochim)/sum(seqtab)
#0.9959795
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="mr16s_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

# #Using package DECIPHER as an alternatie to 'assignTaxonomy'
# 
# #BiocManager::install("DECIPHER")
# library(DECIPHER); packageVersion("DECIPHER")
# #citation("DECIPHER")
# 
# #http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified) file to follow along.
# 
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/Downloads/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE, threshold=50) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#also doing other taxonomy method:
#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr_v132_train_set.fa.gz",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v132.fa.gz",tryRC=TRUE,verbose=TRUE)
# 237 out of 2380 were assigned to the species level.
# Of which 214 had genera consistent with the input table.

saveRDS(taxa.plus, file="mr16s_taxaplus.rds")
saveRDS(taxa, file="mr16s_taxa.rds")
write.csv(taxa.plus, file="mr16s_taxaplus.csv")
write.csv(taxa, file="mr16s_taxa.csv")

saveRDS(seqtab.nochim, file="mr16s_seqtab.nochim.rds")
write.csv(seqtab.nochim, file="mr16s_seqtab.nochim.csv")

##If you need to read in previously saved datafiles
setwd("~/moorea_holobiont/mr_16S/")
seqtab.nochim <- readRDS("mr16s_seqtab.nochim.rds")
taxa <- readRDS("mr16s_taxa.rds")
taxa.plus <- readRDS("mr16s_taxaplus.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library(cowplot)

#import dataframe holding sample information
samdf<-read.csv("mr16s_sampledata.csv")
head(samdf)
rownames(samdf) <- samdf$id

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps

#first look at data
ps_glom <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="site_zone", fill="Family")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences
#path='~/moorea_holobiont/mr_16s/mr16s.fasta'
#uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)

colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#### rarefy #####
library(vegan)

rarecurve(seqtab.nochim,step=100,label=FALSE) #after clustering & trimming

total <- rowSums(seqtab.nochim)
subset(total, total <13603)
#9 samples

row.names.remove <- c("A8","B5","B6","B8","C2","F9","H4","H5")
seqtab.less <- seqtab.nochim[!(row.names(seqtab.nochim) %in% row.names.remove),]
samdf.rare <- samdf[!(row.names(samdf) %in% row.names.remove), ]
#85 samples left

seqtab.rare <- rrarefy(seqtab.less,sample=13603)
rarecurve(seqtab.rare,step=100,label=FALSE)

#phyloseq object but rarefied
ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
               sample_data(samdf.rare), 
               tax_table(taxa2))
ps.rare

#saving
write.csv(seqtab.rare, file="mr16s_seqtab.rare_13.6k.csv")
write.csv(samdf.rare, file="mr16s_samdf.rare_13.6k.csv")

#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).
plot_richness(ps, x="site", measures=c("Shannon", "Simpson"), color="zone") + theme_bw()

##if I had more missing data:
#ps_pruned <- prune_taxa(taxa_sums(ps) > 0, ps)
df <- data.frame(estimate_richness(ps.rare, split=TRUE, measures=c("Shannon","InvSimpson")))
df
df <- data.frame(estimate_richness(ps, split=TRUE, measures=c("Shannon","InvSimpson")))
df

df$id <- rownames(df)
df.div <- merge(df,samdf,by="id") #add sample data
df.div <- merge(df,samdf.rare,by="id") #add sample data

#write.csv(df,file="~/moorea_holobiont/mr_16S/mr16s_diversity.csv") #saving
#df <- read.csv("~/Desktop/mrits/mr16s_diversity.csv") #reading back in 

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("zone","site"))
diver.si <- summarySE(data=df.div,measurevar=c("InvSimpson"),groupvars=c("zone","site"))

quartz()
gg.sh <- ggplot(diver.sh, aes(x=site, y=Shannon,color=zone,fill=zone,shape=zone))+
  geom_errorbar(aes(ymin=Shannon-sd,ymax=Shannon+sd),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  guides(fill=guide_legend(title="Reef zone"),color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.sh

gg.si <- ggplot(diver.si, aes(x=site, y=InvSimpson,color=zone,fill=zone,shape=zone))+
  geom_errorbar(aes(ymin=InvSimpson-sd,ymax=InvSimpson+sd),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=zone, shape=zone),size=4,position=position_dodge(0.5))+
  xlab("Site")+
  ylab("Inverse Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  scale_fill_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_colour_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  guides(fill=guide_legend(title="Reef zone"),color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.si

wilcox.test(Shannon~zone,data=df.div)
#W = 1101, p-value = 0.8815
#W = 944, p-value = 0.7232 rarefied
wilcox.test(InvSimpson~zone,data=df.div)
#W = 1138, p-value = 0.6656
#W = 956, p-value = 0.6461 rarefied

mnw <- subset(df.div,site=="MNW")
mse <- subset(df.div,site=="MSE")
tah <- subset(df.div,site=="TNW")

wilcox.test(Shannon~zone,data=mnw)
#W = 148, p-value = 0.1485
#W = 127, p-value = 0.1936 rarefied
wilcox.test(InvSimpson~zone,data=mnw)
#W = 158, p-value = 0.06128 . 
#W = 137, p-value = 0.07666 rarefied
wilcox.test(Shannon~zone,data=mse)
#W = 63, p-value = 0.02398 *
#W = 62, p-value = 0.03672 * rarefied
wilcox.test(InvSimpson~zone,data=mse)
#W = 78, p-value = 0.1015
#W = 74, p-value = 0.116 rarefied
wilcox.test(Shannon~zone,data=tah)
#W = 157, p-value = 0.2871
#W = 115, p-value = 0.2588 rarefied
wilcox.test(InvSimpson~zone,data=tah)
#W = 160, p-value = 0.2388
#W = 115, p-value = 0.2588 rarefied

#### trim low abundance ASVs ####
ps.rare.relabun <- transform_sample_counts(ps.rare, function(x) x / sum(x) )
ps.rare.less <- filter_taxa(ps.rare.relabun, function(x) mean(x) > 1e-5, TRUE)

#### Bar-plots ####
top30 <- names(sort(taxa_sums(ps.rare.less), decreasing=TRUE))[1:30]
ps.top30 <- prune_taxa(top30, ps.rare.less)
plot_bar(ps.top30, x="site_zone", fill="Class") 

ps_glom <- tax_glom(ps.rare.less, "Family")
#this step gets rid of NAs which we don't want
ps1 <- merge_samples(ps_glom, "site_zone")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="site_zone", fill="Family")

#### reorganized up to this point ####
totalsums <- colSums(seqtab.rare)
summary(totalsums)

ps.rare.13 <- prune_taxa(taxa_sums(ps.rare)>13,ps.rare) #remove taxa with <0.1% of reads total
ps.rare.13

#### DESEQ to find differentially abundant ASVs ####
library(DESeq2)

ps.mnw = subset_samples(ps.rare.less, site=="MNW")
ps.mnw.0 <- prune_taxa(taxa_sums(ps.mnw)>0,ps.mnw) #remove samples with 0 total

ds.mnw = phyloseq_to_deseq2(ps.mnw.0, ~ zone)
dds.mnw <- estimateSizeFactors(ds.mnw,type="poscounts")
stat.mnw = DESeq(dds.mnw, test="Wald", fitType="parametric")

res = results(stat.mnw, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mnw = res[which(res$padj < alpha), ]
sigtab.mnw = cbind(as(sigtab.mnw, "data.frame"), as(tax_table(ps.mnw)[rownames(sigtab.mnw), ], "matrix"))
head(sigtab.mnw)
dim(sigtab.mnw)
write.csv(sigtab.mnw,"sigtab.mnw.csv")
#6 for MNW after rarefying & removing low quality samples! 
#1 is just a chloroplast though
#removing chloroplasts
row.names.remove <- c("sq5")
sigtab.mnw <- sigtab.mnw[!(row.names(sigtab.mnw) %in% row.names.remove), ]

goodtaxa <- c(row.names(sigtab.mnw))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
ps.mnw = subset_samples(ex1, site=="MNW")

quartz()
plot_bar(ps.mnw, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea NW")+
  facet_wrap(~class)

ps.mse = subset_samples(ps.rare, site=="MSE")
ds.mse = phyloseq_to_deseq2(ps.mse, ~ zone)
dds.mse <- estimateSizeFactors(ds.mse,type="poscounts")
stat.mse = DESeq(dds.mse, test="Wald", fitType="parametric")

res = results(stat.mse, cooksCutoff = FALSE)
alpha = 0.05
sigtab.mse = res[which(res$padj < alpha), ]
sigtab.mse = cbind(as(sigtab.mse, "data.frame"), as(tax_table(ps.mse)[rownames(sigtab.mse), ], "matrix"))
dim(sigtab.mse)
write.csv(sigtab.mse,"sigtab.mse.csv")
#9 for MSE 
#removing 2 chloroplasts
row.names.remove <- c("sq5","sq29")
sigtab.mse <- sigtab.mse[!(row.names(sigtab.mse) %in% row.names.remove), ]

goodtaxa <- c(row.names(sigtab.mse))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
ps.mse = subset_samples(ex1, site=="MSE")

quartz()
plot_bar(ps.mse, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea SE")+
  facet_wrap(~class)

ps.t = subset_samples(ps.rare, site=="TNW")
ds.t = phyloseq_to_deseq2(ps.t, ~ zone)
dds.t <- estimateSizeFactors(ds.t,type="poscounts")
stat.t = DESeq(dds.t, test="Wald", fitType="parametric")

res = results(stat.t, cooksCutoff = FALSE)
alpha = 0.05
sigtab.t = res[which(res$padj < alpha), ]
sigtab.t = cbind(as(sigtab.t, "data.frame"), as(tax_table(ps.t)[rownames(sigtab.t), ], "matrix"))
dim(sigtab.t)
write.csv(sigtab.t,"sigtab.tah.csv")
#8 for tahiti finally
row.names.remove <- c("sq3") #doesn't even blast to bacteria
sigtab.t <- sigtab.t[!(row.names(sigtab.t) %in% row.names.remove), ]
tax_table(ps.rare)
#visualizing deseq results - can't do with no results
theme_set(theme_cowplot())
#scale_fill_discrete <- function(palname = "Set1", ...) {
#  scale_fill_brewer(palette = palname, ...)
#}
# Phylum order
#x = tapply(sigtab.mnw$log2FoldChange, sigtab.mnw$phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab.mnw$Phylum = factor(as.character(sigtab.mnw$phylum), levels=names(x))
# class 

x = tapply(sigtab.mnw$log2FoldChange, sigtab.mnw$class, function(x) max(x))
x = sort(x, TRUE)
sigtab.mnw$class = factor(as.character(sigtab.mnw$class), levels=names(x))
quartz()

sigtab.mnw$family <- sub("Xanthobacteraceae","Xanthobacteraceae         ",sigtab.mnw$family)
#gg.mnw <- ggplot(sigtab.mnw, aes(x=class, y=log2FoldChange,color=class)) + 
geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-45,hjust = 0, vjust=0.5),legend.position="none")+
  xlab("Class")+
  ylab("Log2 fold change")+
  ggtitle("Moorea NW")#+
gg.mnw <- ggplot(sigtab.mnw, aes(x=family, y=log2FoldChange,color=family)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-90,hjust = 0, vjust=0.5), text=element_text(family="Gill Sans MT"),legend.position="none")+
  xlab("Family")+
  ylab("Log2 fold change")+
  ggtitle("Moorea NW")#+
#  scale_color_manual(name="Reef zone",aesthetics = "color",values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))
gg.mnw

x = tapply(sigtab.mse$log2FoldChange, sigtab.mse$class, function(x) max(x))
x = sort(x, TRUE)
sigtab.mse$class = factor(as.character(sigtab.mse$class), levels=names(x))

gg.mse <- ggplot(sigtab.mse, aes(x=class, y=log2FoldChange,color=class)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-45,hjust = 0, vjust=0.5),legend.position="none")+
  xlab("Class")+
  ggtitle("Moorea SE")+
  ylab("Log2 fold change")
gg.mse <- ggplot(sigtab.mse, aes(x=family, y=log2FoldChange,color=family)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-90,hjust = 0, vjust=0.5),text=element_text(family="Gill Sans MT"),legend.position="none")+
  xlab("Family")+
  ggtitle("Moorea SE")+
  ylab("Log2 fold change")
gg.mse

x = tapply(sigtab.t$log2FoldChange, sigtab.t$class, function(x) max(x))
x = sort(x, TRUE)
sigtab.t$class = factor(as.character(sigtab.t$class), levels=names(x))

sub()

gg.t <- ggplot(sigtab.t, aes(x=class, y=log2FoldChange,color=class)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-45,hjust = 0, vjust=0.5),legend.position="none")+
  xlab("Class")+
  ggtitle("Tahiti NW")+
  ylab("Log2 fold change")
gg.t <- ggplot(sigtab.t, aes(x=family, y=log2FoldChange,color=family)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle=-90,hjust = 0, vjust=0.5),text=element_text(family="Gill Sans MT"),legend.position="none")+
  xlab("Family")+
  ggtitle("Tahiti NW")+
  ylab("Log2 fold change")
gg.t

quartz()
library(ggpubr)
ggarrange(gg.mnw,gg.mse,gg.t,nrow=1)

goodtaxa = c("sq4","sq16","sq26")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)
mnw1 <- subset_samples(ex1,site=="MNW")
quartz()
plot_bar(mnw1, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45))+
  ggtitle("MNW")

goodtaxa = c("sq16","sq26","sq54","sq56","sq65")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex2 = prune_taxa(allTaxa, ps.rare)
mnw <- subset_samples(ex2,site=="MNW")

goodtaxa = c("sq16","sq26","sq54","sq56","sq65")
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex2 = prune_taxa(allTaxa, ps.rare)
mnw <- subset_samples(ex2,site=="MNW")

quartz()
plot_bar(mnw, x="zone", fill="genus")+
  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  ggtitle("Moorea NW")+
  facet_wrap(~genus)


#### DESEQ to normalize counts ####
#BiocManager::install("DESeq2")
library("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("DEFormats")
library(edgeR)
library(DEFormats)
library("MCMC.OTU")
library("cowplot")

#count data
seq <- read.csv("mr16s_AllOTUs_mcmc.csv")
#getting rid of crap ones
#mcmc.otu requires a blank column with sample names first, did manually in excel

seq2 <- purgeOutliers(seq,count.columns=10:2285,sampleZcut=-2.5,otu.cut=0.0001)
#bad samples: 304 (F9)
#50 OTUs pass frequency cutoff

row.names(seq2) <- seq2$id #setting my sample names as row names so I can remove the sample name column
seq2 <- seq2[,10:263] #now removing sample name columns
seqt <- t(seq2) #deseq requires count data in columns

#sample data
sam <- read.csv("mr16s_sampledata.csv")
row.names(sam) <- sam$id
row.names.remove <- c("F9")
sam <- sam[!(row.names(sam) %in% row.names.remove), ]

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

ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
colnames(vstt)<-ids
##replace sequences with shorter names (correspondence table output below)
row.names(taxid)<-ids

ps.norm <- phyloseq(otu_table(vstt, taxa_are_rows=FALSE), 
                    sample_data(sam), 
                    tax_table(taxid))
ps.norm
tax_table(ps.norm)
# sam$new.norm <- sample_sums(ps.norm)
# plot(new.norm~site_zone,data=sam) #more normalized now
# a1 <- aov(new.norm~site_zone,data=sam)
# summary(a1) #p = 0.0646 lol
# plot(new.norm~island,data=sam) #even
# plot(new.norm~site,data=sam) #MSE still a little more
# plot(new.norm~in_off,data=sam) #outer more
# a1 <- aov(new.norm~in_off,data=sam)
# summary(a1) # p = significant

#### rarefy ####
library(vegan)
library(MCMC.OTU)
library(phyloseq)
setwd("~/moorea_holobiont/mr_16S/")
samdf<-read.csv("mr16s_sampledata.csv")
rownames(samdf) <- samdf$id

rarecurve(seqtab.nochim,step=100,label=FALSE)
total <- rowSums(seqtab.nochim)
subset(total, total <10000)
row.names.remove <- c("A5","B5","B6","B8","C2","C8","D12","F8","F9","G9","H3","H4","H5","H6")
seqtab.nochim <- seqtab.nochim[!(row.names(seqtab.nochim) %in% row.names.remove), ]
samdf <- samdf[!(row.names(samdf) %in% row.names.remove), ]

seq.rare <- rrarefy(seqtab.nochim,sample=10000)
rarecurve(seq.rare,step=100,label=FALSE)

write.csv(seq.rare,"~/moorea_holobiont/mr_16S/seq.rare10k_2.csv")
#added the sample column that mcmcotu requires
seq.rare <- read.csv("~/moorea_holobiont/mr_16S/seq.rare10k_2.csv")
str(seq.rare)

seq2 <- purgeOutliers(seq.rare,count.columns=3:2278,sampleZcut=-2.5,otu.cut=0.001)
#keeps 50 OTUs

row.names(seq2) <- seq2$sample #setting my sample names as row names so I can remove the sample name column
seq2 <- seq2[,3:52] #now removing sample name columns
row.names(samdf) <- samdf$id

ps.rare <- phyloseq(otu_table(seq2, taxa_are_rows=FALSE), 
                    sample_data(samdf), 
                    tax_table(taxid))

#making ids easier
ids <- paste0("sq", seq(1, length(colnames(seq2))))
colnames(seq2)<-ids
##replace sequences with shorter names (correspondence table output below)
taxa_names(ps.rare)<-ids

#### Pcoa #####
write.csv(vstt,"~/moorea_holobiont/mr_16S/vstt.16s.csv")
write.csv(seq.rare,"~/moorea_holobiont/mr_16s/seq.rare10k.csv")

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
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  labs(color="Reef zone",shape="Reef zone")
#geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
#coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")
#scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
#theme(legend.position="none")

#after rarefying
rare_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#vst_pcoa <- ordinate(ps.rare, method="MDS", distance="euclidean")
#^this looks horrible
eigen_vals <- rare_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
quartz()
plot_ordination(ps.rare, rare_pcoa, color="zone",shape="zone") + 
  geom_point(size=1) +
  facet_wrap(~site)+
  stat_ellipse(level=0.95)+
  theme_cowplot()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(16,15),labels=c("Backreef","Forereef"))+
  labs(color="Reef zone",shape="Reef zone")
#geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
#coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")
#scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
#theme(legend.position="none")

#trying to select only the variables I want
#resource for analysis here:
#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

ds.all = phyloseq_to_deseq2(ps.rare, ~ site*zone)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")

res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all= cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare)[rownames(sigtab.all), ], "matrix"))
head(sigtab.all)
dim(sigtab.all)

goodtaxa <- c(row.names(sigtab.all))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)

quartz()
plot_bar(ex1, x="zone", fill="class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45))+
  ggtitle("MNW")+
  facet_wrap(~site+genus)


#### bar plot
samdf$newname <- paste(samdf$site_zone,samdf$sample,sep="_")
row.names(samdf) <- samdf$newname
row.names(seq2) <- samdf$newname

taxid2 <- tax_table(ps.rare)

ps.rare2 <- phyloseq(otu_table(seq2, taxa_are_rows=FALSE), 
                     sample_data(samdf), 
                     tax_table(taxid2))

top3 <- names(sort(taxa_sums(ps.rare2), decreasing=TRUE))[1:118]
ps.top3 <- transform_sample_counts(ps.rare2, function(OTU) OTU/sum(OTU))
ps.top3 <- prune_taxa(top3, ps.top3)
quartz()
plot_bar(ps.top3, x="newname",fill="class") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#summed bar plot
  ps <- tax_glom(ps.rare, "family")
  ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
  ps1 <- merge_samples(ps0, "site_zone")
  ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))

ps.melt <- psmelt(ps2)
  
  ps.melt$site_zone <- gsub(1,"MNW-B",ps.melt$site_zone)
  ps.melt$site_zone <- gsub(2,"MNW-F",ps.melt$site_zone)
  ps.melt$site_zone <- gsub(3,"MSE-B",ps.melt$site_zone)
  ps.melt$site_zone <- gsub(4,"MSE-F",ps.melt$site_zone)
  ps.melt$site_zone <- gsub(5,"TNW-B",ps.melt$site_zone)
  ps.melt$site_zone <- gsub(6,"TNW-F",ps.melt$site_zone)
  
ggplot(ps.melt,aes(x=site_zone,y=Abundance,fill=family))+
    geom_bar(stat="identity", colour="black")+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle=0),text=element_text(family="Gill Sans MT"))+
  xlab('Site')+
  labs(fill="Family")
  
quartz()
tb2 <- tb%>%
  group_by(site_zone,family)%>%
  summarize_at("Abundance",sum)
plot_bar(ps.all, x="site_zone",fill="family") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#############################################

# Principal coordinate analysis 

library(vegan)
#install.packages("MCMC.OTU")
library(MCMC.OTU)
#install.packages("ggfortify")
library(ggfortify)
library(cluster)
#install.packages('labdsv')
library(labdsv)

#have to add a first column with no name numbered 1-93
alldat<-read.csv("mr16s_AllOTUs_mcmc.csv")
names(alldat)
str(alldat)
head(alldat) 
#dat<-alldat[c(1:3,5:110),] #to remove rows with missing data, but I don't have any

# purging under-sequenced samples; and OTUs represented in less than 3% of all samples (~must occur at least 3 separate times)
goods2=purgeOutliers(alldat,count.columns=10:2285,sampleZcut=-2.5,otu.cut=0.032)

# creating a log-transfromed normalized dataset for PCoA:
goods.log=logLin(data=goods2,count.columns=10:length(names(goods2)))

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist=vegdist(goods.log,method="manhattan")
goods.pcoa=pcoa(goods.dist)

scores=goods.pcoa$vectors
conditions=goods2[,1:9]
margin=0.01
quartz()
plot(scores[,1], scores[,2],type="n",
     xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     #mgp=c(2.3,1,0),
     xlab=paste("Axis1 (",round(goods.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
     ylab=paste("Axis2 (",round(goods.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
     main="")
points(scores[conditions$zone=="in",1],scores[conditions$zone=="in",2],pch=1)
points(scores[conditions$zone=="out",1],scores[conditions$zone=="out",2],pch=19)
legend("bottomright", c("Inshore","Offshore"), pch=c(1, 19), cex=0.8)

#prettier PCA
scores=goods.pcoa$vectors
conditions=goods[,1:9]
zone <- conditions$zone

p1 <- ggordiplots::gg_ordiplot(ord=goods.pcoa$vectors,groups=zone,scaling=0.8,choices=c(1,2),conf=0.95,spiders=TRUE,pt.size=2)
plot <- p1$plot
plot + xlab('Axis1 (19.%)')+ 
  ylab('Axis2 (14.42%)')+
  scale_shape_manual(values=c(1,5))
names(p1)
spid <- p1$df_spiders
p1$df_spiders
gg <- p1$df_ord
quartz()
ggplot(gg,aes(x=x,y=y,group=Group,color=Group,shape=Group,fill=Group))+
  geom_point(size=2)+
  geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.95)+
  xlab('Axis 1 (19.5%)')+
  ylab('Axis 2 (14.4%)')+
  theme_classic()+
  theme(text=element_text(family="Gill Sans MT"))+
  scale_shape_manual(values=c(21,24),labels=c("Inshore","Offshore"))+
  scale_fill_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  scale_color_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
  labs(shape="Reef zone",color="Reef zone",fill="Reef zone")
