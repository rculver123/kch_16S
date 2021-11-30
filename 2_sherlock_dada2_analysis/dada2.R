library(dada2); packageVersion("dada2")

# Edit the following lines
per_sample_folder_fp <- "path" #THIS SHOULD BE THE FOLDER CONTAINING THE FWD AND REV_READS FOLDERS


fwd_path <- file.path(per_sample_folder_fp,"fwd_reads/samples")
rev_path <- file.path(per_sample_folder_fp,"rev_reads/samples")

filt_pathF  <- file.path(fwd_path, "filtered")
filt_pathR <- file.path(rev_path, "filtered")

fnFs <- sort(list.files(fwd_path, pattern="fastq"))
fnRs <- sort(list.files(rev_path, pattern="fastq"))


filtFs <- file.path(filt_pathF, fnFs)
filtRs <- file.path(filt_pathR, fnRs)

if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")
out <- filterAndTrim(fwd = file.path(fwd_path, fnFs), filt =filtFs, rev = file.path(rev_path, fnRs), filt.rev =filtRs,
                      truncLen = c(250, 180),
                      trimLeft = c(2, 2),
                      maxN = 0,
                      maxEE = c(2,2),
                      truncQ = 2,
                      rm.phix = TRUE,
                      compress = TRUE,
                      multithread = TRUE,
                      verbose = TRUE)

saveRDS(out, "out.RDS")

filtpathF <- file.path(per_sample_folder_fp,"fwd_reads/samples/filtered")
filtpathR <- file.path(per_sample_folder_fp,"rev_reads/samples/filtered")

filtFs <- list.files(filtpathF, pattern=".fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern=".fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "[.]"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "[.]"), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Learn forward error rates
set.seed(357167)
errF <- learnErrors(filtFs, nbases=1e9, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 15)
saveRDS(errF, "errF.RDS")

# Learn reverse error rates
set.seed(801295)
errR <- learnErrors(filtRs, nbases = 1e9, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 15)
saveRDS(errR, "errR.RDS")


derepF<-derepFastq(filtFs)
derepR<-derepFastq(filtRs)

#run dada2
dadaFs <- dada(derepF, err=errF, multithread=TRUE, pool = "pseudo")
dadaRs <- dada(derepR, err=errR, multithread=TRUE, pool = "pseudo")

saveRDS(dadaFs,"dadaFs.RDS")
saveRDS(dadaRs,"dadaRs.RDS")
# 
mergers <- mergePairs(dadaFs, derepF, dadaRs, derepR, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



############# Make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "seqtab.RDS")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, "seqtab.nochim.RDS")





############# Assign taxa
library(dplyr)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/groups/kchuang/16S_version_Becca/silva_nr_v132_train_set.fa", multithread=TRUE)

#make species level assginments based on EXACT MATCHING between ASVs and reference strains
taxa <- addSpecies(taxa, "/home/groups/kchuang/16S_version_Becca/silva_species_assignment_v132.fa")

saveRDS(taxa, "taxa.rds")





########### Make Tree

library (DECIPHER)
library(phangorn)
library(phyloseq)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

write.tree(phyloseq(phy_tree(fitGTR$tree)),'phy_tree.tree')




