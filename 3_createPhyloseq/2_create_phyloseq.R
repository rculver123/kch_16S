############################################
#
# R scripts!!
#
############################################
# A phyloseq is an R object that is commonly used in sequencing analysis
#
# Here is a toy example that goes through some of its uses: https://joey711.github.io/phyloseq/preprocess.html
#
# Briefly, it is an object that stores:
# otu_table <- table of read counts for each ASV
# tax_table <- taxonomic assignment of ASVs
# sample_data <- metadata for all your experiments
# phy_tree <- phylogenetic tree from your data
#
############################################
# Step 1: Create a .csv table of all your metadata
# ESSENTIAL: There must be a column that matches the sample names of seq files on Sherlock (typicallly S followed by a number; S1, S2, etc.)
############################################

# View sample.table.mouse.csv for an example
# Note the SeqName column has the sample #s (S1, S2, etc)

############################################
# Step 2: Install & load phyloseq package
############################################

# These are R commands!!
# If it doesn't work, google it. There's lots of install options.

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

# Load the library
library(phyloseq)

############################################
# Step 3: Create phyloseq
############################################

# Change raw_data_dir to wherever your folder is of your raw data.

raw_data_dir <- 'sherlock_data'

# Load seqtab, taxa, and phylogenetic tree
seqtab.nochim <- readRDS(paste0(raw_data_dir,'seqtab.nochim.RDS'))
taxa <- readRDS(paste0(raw_data_dir,'taxa.rds'))
phy_tree <- read_tree(paste0(raw_data_dir,'phy_tree.tree'))

# Option 1: Simple load sample data frame simple 

samdf <- read.csv(paste0(raw_data_dir,'SImouse_exp1_samples.csv')) # read csv
rownames(samdf) <- samdf$SeqName #change row names to be the sequencing name!!

# Option 2: Advanced load sample data frame
library(reshape2) #install it if you don't have it
library(dplyr) #install it if you don't have it
samdf <- read.csv(paste0(raw_data_dir,'SImouse_exp1_samples.csv')) %>%
    filter(!SeqName == '') %>% # remove samples that were not sequenced
    mutate(SampleType = factor(SampleType, levels=c('Proximal','Mid','Distal','Cecum','Stool'))) %>% # tells R order to plot things in
    
    mutate(Location = factor(Location, levels=c('Small Intestine','Cecum','Stool','SI Community','Stool Community'))) %>% # tells R order to plot things in
    mutate(finalDaySamples = ifelse(Day %in% '36' & Community %in% c('Small Intestine','Stool'),1,
                             ifelse(Day %in% '49' & Community %in% c('SI2Stool','Stool2SI'),1,0))) %>% # Creates a new column of my final days of samples
    mutate(Community = factor(Community, levels=c('Small Intestine','SI2Stool','Stool','Stool2SI'))) # tells R order to plot in
rownames(samdf) <- samdf$SeqName


# Create the phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxa), phy_tree(phy_tree))
#store DNA sequences in refseq slot of phyloseq object and rename taxa to a short taxa names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps

saveRDS(ps, paste0(raw_data_dir, 'phyloseq_raw.rds'))

