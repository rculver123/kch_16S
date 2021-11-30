############################################
# There are lots of options for normalizing data
# You only ever need to normalilze ONCE,
# in other words, don't log2 fold change your data AND plot realtive abundance

############################################

# I make a folder called clean_data where I store all my filtered/modified phyloseq objects

clean_data_dir <- 'clean_data'


############################################
# Option 1: Relative abundance.
# Advantages: easy to interprut, can deal with samples with a wide range of read depths
# Disadvantages: does not do well with unique singletons or PCR bias

############################################

ps.relabs.raw <- transform_sample_counts(ps.filt, function(x) x / sum(x))                          

# You can try to filter out random singletons by doing something like one of these two options:
# PS: you'll need library(dplyr) to use the %>% function

# Filter out low prevalencesingletons
ps.relabs <- filter_taxa(ps.filt, function(x) sum(x > 5) > (0.05*length(x)), TRUE) %>% #Must have 5 reads in 5% of samples
    transform_sample_counts(., function(x) x / sum(x)) 

# Deal with singletons by merging highly similar ASVs
ps.relabs.glom <- tip_glom(ps.filt, h=0.1) %>% # Merges highly similar ASVs
    transform_sample_counts(., function(x) x / sum(x)) 

############################################
# Option 2: Log2 fold change
############################################

ps.relabs.raw <- transform_sample_counts(ps.filt, function(x) log2(x))                          

 