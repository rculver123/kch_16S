############################################
# So you want to make a barplot???
############################################


############################################
# The basics
############################################

# I make a folder called figues, then another follder within it called barplots where I'm going to save all my barplots
bar_dir <- 'figures/barplots'

# Load your libraries!
library(phyloseq)
library(dplyr)
library(ggplot2)

# First I'm going to subset the data that I want to plot (df2plot = data frame to plot)
df2plot <- subset_samples(ps, finalDays %in% 1) %>% # choosing my final day samples
    filter_taxa(., function(x) sum(x > 5) > (0.05*length(x)), TRUE) %>% # I've decided to filter out low abundance singletons here
    transform_sample_counts(., function(x) x / sum(x)) %>% # Changed everything to relative abundance
    psmelt() # this is a melt function that turns everything into an giant data table

head(df2plot) # this will print out the top 10 lines

# Now I'm going to plot!
ggplot(final.days, aes(x=SeqName,y=Abundance, fill=Phylum)) + # you can change fill to be family, OTU, etc.
    geom_bar(stat='identity') + 
    theme_minimal() +
    facet_wrap(~Community+SampleType, scales='free_x' # this will break my barplots up into categoires
ggsave(paste0(bar_dir,'first_barplot.pdf'),width=6,height=12) # save it to a pdf

############################################
# Filtering
############################################

# This time, I only want Bacteroides in greater than 0.1 abundance!
df2plot <- subset_samples(ps, finalDays %in% 1) %>% # choosing my final day samples
    filter_taxa(., function(x) sum(x > 5) > (0.05*length(x)), TRUE) %>% # I've decided to filter out low abundance singletons here
    transform_sample_counts(., function(x) x / sum(x)) %>% # Changed everything to relative abundance
    psmelt() %>% # this is a melt function that turns everything into an giant data table
    filter(Genus %in% 'Bacteroides') %>% # only include Bacteroides!
    filter(Abundance > 0.1) # only include things >0.1 abundance

head(df2plot) # this will print out the top 10 lines

# Now I'm going to plot!
ggplot(final.days, aes(x=SeqName,y=Abundance, fill=Phylum)) + # you can change fill to be family, OTU, etc.
    geom_bar(stat='identity') + 
    theme_minimal() +
    facet_wrap(~Community+SampleType, scales='free_x' # this will break my barplots up into categoires
ggsave(paste0(bar_dir,'top_bacteroides_barplot.pdf'),width=6,height=12) # save it to a pdf
