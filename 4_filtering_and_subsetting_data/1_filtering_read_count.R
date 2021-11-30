############################################
# Step 1:  Filter out samples with low reads
############################################

# I make a folder called clean_data where I store all my filtered/modified phyloseq objects

clean_data_dir <- 'clean_data'
# Look at reads distribution
sample_sum_df <- data.frame(sum = sample_sums(ps))
ggplot(sample_sum_df, aes(x = sum)) +
  geom_histogram(color = "black", fill = "indianred", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# In this example, I set my read threshold to be 5000. You can change this to whatever you like
# Remove reads and samples that were messed up during prep
ps.filt <- prune_samples(sample_sums(ps) > 5000, ps.study)
ps.filt
saveRDS(ps.filt, paste0(clean_data_dir,'filtered_phlyoseq.RDS'))                                  

