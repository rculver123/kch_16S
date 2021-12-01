############################################
# So you want to make an ordination
############################################


############################################
# Load in data
############################################

# I make a folder called figues, then another follder within it called barplots where I'm going to save all my barplots
ord_dir <- 'figures/ordinations'

# Load your libraries!
library(phyloseq)
library(dplyr)
library(ggplot2)

# Load in your phyloseq object
ps <- readRDS(paste0(clean_data_dir,'filtered_phlyoseq.RDS'))  


############################################
# Define functions needed - alter "SeqName" to be whatever you call
############################################

get_evals <- function(pcoa_out) {
  evals <- pcoa_out$values[,1]
  var_exp <- 100 * evals/sum(evals)
  return(list("evals" = evals, "variance_exp" = var_exp))
}

############################################
# Unifrac - basics
############################################


ps2ordinate <- subset_samples(ps, !SampleType %in% 'Cecum') # Choose which samples you want to ordinate

# here weighted=TRUE, change to FALSE if you want unweighted
pcoa_weighted <- ordinate(ps2ordinate,  method = "MDS", distance = "unifrac", weighted=TRUE)

var_exp <- get_evals(pcoa_weighted)$variance_exp

# Define colors
colors.community <- c('red','black','blue','green','yellow')

# Plot
plot_ordination(ps, pcoa_canberra, color="Community", shape="Community") +
facet_wrap(~Location) + # divides into multiple plots by variable
scale_color_manual(values=colors.community)+
theme_minimal()+
coord_fixed(sqrt(var_exp[2] / var_exp[1])) #adjust axes to be proportional to variance explained

ggsave(paste0(ord_dir, 'weighted_unifrac.pdf'),width=8,height=4)
