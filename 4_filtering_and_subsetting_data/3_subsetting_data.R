############################################
# You don't always want to plot everything,
# so here are ways of subsetting data
############################################

# the following column names are based on the sample data, but you should
# change it to whatever your column names are/ catgories you want to subset by

# Subset by column named Project
ps.study<-subset_samples(ps, Project %in% 'MouseExp1') %>%
filter_taxa(., function(x) sum(x) >0, TRUE) # Get rid of taxa that are now 0

# Subset by colmn named Project but you want multiple projects
ps.study<-subset_samples(ps, Project %in% c('MouseExp1','Cultures')) %>%
filter_taxa(., function(x) sum(x) >0, TRUE) # Get rid of taxa that are now 0

# Using 'OR' operator
final.days<- subset_samples(ps, (Day %in% '36' & Community %in% c('Small Intestine','Stool')) | Day %in% c('0','49')) %>%
filter_taxa(., function(x) sum(x) >0, TRUE) # Get rid of taxa that are now 0
