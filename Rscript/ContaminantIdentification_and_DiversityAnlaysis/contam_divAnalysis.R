# MICROBIOME ANALYSIS WITH PHYLOSEQ
## written by: Elise Hickman (ehickman@email.unc.edu)
# To install and load phyloseq and microboime
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

library(BiocManager)

BiocManager::install("phyloseq")
BiocManager::install("microbiome")
BiocManager::install("decontam")

# Load libraries:

library(phyloseq)
library(decontam)
library(ggplot2)
library(microbiome)
library(RColorBrewer)
library(vegan)

## DATA IMPORT

# Create phyloseq class object.

biom <- import_biom("Data/raw-OTUtable-with-controls-hdf5.biom", parseFunction = parse_taxonomy_greengenes) 
map <- import_qiime_sample_data("Data/ELFmetadata.txt")
tree <- read_tree_greengenes("Data/rooted-tree.nwk")
data_original <- merge_phyloseq(biom, map, tree)
data <- subset_samples(data_original, 
                       X.SampleID != "SA124V1R" & 
                         X.SampleID != "SA092V1L" & 
                         X.SampleID != "SA075V1R" & 
                         X.SampleID != "SA119V1L" & 
                         X.SampleID != "SA111V1R")

## FILTERING CONTAMINANTS & CONTROLS USING DECONTAM

# Transform data to relative abundance to examine histogram and set threshold.
data_relative <- transform_sample_counts(data, function(OTU) OTU/sum(OTU))
sample_data(data_relative)$is.neg <- sample_data(data_relative)$SampleDescription == "Negative Control"
contamdf.prev <- isContaminant(data_relative, method="prevalence", neg="is.neg", threshold = 0.3)
hist(contamdf.prev$p, breaks = 20)
table(contamdf.prev$contaminant) # This identifies 4677 taxa as non-contaminant and 669 taxa as contaminant

# Make phyloseq object of presence-absence in negative controls and true samples.
ps.pa <- transform_sample_counts(data_relative, function(abund) 1*(abund>0))

# Make data.frame of prevalence in positive and negative samples and graph the data.
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleDescription == "Negative Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleDescription == "ELF", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev$contaminant)

# Plot
th <- theme_set(theme_bw())
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Generate phyloseq object that has contaminants and controls filtered out. 
keep_taxa <- contamdf.prev[contamdf.prev$contaminant == FALSE, ]
new_keep_taxa <- setNames(cbind(rownames(keep_taxa), keep_taxa, row.names = NULL), 
                          c("ID", "freq", "prev", "p.freq", "p.prev", "p", "contaminant"))
new_keep_taxa$ID <- as.character(new_keep_taxa$ID)
keep_taxa_vector <- new_keep_taxa$ID

contaminants_removed_data <- prune_taxa(keep_taxa_vector, data)
contaminants_and_controls_removed_data <- subset_samples(contaminants_removed_data, 
                                                         X.SampleID != "PCRNegativeA" & 
                                                           X.SampleID != "PseudomonasPositive" & 
                                                           X.SampleID != "KITNEG01" & 
                                                           X.SampleID != "KITNEG02" & 
                                                           X.SampleID != "ELFNEG01" & X.SampleID != "ELFNEG02")


# Set theme for plotting alpha diversity
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

# Graphing alpha diversity using non-rarefied data (as recommended by phyloseq). 
# To change the group it is graphed by, change the term for x = to match your group(s) of interest.
plot_richness(contaminants_and_controls_removed_data, x = "Sex") + geom_boxplot()

# Creating metadata object to bind to alpha diversity data
metadata <- data.frame(map)

# Calculate alpha diversity results data frame
alpha_diversity_step1 <- estimate_richness(contaminants_and_controls_removed_data)

# Make sample ID a column instead of a row name
alpha_diversity_step2 <- setNames(cbind(rownames(alpha_diversity_step1), alpha_diversity_step1, row.names = NULL), 
                                  c("X.SampleID", "Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))

# Merged alpha diversity and metadata table
alpha_diversity_data <- merge(alpha_diversity_step2, metadata, by = "X.SampleID")

# This is the table that was exported for further analysis in JMP Pro and Prism. 
# You could also perform statistics in R
write.table(alpha_diversity_data, "Output/alpha_diversity_raw_with_meta.txt", sep = "\t", row.names = F)

## DATA FILTERING

# Data filtering by prevalence, removing all taxa that have fewer than 5 reads total. 
# This data table is what I exported and gave as input data for further analyses to Andrew Hinton. 
filtered_data1 <- filter_taxa(contaminants_and_controls_removed_data, function (x) {sum(x > 0) > 4}, prune=TRUE)
### This is the same as the renamed "OTU_Table_contaminantsFilteredOut.txt" in the \Data directory                                 
                                 
                            

## PLOTTING RELATIVE ABUNDANCE BAR GRAPHS

# Create column in metadata for sex and device combined
sample_data(filtered_data1)$SexBySubjGroup <- paste(sample_data(filtered_data1)$Sex, 
                                                    sample_data(filtered_data1)$SubjectGroup,
                                                    sep = " ")



# Specify the order you want each of the groups to appear for variables of interest
SubjGroupLevelOrder <- c("Nonsmoker", "Ecig", "Smoker")
SexLevelOrder <- c("Male", "Female")
SexBySubjGroupLevelOrder <- c("Male Nonsmoker", "Female Nonsmoker", "Male E Cig", "Female E Cig", "Male Smoker", "Female Smoker")

# Phylum by group (example - can exchange different groups, taxa levels, or top # of taxa)

# Merging "other" taxa that you don't want to show on graph. 
# First use merge_less_than_top function I found on the internet. 
merge_less_than_top <- function(pobject, top=10){
  transformed <- transform_sample_counts(pobject, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Other"
      tax_table(merged)[i,1:7] <- "Other"}
  }
  return(merged)
}

# Then run the function, specifying the number of slots to fill and taxa.
ps3 <- tax_glom(filtered_data1, "Phylum")
ps3.top <- merge_less_than_top(ps3, top=4)
ps4 <- merge_samples(ps3.top, "SubjectGroup")
ps5 <- transform_sample_counts(ps4, function(x) x/sum(x))

ps6 <- psmelt(ps5)
ps6$Phylum <- as.character(ps6$Phylum)
colourCount = length(unique(ps6$Phylum))

# Note, for Genus set brewer pal to 9 and Set1, for Phylum, 5 and Dark2
getPalette = colorRampPalette(brewer.pal(5, "Dark2"))

plot <- plot_bar(ps5, fill = "Phylum") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", size = 1.5),
        axis.title = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 360, hjust = 0.5, size = 14, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color = "black", size = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  scale_fill_manual(values = getPalette(colourCount))
plot$data$Sample <- factor(plot$data$Sample, levels = SubjGroupLevelOrder)
plot 

# These plots then exported and labels "prettied up" using Powerpoint for final figure in paper.
