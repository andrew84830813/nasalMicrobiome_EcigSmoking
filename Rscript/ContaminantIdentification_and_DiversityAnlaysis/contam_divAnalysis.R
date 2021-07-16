# MICROBIOME ANALYSIS WITH PHYLOSEQ
## written by: Elise Hickman (ehickman@email.unc.edu)


## Load/Install Required Packages
library(phyloseq)
library(ggplot2)
library(microbiome)
library(RColorBrewer)
library(vegan)
library(data.table)

setwd("16sorganizedanalysisR")

# Create phyloseq class object.

biom <- import_biom("Data/raw-OTUtable-hdf5.biom", parseFunction = parse_taxonomy_greengenes)
map <- import_qiime_sample_data("Data/ELFmetadata.txt")
tree <- read_tree_greengenes("Data/rooted-tree.nwk")
data_original <- merge_phyloseq(biom, map, tree)
data <- subset_samples(data_original,
                       X.NAME != "SA124V1R" &
                         X.NAME != "SA092V1L" &
                         X.NAME != "SA075V1R" &
                         X.NAME != "SA119V1L" &
                         X.NAME != "SA111V1R")

# Data Filtering by Contaminants

data_relative <- transform_sample_counts(data, function(OTU) OTU/sum(OTU))
sample_data(data_relative)$is.neg <- sample_data(data_relative)$SampleDescription == "Negative Control"
contamdf.prev <- isContaminant(data_relative, method="prevalence", neg="is.neg", threshold = 0.3)
table(contamdf.prev$contaminant)

ps.pa <- transform_sample_counts(data_relative, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleDescription == "Negative Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleDescription == "ELF", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

keep_taxa <- contamdf.prev[contamdf.prev$contaminant == TRUE, ]
new_keep_taxa <- setNames(cbind(rownames(keep_taxa), keep_taxa, row.names = NULL),
                          c("ID", "freq", "prev", "p.freq", "p.prev", "p", "contaminant"))
new_keep_taxa$ID <- as.character(new_keep_taxa$ID)
keep_taxa_vector <- new_keep_taxa$ID

contaminants_removed_data <- prune_taxa(keep_taxa_vector, data)
contaminants_and_controls_removed_data <- subset_samples(contaminants_removed_data,
                                                         +                        X.SampleID != "PCRNegativeA" &
                                                           +                            X.SampleID != "PseudomonasPositive" &
                                                           +                            X.SampleID != "KITNEG01" &
                                                           +                            X.SampleID != "KITNEG02" &
                                                           +                            X.SampleID != "ELFNEG01" & X.SampleID != "ELFNEG02")


# To return list of top x number of [taxa] in a list form.

filtered_data1_genus <- tax_glom(filtered_data1, "Genus")
genustop20 <- names(sort(taxa_sums(filtered_data1_genus), TRUE)[1:20])
filtered_data1_genustop20 <- prune_taxa(genustop20, filtered_data1)

# To relevel a factor.

sample_data(data)$SexBySubjGroup <- factor(sample_data(data)$SexBySubjGroup, levels = c("Male Nonsmoker", "Female Nonsmoker", "Male E Cig", "Female E Cig", "Male Smoker", "Female Smoker"))
print(levels(sample_data(data)$SexBySubjGroup))

# Plot alpha diversity.

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

plot_richness(contaminants_and_controls_removed_data, x = "Sex") + geom_boxplot()

# Saving alpha diversity measures.

metadata <- read.table("ELFmetadata.txt", sep = "\t", header = TRUE)
alpha_diversity <- estimate_richness(contaminants_and_controls_removed_data)
new_alpha_diversity <- setNames(cbind(rownames(alpha_diversity), alpha_diversity, row.names = NULL),
                                c("X.SampleID", "Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alpha_diversity_data <- merge(new_alpha_diversity, metadata, by = "X.SampleID")
write.table(alpha_diversity_data, "alpha_diversity_data_with_meta.txt", sep = "\t", row.names = F)

# Single Graph Plotting with ggplot2.

p <- ggplot(alpha_diversity_data, aes (x = Sex, y = Shannon)) + geom_jitter(height = 0.15, width = 0.15)

# Statistics

kruskal.test(Shannon ~ Sex, data = alpha_diversity_data)
pairwise.wilcox.test(alpha_diversity_data$Shannon, alpha_diversity_data$Sex, p.adjust.method = "BH")

# Set spectral palette with RColorBrewer or make your own palette.
pal <- "Spectral"
axis.label.14.text <- element_text(face = "bold", color = "black", size = 14)

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

# To print out glommed abundances.
ps3 <- tax_glom(contaminants_and_controls_removed_data, "Genus")
ps3.top <- merge_less_than_top(ps3, top=10)
ps4 <- subset_samples(ps3.top, SexBySubjGroup=="Female E Cig")
Female_EC_Genus_Abundance <- summarize_taxa(ps4, "Genus")
write.table(Female_EC_Genus_Abundance, file = "F_EC_G_Abun.txt", sep = "\t", row.names = FALSE)

# Then run the function to make the plot, specifying the number of slots to fill and taxa.
ps3 <- tax_glom(contaminants_and_controls_removed_data, "Phylum")
ps3.top <- merge_less_than_top(ps3, top=4)
ps4 <- merge_samples(ps3.top, "Sex")
ps5 <- transform_sample_counts(ps4, function(x) x/sum(x))

# Set color palette. Note,for Genus run brewer.pal(9, "Set1"); for Phylum brewer.pal(5, "Dark2")
ps6 <- psmelt(ps5)
ps6$Phylum <- as.character(ps6$Phylum)
colourCount = length(unique(ps6$Phylum))
getPalette = colorRampPalette(brewer.pal(6, "Set1"))

# Generate Plot
plot <- plot_bar(ps5, fill = "Phylum") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", size = 1.5),
        axis.title = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(angle = 360, hjust = 0.5, size = 14, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black", size = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + labs (y = "Relative Abundace") +
  scale_fill_manual(values = getPalette(colourCount))
plot$data$Sample <- factor(plot$data$Sample, levels = SexBySubjGroup_Level_Order)
plot$data$Sample <- factor(plot$data$Sample, levels = SubjGroup_Level_Order)
plot

# Data filtering by prevalence, removing all taxa that have fewer than 5 reads total. Can plot histogram of relative abundance, too.

filtered_data1 <- filter_taxa(contaminants_and_controls_removed_data, function (x) {sum(x > 0) > 4}, prune=TRUE)
filtered_data1_relabun  = transform_sample_counts(filtered_data1, function(x) x / sum(x) )
asv_means  <- taxa_sums(filtered_data1_relabun) / nsamples(filtered_data1_relabun)
qplot(asv_means, log = "x")

# Beta Diversity Graphing and Analysis w Standard Methods.
set.seed(4)
filtered_data1_rarefied <- rarefy_even_depth(filtered_data1, rngseed = FALSE, trimOTUs = TRUE)

fd1r_ecig <- subset_samples(filtered_data1_rarefied, SubjectGroup=="Ecig")

unifrac_dist = distance(filtered_data1_rarefied, method = "unifrac")
ordination = ordinate(filtered_data1_rarefied, method="PCoA", distance=unifrac_dist)

plot <- plot_ordination(filtered_data1_rarefied, ordination, color="Sex") + theme(
  aspect.ratio=1,
  plot.title = element_text(hjust = 0.5, size = 16),
  axis.title = element_text(size = 12, color = "black"),
  axis.text.x = element_text(size = 12, color = "black"),
  axis.text.y = element_text(size = 12, color = "black"),
  axis.ticks = element_line(color = "black", size = 1),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12)) +
  ggtitle("Unweighted Unifrac PCoA")

# Plot Beta Diversity for Poster
plot <- plot_ordination(filtered_data1_rarefied, ordination, color="Sex") + theme_bw() + theme(
  aspect.ratio=1,
  plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 20, color = "black", face = "bold", vjust = -0.5),
  axis.text.x = element_text(size = 14, color = "black", face = "bold"),
  axis.text.y = element_text(size = 14, color = "black", face = "bold"),
  axis.ticks = element_line(color = "black", size = 1),
  legend.title = element_text(size = 20, face = "bold"),
  legend.text = element_text(size = 18, face = "bold")) +
  ggtitle("Unweighted Unifrac - Genus") +
  geom_point(size = 3, aes(shape = SubjectGroup, color = SubjectGroup)) +
  stat_ellipse(type = "norm", linetype = 2) +
  scale_colour_manual(values=group_pal) +
  scale_shape_manual(16, 17) +
  guides(shape = FALSE)


adonis(wunifrac_dist ~ sample_data(filtered_data1_rarefied)$SexBySubjGroup)
adonis(wunifrac_dist ~ sample_data(filtered_data1_rarefied)$Sex*sample_data(filtered_data1_rarefied)$SubjectGroup
       + sample_data(filtered_data1_rarefied)$NoseSide)

# Beta Diversity with Tax Glomming
filtered_data1_rarefied_phylum <- tax_glom(filtered_data1_rarefied, "Phylum")

# Getting Taxonomy File in right format for qiime2. First get out tax table from phyloseq and import qiime taxa table.
taxtable = as(tax_table(filtered_data1), "matrix")
taxtabledf <- as.data.frame(taxtable)
taxtabledf2 <- data.frame(as.character(rownames(taxtabledf)), taxtabledf)

colnames(taxtabledf2)[1]="Feature.ID"
taxtabledf2$Feature.ID <- as.character(taxtabledf2$Feature.ID)
keep_taxa_vector_2 <- taxtabledf2$Feature.ID

qiime2taxfiltered <- subset(qiime2tax, Feature.ID %in% keep_taxa_vector_2)

# Plot Relative Abundance for Poster
plot <- plot_bar(ps5, fill = "Genus") + theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              panel.background = element_rect(color = "black", size = 1.5),
                                              axis.title = element_text(size = 20, color = "black", face = "bold"),
                                              axis.text.x = element_text(angle = 315, vjust = 1, size = 20, color = "black", face = "bold"),
                                              axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                                              axis.ticks = element_line(color = "black", size = 1), axis.ticks.length.x = unit(0.25, "inch"),
                                              axis.title.x = element_blank(),
                                              legend.title = element_text(size = 20, face = "bold"),
                                              legend.text = element_text(size = 16, face = "bold")) + scale_fill_manual(values = getPalette(colourCount))
