# Ryan Ward
# Feb 10, 2023
#
# This code filters a binary matrix of gene presence/absence in 
# Alphaproteobacteria to identify unique genes in individual 
# organisms. The final result is a vector of indices for rows 
# with only one "1" value, indicating the presence of unique 
# genes in a single Alphaproteobacteria.

# This code creates a data table that summarizes the presence 
# of genes in a set of Alphaproteobacteria. It loops through 
# a vector of indices (uniqueome_indices), which represent 
# unique combinations of gene presence/absence in the 
# Alphaproteobacteria. For each index, the code uses the 
# "extract_comb" and "melt" functions to convert the binary 
# matrix of gene presence/absence into a more usable format. 
# The resulting data table is filtered to include only those 
# genes that are present in the current Alphaproteobacteria. 
# The filtered data table is then added to the final 
# "genes_per_genome" data table, which will contain 
# information about the presence of all genes in all 
# Alphaproteobacteria.

# This code then summarizes the number of unique and shared 
# genes among the Alphaproteobacteria and creates a bar plot 
# to visualize the results. The plot shows the number of genes 
# present in all or one of the Alphaproteobacteria species.

library(conflicted)
library(pacman)

source("upset_alphas.R")

# Use pacman to load the following packages
p_load(data.table, tidyverse, broom, modelr, Hmisc, conflicted)
p_load_current_gh("DoseResponse/drcData", "ryandward/drc", "hrbrmstr/hrbrthemes")

# Conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarize", "dplyr")

# Define a custom theme for the ggplot visualization
doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12, 
	axis_title_size = 12, 
	axis_col = "black")

uniqueome_indices <- comb_name(presence.mixed.matrix) %>% 
	as_tibble %>% 
	mutate(comb_level = lengths(regmatches(value, gregexpr("1", value)))) %>% 
	mutate(index = 1:nrow(.)) %>% 
	filter(comb_level == 1) %>% 
	pull(index)

genes_per_genome <- data.table()
for (i in uniqueome_indices) {
	one_result <- presence[extract_comb(presence.mixed.matrix, comb_name(presence.mixed.matrix[, i])),] %>% 
		data.table %>% 
		melt(id.vars = "HOG", variable.name = "species", value = "presence") %>% 
		filter(presence == TRUE)
	setDT(one_result)
	set(one_result, j = "presence_index", value = i)
	genes_per_genome <- rbind(genes_per_genome, one_result)
}

orthos_regroup <- orthos %>% 
	separate(
		genome, c("Genus", "Species", "Strain"), 
		remove = FALSE, 
		extra = "merge") %>%
	mutate(species = paste(Genus, Species)) %>% 
	data.table

genes_per_genome <- orthos_regroup[genes_per_genome, on = .(HOG, species)]

alpha_summary <- genes_per_genome %>% group_by(species) %>% summarise(unique_genes = n())

alpha_summary <- orthos_regroup %>% 
	filter(species %in% alpha_summary$species) %>% 
	group_by(species) %>% 
	summarise(all_genes = n()) %>% 
	inner_join(alpha_summary)

alpha_summary <- presence[extract_comb(presence.mixed.matrix, comb_name(presence.mixed.matrix[, 1])), ] %>% 
	data.table %>% 
	melt(id.vars = "HOG", variable.name = "species", value = "presence") %>% 
	filter(presence == TRUE) %>% 
	select(HOG, species) %>% 
	unique %>% 
	inner_join(orthos_regroup, multiple = "all") %>% 
	group_by(species) %>% 
	summarise(shared_genes = n()) %>% 
	inner_join(alpha_summary)

alpha_summary <- alpha_summary %>% data.table %>% melt(id.vars = "species", variable.name = "set", value.name = "gene count") %>%
	mutate(set = factor(set, levels = c("all_genes", "shared_genes", "unique_genes")))

alpha_order <- alpha_summary %>% filter(set == "all_genes") %>% arrange(`gene count`) %>% pull(species)

alpha_summary$species <- factor(alpha_summary$species, levels = alpha_order)

alpha_summary %>% ggplot(aes(fill = set, y = `gene count`, x = species)) + 
	geom_bar(position = "dodge", stat = "identity") +
	doc_theme +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

alpha_summary %>% group_by(species) %>% 
	mutate(`gene count` = `gene count`/`gene count`[set == "all_genes"]) %>%
	filter(set != "all_genes") %>%
	rename("genome proportion" = "gene count") %>%
	ggplot(aes(fill = set, y = `genome proportion`, x = species)) + 
	geom_bar(position = "dodge", stat = "identity") +
	doc_theme +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

