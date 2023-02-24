library(pacman)

p_load(data.table, viridis, grid, tidyverse)
p_load_current_gh("jokergoo/ComplexHeatmap")

ecoli_genes <- fread("Escherichia_coli_BW25113.tsv")

alpha_ecoli_presence <- alpha_ecoli_orthos_wider %>% 
	data.table %>% 
	melt(id.vars = c("genus", "genome"), variable.name = "HOG", value.name = "N") %>% 
	mutate(present_orthogroup = case_when(N > 0 ~ TRUE, N == 0 ~ FALSE)) %>% 
	dcast(HOG ~ genome, value.var = "present_orthogroup")

alpha_ecoli_presence_combomat <- alpha_ecoli_presence %>% make_comb_mat

##########################################################################################
# This function takes a list of genome names and returns a concatenated string 
# representing the presence or absence of each genome in the list. The string is 
# constructed by checking each genome in the input list against a reference 
# list of all genomes, and setting a "1" if the genome is present, or "0" if it 
# is absent. Reference list is obtained from the matrix: "alpha_ecoli_presence_combomat"

genome_presence <- function(genomes, all_genomes) {
	presence_vector <- rep(0, length(all_genomes))
	for (genome in genomes) {
		idx <- match(genome, all_genomes)
		if (!is.na(idx)) {
			presence_vector[idx] <- 1
		}
	}
	return(paste(presence_vector, collapse = ""))
}

# genome_presence("Zymomonas mobilis", alpha_ecoli_presence_combomat %>% rownames)
# genome_presence(c("Zymomonas mobilis", "Escherichia coli BW25113"), alpha_ecoli_presence_combomat %>% rownames)

all_genomes <- alpha_ecoli_presence_combomat %>% rownames

##########################################################################################

#All Alphas + E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bjUeidx5OYaE
all_present_idx <- genome_presence(
	all_genomes, 
	all_genomes)

all_present <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, all_present_idx), ]

#All Alphas - E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bzNJgxu0NyeQ 
all_alphas_only_idx <- genome_presence(
	all_genomes[!(all_genomes %in% ("Escherichia coli BW25113"))], 
	all_genomes)

all_alphas_only <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, all_alphas_only_idx), ]

#Zymomonas mobilis only - https://version-11-5.string-db.org/cgi/network?networkId=bDpAcDlAAZj
Zm_only <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, genome_presence("Zymomonas mobilis", all_genomes)), ]

##########################################################################################













# 
# 
# never_empty_alpha_ecoli <- alpha_ecoli_empty[
# 	extract_comb(
# 		alpha_ecoli_empty_combomat, 
# 		comb_name(alpha_ecoli_empty_combomat)[length(comb_name(alpha_ecoli_empty_combomat))]), ]
# 
# # orthos that are never missing in alphas
# core_alpha_ecoli_genome <- 
# 	orthos %>% 
# 	filter(HOG %in% (never_empty %>% pull(HOG))) %>%
# 	separate(genome, c("Genus", "Rest"), remove = FALSE, extra = "merge") %>%
# 	filter(Genus %in% c(
# 		"Escherichia",
# 		"Zymomonas",
# 		"Caulobacter",
# 		"Agrobacterium",
# 		"Brevundimonas",
# 		"Sphingomonas",
# 		"Rhodobacter",
# 		"Bradyrhizobium",
# 		"Rhodopseudomonas",
# 		"Novosphingobium"))
# 
# # fwrite(core_alpha_genome, "core_alpha_genome.tsv", sep = "\t")
# 
# ################################################################################
# 
# presence <- alpha_ecoli_orthos_wider %>% 
# 	pivot_longer(
# 		!c(genus, genome), 
# 		values_to = "N", 
# 		names_to = "HOG") %>%
# 	mutate(nonempty_orthogroup = case_when(N > 0 ~ TRUE, 1 == 1 ~ FALSE))
# 
# presence <- presence %>% 
# 	pivot_wider(id_cols = HOG, names_from = genome, values_from = nonempty_orthogroup)
# 
# presence.mixed.matrix <- presence  %>%
# 	make_comb_mat
# 
# single_species <- comb_name(presence.mixed.matrix) %>% 
# 	as_tibble %>% 
# 	mutate(comb_level = lengths(regmatches(value, gregexpr("1", value)))) %>% 
# 	mutate(index = 1:nrow(.)) %>% 
# 	filter(comb_level == 1 ) %>% 
# 	pull(value)
# 
# ################################################################################
# # BELOW THIS DOESN"T NEECSSARILY WORK
# 
# zymo_only <- presence[extract_comb(presence.mixed.matrix, comb_name(presence.mixed.matrix)[392]),]
# 
# zymo_only_genome <-
# 	orthos %>%
# 	filter(HOG %in% (zymo_only %>% pull(HOG))) %>%
# 	separate(genome, c("Genus", "Rest"), remove = FALSE, extra = "merge") %>%
# 	filter(Genus %in% c(
# 		"Zymomonas",
# 		"Caulobacter",
# 		"Agrobacterium",
# 		"Brevundimonas",
# 		"Sphingomonas",
# 		"Rhodobacter",
# 		"Bradyrhizobium",
# 		"Rhodopseudomonas",
# 		"Novosphingobium"))
# 
# fwrite(zymo_only_genome, "zymo_only_genome.tsv", sep = "\t")
# 
# z <- comb_name(presence.mixed.matrix) %>% 
# 	as_tibble %>% 
# 	mutate(comb_level = lengths(regmatches(value, gregexpr("1", value)))) %>% 
# 	mutate(index = 1:nrow(.)) %>% 
# 	filter(comb_level == 1 | comb_level == 9) %>% pull(index)
# 
# p <- UpSet(
# 	presence.mixed.matrix[, z], 
# 	top_annotation = HeatmapAnnotation(
# 		"Unique Genes" = anno_barplot(
# 			comb_size(presence.mixed.matrix[, z]),
# 			border = FALSE, 
# 			height = unit(10, "cm"),
# 			add_numbers = T),
# 		annotation_name_side = "left"),
# 	right_annotation = rowAnnotation(
# 		"Genome Size" = anno_barplot(
# 			set_size(presence.mixed.matrix[, z]),
# 			border = FALSE, 
# 			width = unit(10, "cm"),
# 			add_numbers = T)))
# 
# print(p)
# 
