library(pacman)

p_load(data.table, viridis, grid, tidyverse)
p_load_current_gh("jokergoo/ComplexHeatmap")

digitsum <- function(x) sum(floor(x / 10^(0:(nchar(x) - 1))) %% 10)

source("ortho_wrangler.R")


alpha_ecoli_empty <- alpha_ecoli_orthos_wider %>% 
	pivot_longer(
		!c(genus, genome), 
		values_to = "N", 
		names_to = "HOG") %>% mutate(empty_orthogroup = case_when(N == 0 ~ TRUE, 1 == 1 ~ FALSE))

alpha_ecoli_empty <- alpha_ecoli_empty %>% 
	pivot_wider(id_cols = HOG, names_from = genome, values_from = empty_orthogroup)

alpha_ecoli_mixed.matrix <- alpha_ecoli_empty %>%
	make_comb_mat

p <- UpSet(
	alpha_ecoli_mixed.matrix, 
	top_annotation = HeatmapAnnotation(
		"Intersection" = anno_barplot(
			comb_size(alpha_ecoli_mixed.matrix),
			border = FALSE, 
			height = unit(8, "cm"),
			add_numbers = T)))

alpha_ecoli_orthogroups_per_genome <- alpha_ecoli_orthos_wider %>% 
	pivot_longer(
		!c(genus, genome), 
		values_to = "N", 
		names_to = "HOG") %>% 
	filter(N > 0) %>%
	pivot_wider(
		id_cols = HOG, 
		names_from = genome, 
		values_from = N, 
		values_fill = 0)


# hogs that are never missing in alphas
never_empty <- empty[extract_comb(mixed.matrix, comb_name(mixed.matrix)[length(comb_name(mixed.matrix))]), ]

# orthos that are never missing in alphas
core_alpha_ecoli_genome <- 
	orthos %>% 
	filter(HOG %in% (never_empty %>% pull(HOG))) %>%
	separate(genome, c("Genus", "Rest"), remove = FALSE, extra = "merge") %>%
	filter(Genus %in% c(
		"Escherichia",
		"Zymomonas",
		"Caulobacter",
		"Agrobacterium",
		"Brevundimonas",
		"Sphingomonas",
		"Rhodobacter",
		"Bradyrhizobium",
		"Rhodopseudomonas",
		"Novosphingobium"))

# fwrite(core_alpha_genome, "core_alpha_genome.tsv", sep = "\t")

################################################################################

presence <- alpha_ecoli_orthos_wider %>% 
	pivot_longer(
		!c(genus, genome), 
		values_to = "N", 
		names_to = "HOG") %>%
	mutate(nonempty_orthogroup = case_when(N > 0 ~ TRUE, 1 == 1 ~ FALSE))

presence <- presence %>% 
	pivot_wider(id_cols = HOG, names_from = genome, values_from = nonempty_orthogroup)

presence.mixed.matrix <- presence  %>%
	make_comb_mat

single_species <- comb_name(presence.mixed.matrix) %>% 
	as_tibble %>% 
	mutate(comb_level = lengths(regmatches(value, gregexpr("1", value)))) %>% 
	mutate(index = 1:nrow(.)) %>% 
	filter(comb_level == 1 ) %>% 
	pull(value)

################################################################################
# BELOW THIS DOESN"T NEECSSARILY WORK

zymo_only <- presence[extract_comb(presence.mixed.matrix, comb_name(presence.mixed.matrix)[392]),]

zymo_only_genome <-
	orthos %>%
	filter(HOG %in% (zymo_only %>% pull(HOG))) %>%
	separate(genome, c("Genus", "Rest"), remove = FALSE, extra = "merge") %>%
	filter(Genus %in% c(
		"Zymomonas",
		"Caulobacter",
		"Agrobacterium",
		"Brevundimonas",
		"Sphingomonas",
		"Rhodobacter",
		"Bradyrhizobium",
		"Rhodopseudomonas",
		"Novosphingobium"))

fwrite(zymo_only_genome, "zymo_only_genome.tsv", sep = "\t")

z <- comb_name(presence.mixed.matrix) %>% 
	as_tibble %>% 
	mutate(comb_level = lengths(regmatches(value, gregexpr("1", value)))) %>% 
	mutate(index = 1:nrow(.)) %>% 
	filter(comb_level == 1 | comb_level == 9) %>% pull(index)

p <- UpSet(
	presence.mixed.matrix[, z], 
	top_annotation = HeatmapAnnotation(
		"Unique Genes" = anno_barplot(
			comb_size(presence.mixed.matrix[, z]),
			border = FALSE, 
			height = unit(10, "cm"),
			add_numbers = T),
		annotation_name_side = "left"),
	right_annotation = rowAnnotation(
		"Genome Size" = anno_barplot(
			set_size(presence.mixed.matrix[, z]),
			border = FALSE, 
			width = unit(10, "cm"),
			add_numbers = T)))

print(p)

