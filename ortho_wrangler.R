library(data.table)
library(tidyverse)
library(ggrepel)
library(hrbrthemes)
library(ggallin)

essentials <- fread(
	"Essentials.tsv", 
	header = F, 
	col.names = c("genome", "locus_tag"))

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

orthos <- fread("N0.tsv")

orthos <- melt(
	orthos, 
	id.vars = c("HOG", "OG", "Gene Tree Parent Clade"), 
	variable.name = "genome", 
	value.name = "orthologs") %>%
	filter(!genome %like% "vibrioides") %>%
	mutate(orthologs = strsplit(as.character(orthologs), ", ")) %>% 
	unnest(cols = "orthologs") %>%
	data.table %>%
	rename(locus_tag = orthologs)

orthos_summary <- orthos %>% 
	group_by(HOG, OG, `Gene Tree Parent Clade`, genome) %>% 
	summarise(N = n()) %>% 
	separate(genome, c("Genus", "Rest"), remove = FALSE, extra = "merge")

orthos_summary <- orthos_summary %>% mutate(genome = paste(Genus, Rest))

orthos_wider <- orthos_summary %>% 
	pivot_wider(id_cols = c(genome, Genus), names_from = HOG, values_from = N, values_fill = 0)

# alpha_orthos_wider %>% 
# 	pivot_longer(
# 		!c(Genus, genome), 
# 		values_to = "N", 
# 		names_to = "HOG") %>% 
# 	filter(N > 0) %>%
# 	pivot_wider(
# 		id_cols = c(genome, Genus), 
# 		names_from = HOG, 
# 		values_from = N, 
# 		values_fill = 0)

alpha_orthos_wider <- 	
	orthos_wider %>% 
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

alpha_orthos_wider <- alpha_orthos_wider %>% 
	pivot_longer(
		!c(Genus, genome), 
		values_to = "N", 
		names_to = "HOG") %>% 
	filter(N > 0) %>%
	mutate(N = ifelse(N > 5, 5, N)) %>% pivot_wider(
		id_cols = c(genome, Genus), 
		names_from = HOG, 
		values_from = N, 
		values_fill = 0)
	

d <- 
	alpha_orthos_wider %>%
	select(-genome, -Genus) %>%
	dist(method = "canberra")

d_fit <- cmdscale(
	d,
	eig = TRUE,
	k = 2)

d_scale <- d_fit$points %>% data.table

d_scale$genome <- alpha_orthos_wider$genome
d_scale$Genus <- alpha_orthos_wider$Genus

d_scale %>% ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(fill = Genus), size = 7, shape = 21, alpha = 0.25) +
	geom_text_repel(aes(label = genome)) +
	doc_theme


#########

single_copy_essentials <- essentials %>% 
	inner_join(orthos) %>% 
	group_by(HOG, OG, genome) %>% 
	tally %>% 
	pivot_wider(
		id_cols = c(HOG, OG), 
		names_from = genome, 
		values_from = n, 
		values_fill = 0) %>% 
	pivot_longer(!c(OG, HOG)) %>% 
	filter(value == 1) %>% 
	group_by(HOG, OG) %>% 
	tally %>% 
	filter(n >= 6) %>% 
	select(HOG, OG) %>% 
	inner_join(orthos) %>% 
	filter(genome == "Zymomonas_mobilis")

############################################################################
# 
# median_melted_results %>% 
# 	filter(
# 		condition %in% c(
# 			"aerobic_T3 - aerobic_T0", 
# 			"anaerobic_T3 - anaerobic_T0")) %>% 
# 	ggplot(aes(x = medLFC)) + geom_density() + doc_theme
# 
# 
# median_melted_results %>% 
# 	inner_join(single_copy_essentials) %>% 
# 	filter(
# 		condition %in% c(
# 			"aerobic_T3 - aerobic_T0", 
# 			"anaerobic_T3 - anaerobic_T0")) %>% 
# 	ggplot(aes(x = medLFC)) + geom_density() + doc_theme


# 
# median_melted_results_all <-
# 	median_melted_results %>% 
# 	filter(
# 		type == "perfect",
# 		condition %in% c(
# 			"aerobic_T3 - aerobic_T0", 
# 			"anaerobic_T3 - anaerobic_T0")) %>%
# 	mutate(set = "All")
# 
# median_melted_results_orthoparsed <-
# 	median_melted_results %>% 
# 	inner_join(single_copy_essentials) %>% 
# 	filter(
# 		type == "perfect",
# 		condition %in% c(
# 			"aerobic_T3 - aerobic_T0", 
# 			"anaerobic_T3 - anaerobic_T0")) %>%
# 	mutate(set = "Orthos = 6/6 Alphas")
# 
# setwise_results <- 
# 	median_melted_results_orthoparsed %>% 
# 	rbind(median_melted_results_all, fill = T)
# 
# 
# setwise_results %>%
# 	arrange(medLFC) %>%
# 	ggplot(aes(x = medLFC, y = FDR)) + 
# 	geom_point(aes(colour = condition), alpha = 0.5) + 
# 	doc_theme + 
# 	scale_y_continuous(
# 		trans = scales::reverse_trans() %of% scales::log10_trans()) +
# 	facet_grid(facets = c("condition", "set"))

