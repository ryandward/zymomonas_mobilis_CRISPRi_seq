library(data.table)
library(tidyverse)
library(ggrepel)
library(hrbrthemes)
library(ggallin)

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")

################################################################################
# read and wrangle

essentials <- fread(
	"Essentials.tsv", 
	header = F, 
	col.names = c("genome", "locus_tag"))

genome_locus_tags <- fread(
	"genome_locus_tags.tsv.gz", 
	header = F, 
	col.names = c("genome", "locus_tag"))

orthos <- fread("N0.tsv")

orthos <- melt(
	orthos, 
	id.vars = c("HOG", "OG", "Gene Tree Parent Clade"), 
	variable.name = "genome", 
	value.name = "orthologs") %>%
	filter(!genome == "Caulobacter_vibrioides") %>%
	mutate(orthologs = strsplit(as.character(orthologs), ", ")) %>% 
	unnest(cols = "orthologs") %>%
	data.table %>%
	rename(locus_tag = orthologs)

################################################################################
# parse

orthos <- genome_locus_tags %>% full_join(orthos) %>%
	filter(!genome == "Caulobacter_vibrioides") %>%
	separate(genome, c("genus", "species", "rest"), remove = FALSE, extra = "merge") %>%
	mutate(genome = paste(genus, species, rest)) %>% 
	mutate(genome = gsub(" NA", "", genome))

orthos[is.na(HOG), HOG := locus_tag]

model_organisms <- c("Escherichia coli BW25113", "Caulobacter crescentus")

genome_order <- c(model_organisms, sort(setdiff(levels(orthos$genome), model_organisms)))

orthos$genome <- factor(orthos$genome, levels = genome_order)


orthos_summary <- orthos %>% 
	group_by(HOG, OG, `Gene Tree Parent Clade`, genus, genome) %>% 
	summarise(N = n()) 

orthos_wider <- orthos_summary %>% 
	data.table %>% 
	data.table::dcast(genome + genus ~ HOG, value.var = "N", fill = 0) %>%
	as_tibble()

alpha_ecoli_orthos_wider <- 	
	orthos_wider %>% 
	filter(
		genome == "Escherichia coli BW25113" |
		genus %in% c(
		"Zymomonas",
		"Caulobacter",
		"Agrobacterium",
		"Brevundimonas",
		"Sphingomonas",
		"Rhodobacter",
		"Bradyrhizobium",
		"Rhodopseudomonas",
		"Novosphingobium"))

alpha_ecoli_orthos_wider <- alpha_ecoli_orthos_wider %>% 
	data.table %>% 
	melt(id.vars = c("genome", "genus"), variable.name = "HOG", value.name = "N") %>% 
	filter(N > 0) %>% 
	dcast(genome + genus ~ HOG, value.var = 'N', fill = 0) %>% 
	as_tibble()

d <- 
	alpha_ecoli_orthos_wider %>%
	select(-genome, -genus) %>%
	dist(method = "canberra")

d_fit <- cmdscale(d, eig = TRUE, k = 2)

d_scale <- d_fit$points %>% data.table

d_scale$genome <- alpha_ecoli_orthos_wider$genome
d_scale$genus <- alpha_ecoli_orthos_wider$genus

d_scale %>% ggplot(aes(x = V1, y = V2)) + 
	geom_point(aes(fill = genus), size = 5, shape = 21, alpha = 0.5) +
	geom_text_repel(aes(label = genome)) +
	doc_theme

