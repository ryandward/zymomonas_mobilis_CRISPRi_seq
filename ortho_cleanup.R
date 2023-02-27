library(pacman)

p_load(data.table, tidyverse, ggrepel, hrbrthemes, ggallin, RColorBrewer, viridis, grid, clipr)

p_load_current_gh("jokergoo/ComplexHeatmap")

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
	mutate(orthologs = strsplit(as.character(orthologs), ", ")) %>% 
	unnest(cols = "orthologs") %>%
	data.table %>%
	rename(locus_tag = orthologs)

ecoli_genes <- fread("Escherichia_coli_BW25113.tsv")

################################################################################
# parse

orthos <- genome_locus_tags %>% full_join(orthos) %>%
	separate(genome, c("genus", "species", "rest"), remove = FALSE, extra = "merge") %>%
	mutate(genome = paste(genus, species, rest)) %>% 
	mutate(genome = gsub(" NA", "", genome)) %>% 
	data.table

orthos[is.na(HOG), HOG := locus_tag]

model_organisms <- c("Escherichia coli BW25113", "Caulobacter crescentus")

orthos$genome <- factor(orthos$genome)

genome_order <- c(model_organisms, sort(setdiff(levels(orthos$genome), model_organisms)))

orthos$genome <- factor(orthos$genome, levels = genome_order)

orthos_summary <- orthos %>% 
	group_by(HOG, OG, `Gene Tree Parent Clade`, genus, genome) %>% 
	summarise(N = n()) 

essential_orthos <- orthos %>% 
	inner_join(essentials %>% select(locus_tag), by = "locus_tag", multiple = "all") %>% 
	select(HOG) %>% unique %>% mutate(essential = TRUE)

orthos_summary <- orthos_summary %>% left_join(essential_orthos, multiple = "all") %>% data.table

orthos_summary[is.na(essential), essential := FALSE]

orthos_table <- orthos_summary %>% 
	data.table %>% 
	dcast(genome + genus ~ HOG, value.var = "N", fill = 0) %>%
	as_tibble()

alpha_ecoli_orthos_table <- 	
	orthos_table %>% 
	filter(
		genome %in% c(
			"Escherichia coli BW25113",
			"Caulobacter crescentus",
			"Agrobacterium fabrum",
			"Bradyrhizobium diazoefficiens",
			"Brevundimonas subvibrioides",
			"Novosphingobium aromaticivorans",
			"Rhodobacter sphaeroides",
			"Rhodopseudomonas palustris",
			"Sphingomonas wittichii",
			"Zymomonas mobilis"))
			

alpha_ecoli_orthos_table <- alpha_ecoli_orthos_table %>% 
	data.table %>% 
	melt(id.vars = c("genome", "genus"), variable.name = "HOG", value.name = "N") %>% 
	filter(N > 0) %>% 
	dcast(genome + genus ~ HOG, value.var = 'N', fill = 0) %>% 
	as_tibble()

# d <- 
# 	alpha_ecoli_orthos_table %>%
# 	select(-genome, -genus) %>%
# 	dist(method = "canberra")
# 
# d_fit <- cmdscale(d, eig = TRUE, k = 2)
# 
# d_scale <- d_fit$points %>% data.table
# 
# d_scale$genome <- alpha_ecoli_orthos_table$genome
# d_scale$genus <- alpha_ecoli_orthos_table$genus
# 
# d_scale %>% ggplot(aes(x = V1, y = V2)) + 
# 	geom_point(aes(fill = genus), size = 5, shape = 21, alpha = 0.5) +
# 	geom_text_repel(aes(label = genome)) +
# 	doc_theme


essential_orthos_summary <- orthos %>% 
	inner_join(essentials %>% select(locus_tag), by = "locus_tag", multiple = "all") %>% 
	group_by(genome, HOG) %>% tally(name = "essential_ortho_count")

##########################################################################################

genomes_presence <- alpha_ecoli_orthos_table %>% 
	data.table %>% 
	melt(id.vars = c("genus", "genome"), variable.name = "HOG", value.name = "N") %>% 
	mutate(present_orthogroup = case_when(N > 0 ~ TRUE, N == 0 ~ FALSE)) %>% 
	dcast(HOG ~ genome, value.var = "present_orthogroup")

genomes_presence_combo_distinct <- genomes_presence %>% make_comb_mat(mode = "distinct")
genomes_presence_combo_intersect <- genomes_presence %>% make_comb_mat(mode = "intersect")
genomes_presence_combo_union <- genomes_presence %>% make_comb_mat(mode = "union")

##########################################################################################
# This function takes a list of genome names and returns a concatenated string
# representing the presence or absence of each genome in the list. The string is
# constructed by checking each genome in the input list against a reference list
# of all genomes, and setting a "1" if the genome is present, or "0" if it is
# absent. Reference list is obtained from the matrix:
# "genomes_presence_combo_distinct"

extract_idx <- function(genomes, all_genomes) {
	presence_vector <- rep(0, length(all_genomes))
	for (genome in genomes) {
		idx <- match(genome, all_genomes)
		if (!is.na(idx)) {
			presence_vector[idx] <- 1
		}
	}
	return(paste(presence_vector, collapse = ""))
}

lookup_idx <- function(table, idx, combo) {
	return(table[extract_comb(combo, idx), ] %>% pull(HOG))
}

##########################################################################################

all_genomes <- genomes_presence_combo_distinct %>% rownames

##########################################################################################

#All Alphas + E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bjUeidx5OYaE
all_present_idx <- extract_idx(all_genomes, all_genomes)

all_present <- genomes_presence[extract_comb(
	genomes_presence_combo_distinct, all_present_idx), ]

#All Alphas - E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bzNJgxu0NyeQ 
all_alphas_only_idx <- extract_idx(all_genomes[!(all_genomes %in% ("Escherichia coli BW25113"))], all_genomes)

all_alphas_only <- genomes_presence[extract_comb(
	genomes_presence_combo_distinct, all_alphas_only_idx), ]

#Zymomonas mobilis only - https://version-11-5.string-db.org/cgi/network?networkId=bDpAcDlAAZj
Zm_only <- genomes_presence[extract_comb(
	genomes_presence_combo_distinct, extract_idx("Zymomonas mobilis", all_genomes)), ]

##########################################################################################
# create a table with info about each genome

genome_sets <- data.table()
for (i in all_genomes) {
	
	genome_sets <- orthos %>% 
		filter(genome == i) %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "All genes") %>%
		data.table %>% 
		rbind(genome_sets)

	genome_sets <- all_present %>% 
		select(HOG) %>% 
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "Present in all alphas and E. coli") %>% 
		data.table %>% 
		rbind(genome_sets)
	
	genome_sets <- all_alphas_only %>% 
		select(HOG) %>% 
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "Common only to all alphas") %>% 
		data.table %>% 
		rbind(genome_sets)
	
	genome_sets <-	genomes_presence[extract_comb(
		genomes_presence_combo_distinct, 
		extract_idx(i, all_genomes)), ] %>% 		
		select(HOG) %>%
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>%
		mutate(set = "Unique to genome") %>% 
		data.table %>% 
		rbind(genome_sets)
}

genome_sets$genome <- factor(genome_sets$genome, levels = genome_order)

genome_sets_stats <- genome_sets %>% 
	group_by(genome, set) %>% 
	summarise(size = n()) %>% 
	data.table %>% 
	dcast(genome ~ set, value.var = "size", fill = 0) %>%
	melt(id.vars = "genome", value.name = "size", variable.name = "set")

genome_sets_stats$genome <- factor(genome_sets_stats$genome, levels = genome_order)

##########################################################################################
# take info from upset plot, but make it simpler and plot

palette <- brewer.pal(8, "Paired")
my_colors <- palette[c(2, 4, 6, 8)]

genome_sets_stats_plot <- genome_sets_stats %>% 
	ggplot(aes(fill = set, y = `size`, x = genome)) + 
	geom_bar(position = "dodge", stat = "identity") +
	doc_theme +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(values = my_colors)

print(genome_sets_stats_plot)

##########################################################################################
# Hard questions, these are answered using "intersect" and "union" combo matrices.

# All genomes
All_idx <- extract_idx(all_genomes, all_genomes)

alphas_eco_shared <- lookup_idx(
	genomes_presence, 
	All_idx, 
	genomes_presence_combo_distinct) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

# Zymomonas only index
Zmo_idx <- extract_idx("Zymomonas mobilis", all_genomes)

Zmo_genome <- lookup_idx(
	genomes_presence, 
	Zmo_idx, 
	genomes_presence_combo_union) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

Eco_idx <- extract_idx("Escherichia coli BW25113", all_genomes)

Eco_genome <- lookup_idx(
	genomes_presence, 
	Eco_idx, 
	genomes_presence_combo_union) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

Alphas_idx <- extract_idx(all_genomes[!(all_genomes %in% ("Escherichia coli BW25113"))], all_genomes)

Alphas_common_not_Eco <- lookup_idx(
	genomes_presence, 
	Alphas_idx, 
	genomes_presence_combo_distinct) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)


# Now we don't care about E. coli. These HOGS are present everywhere in Alphas
# and may be in E. coli # Generate an index that corresponds to everything m
# (minus) E. coli
all_alphas_idx <- extract_idx(
	all_genomes[!(all_genomes %in% (c("Escherichia coli BW25113")))], 
	all_genomes)

# alpha_common_p_Zmo_m_all <- 
	
# Find HOGs in all alphas, and may be in E. coli
alpha_common <- lookup_idx(
	genomes_presence, 
	all_alphas_idx, 
	genomes_presence_combo_intersect) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

alpha_union <- lookup_idx(
	genomes_presence, 
	all_alphas_idx, 
	genomes_presence_combo_union) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

# Now we don't care about E. coli or Zymomonas
# Find index
alpha_not_Zmo_idx <- extract_idx(all_genomes[!(all_genomes %in% (c("Escherichia coli BW25113", "Zymomonas mobilis")))], all_genomes)

alpha_common_maybe_Zmo <- lookup_idx(
	genomes_presence, 
	alpha_not_Zmo_idx, 
	genomes_presence_combo_intersect) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

alpha_union_maybe_Zmo <- lookup_idx(
	genomes_presence, 
	alpha_not_Zmo_idx, 
	genomes_presence_combo_union) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

alpha_not_Zmo_not_Eco <- lookup_idx(
	genomes_presence, 
	alpha_not_Zmo_idx, 
	genomes_presence_combo_distinct) %>% 
	droplevels %>% 
	data.table(HOG = .) %>% 
	arrange(HOG)

### Genes shared by all Alphas and E. coli
## on E. coli
# https://version-11-5.string-db.org/cgi/network?networkId=bml1LlWvnUED
alphas_eco_shared %>% inner_join(orthos %>% filter(genome == "Escherichia coli BW25113"), multiple = "all") %>%
	select(locus_tag) %>% inner_join(ecoli_genes) %>% select(gene) %>%
	write_clip()

## on Caulobacter
# https://version-11-5.string-db.org/cgi/network?networkId=bS8ooauefqaO
alphas_eco_shared %>% inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>%
	select(locus_tag) %>% write_clip()

## on Zymomonas
# https://version-11-5.string-db.org/cgi/network?networkId=bBC7ZBZczODY
alphas_eco_shared %>% inner_join(orthos %>% filter(genome == "Zymomonas mobilis"), multiple = "all") %>%
	select(locus_tag) %>% 	mutate(locus_tag = gsub("ZMO1_", "", locus_tag)) %>% write_clip()




## Genes that may be in E. coli, but in ALL alphas including Zymomonas
# on Caulobacter
# https://version-11-5.string-db.org/cgi/network?networkId=bOEFxWYpUzvy
alpha_common_on_Ccr <- alpha_common %>%
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>%
	select(locus_tag)

# on Zymomonas
# https://version-11-5.string-db.org/cgi/network?networkId=bcrrXPZjXBlZ
alpha_common_on_Zmo <- alpha_common %>%
	inner_join(orthos %>% filter(genome == "Zymomonas mobilis"), multiple = "all") %>%
	select(locus_tag)

## Genes that may be in E. coli, may be in Zymomonas, but present in all other alphas
# https://version-11-5.string-db.org/cgi/network?networkId=b1UklJYoaHIx
alpha_common_maybe_Zmo_on_Ccr <- alpha_common_maybe_Zmo %>%
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>%
	select(locus_tag)

## Genes in all alphas except Zymomonas, maybe in E. coli
# https://version-11-5.string-db.org/cgi/network?networkId=buYCYUpvXGzW
alpha_common_not_Zmo_on_Ccr <- alpha_common_maybe_Zmo %>% 
	anti_join(alpha_common) %>% 
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>%
	select(locus_tag)

## Genes in Zymomonas, maybe E. coli, not other alphas
# https://version-11-5.string-db.org/cgi/network?networkId=b0tJaZaTX72r
Zmo_genome %>%
	anti_join(alpha_union_maybe_Zmo) %>%
	inner_join(orthos %>% filter(genome == "Zymomonas mobilis"), multiple = "all") %>%
	select(locus_tag) %>%
	mutate(locus_tag = gsub("ZMO1_", "", locus_tag)) %>%
	clipr::write_clip()

## Genes in alphas, maybe Zmo, not E. coli
# https://version-11-5.string-db.org/cgi/network?networkId=boJkAooCE4nu
alpha_common_maybe_Zmo %>% 
	anti_join(Eco_genome) %>% 
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>% 
	select(locus_tag) %>% clipr::write_clip()

## Genes in alphas, including Zmo, not E. coli
# on Caulobacter
# https://version-11-5.string-db.org/cgi/network?networkId=buOIfJpqAKIg
Alphas_common_not_Eco %>% 
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>% 
	select(locus_tag) %>% clipr::write_clip()

# on Zymomonas
# https://version-11-5.string-db.org/cgi/network?networkId=bMJ5kCiTbfrD
Alphas_common_not_Eco %>% 
	inner_join(orthos %>% filter(genome == "Zymomonas mobilis"), multiple = "all") %>% 
	select(locus_tag) %>%
	mutate(locus_tag = gsub("ZMO1_", "", locus_tag)) %>%
	clipr::write_clip()




