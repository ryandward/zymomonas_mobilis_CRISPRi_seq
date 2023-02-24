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

orthos_table <- orthos_summary %>% 
	data.table %>% 
	data.table::dcast(genome + genus ~ HOG, value.var = "N", fill = 0) %>%
	as_tibble()

alpha_ecoli_orthos_table <- 	
	orthos_table %>% 
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

##########################################################################################

p_load(data.table, viridis, grid, tidyverse)
p_load_current_gh("jokergoo/ComplexHeatmap")

ecoli_genes <- fread("Escherichia_coli_BW25113.tsv")

alpha_ecoli_presence <- alpha_ecoli_orthos_table %>% 
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

idx_to_hog <- function(ortho_table, idx, combo_matrix) {
	return(alpha_ecoli_presence[extract_comb(alpha_ecoli_presence_combomat, idx), ] %>% pull(HOG))
}

# genome_presence("Zymomonas mobilis", alpha_ecoli_presence_combomat %>% rownames)
# genome_presence(c("Zymomonas mobilis", "Escherichia coli BW25113"), alpha_ecoli_presence_combomat %>% rownames)

all_genomes <- alpha_ecoli_presence_combomat %>% rownames

##########################################################################################

#All Alphas + E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bjUeidx5OYaE
all_present <- genome_presence(all_genomes, all_genomes)

all_present <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, all_present_idx), ]

#All Alphas - E. coli: https://version-11-5.string-db.org/cgi/network?networkId=bzNJgxu0NyeQ 
all_alphas_only_idx <- genome_presence(all_genomes[!(all_genomes %in% ("Escherichia coli BW25113"))], all_genomes)

all_alphas_only <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, all_alphas_only_idx), ]

#Zymomonas mobilis only - https://version-11-5.string-db.org/cgi/network?networkId=bDpAcDlAAZj
Zm_only <- alpha_ecoli_presence[extract_comb(
	alpha_ecoli_presence_combomat, genome_presence("Zymomonas mobilis", all_genomes)), ]

##########################################################################################
# create a table with info about each genome

genome_sets <- data.table()
for (i in all_genomes) {
	
	genome_sets <- orthos %>% 
		filter(genome == i) %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "All genes") %>%
		data.table %>% rbind(genome_sets)

	genome_sets <- all_present %>% 
		select(HOG) %>% 
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "Present in all alphas and E. coli") %>% 
		data.table %>% rbind(genome_sets)
	
	genome_sets <- all_alphas_only %>% 
		select(HOG) %>% 
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>% 
		mutate(set = "Common only to all alphas") %>% 
		data.table %>% rbind(genome_sets)
	
	genome_sets <-	alpha_ecoli_presence[extract_comb(
		alpha_ecoli_presence_combomat, 
		genome_presence(i, all_genomes)), ] %>% 		
		select(HOG) %>%
		inner_join(orthos %>% filter(genome == i), multiple = "all") %>% 
		select(genome, genus, locus_tag) %>%
		mutate(set = "Unique to genome") %>% 
		data.table %>% rbind(genome_sets)
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
# Hard questions



# Now we don't care about E. coli. These HOGS are present everywhere in Alphas and may be in E. coli
## Find index
alphas_p_ec_idx <- genome_presence(all_genomes, all_genomes)
alphas_m_ec_idx <- genome_presence(all_genomes[!(all_genomes %in% (c("Escherichia coli BW25113")))], all_genomes)

# Find HOGs pom (plus or minus) E. coli
alphas_pom_ec <- c(
	idx_to_hog(alpha_ecoli_presence, alphas_p_ec_idx, alpha_ecoli_presence_combomat), 
	idx_to_hog(alpha_ecoli_presence, alphas_m_ec_idx, alpha_ecoli_presence_combomat)) %>% 
	droplevels %>% data.table(HOG = .) %>%
	arrange(HOG)

# Now we don't care about E. coli or Zymomonas
## Find index
alphas_m_zm_idx <- genome_presence(all_genomes[!(all_genomes %in% (c("Zymomonas mobilis")))], all_genomes)
alphas_m_ec_m_ec_idx <- genome_presence(all_genomes[!(all_genomes %in% (c("Escherichia coli BW25113", "Zymomonas mobilis")))], all_genomes)

# Find HOGs pom (plus or minus) E. coli, Zymomonas
alphas_pom_ec_pom_zm <- c(
	idx_to_hog(alpha_ecoli_presence, alphas_p_ec_idx, alpha_ecoli_presence_combomat),
	idx_to_hog(alpha_ecoli_presence, alphas_m_ec_idx, alpha_ecoli_presence_combomat),
	idx_to_hog(alpha_ecoli_presence, alphas_m_zm_idx, alpha_ecoli_presence_combomat), 
	idx_to_hog(alpha_ecoli_presence, alphas_m_ec_m_ec_idx, alpha_ecoli_presence_combomat)) %>% 
	droplevels %>% data.table(HOG = .) %>%
	arrange(HOG)

# Find the difference and map it onto Caulobacter 
alphas_pom_ec_m_zm_mapped_cc <- alphas_pom_ec_pom_zm %>% 
	anti_join(alphas_pom_ec) %>% 
	inner_join(alpha_ecoli_presence) %>% 
	inner_join(orthos %>% filter(genome == "Caulobacter crescentus"), multiple = "all") %>% select(locus_tag)

