source("ortho_wrangler.R")

library(data.table)
library(tidyverse)
library(ggrepel)
library(hrbrthemes)
library(ggallin)

zymo_orthos <- fread("ZMO_N0.tsv")

zymo_orthos <- melt(
	zymo_orthos, 
	id.vars = c("HOG", "OG", "Gene Tree Parent Clade"), 
	variable.name = "genome", 
	value.name = "orthologs") %>%
	mutate(orthologs = strsplit(as.character(orthologs), ", ")) %>% 
	unnest(cols = "orthologs") %>%
	data.table %>%
	rename(locus_tag = orthologs)

zymo_orthos_summary <- zymo_orthos %>% 
	group_by(HOG, OG, `Gene Tree Parent Clade`, genome) %>% 
	summarise(N = n()) %>% 
	separate(genome, c("Genus", "Species", "Strain"), remove = FALSE, extra = "merge")

orthos_summary <- orthos_summary %>% mutate(genome = paste(Genus, Rest))

ZM4_uniques <- zymo_only %>% 
	select(HOG) %>% 
	left_join(orthos) %>% 
	filter(genome == "Zymomonas_mobilis") %>% 
	select(locus_tag) %>% 
	inner_join(zymo_orthos) %>% 
	select(HOG, OG) %>% 
	inner_join(zymo_orthos) %>% 
	filter(genome == "Zymomonas_mobilis_ZM4") %>% 
	select(locus_tag) 

fwrite(ZM4_uniques, "ZM4_uniques.tsv", sep = "\t")

ZMATCC10988_uniques <- zymo_only %>% 
	select(HOG) %>% 
	left_join(orthos) %>% 
	filter(genome == "Zymomonas_mobilis") %>% 
	select(locus_tag) %>% 
	inner_join(zymo_orthos) %>% 
	select(HOG, OG) %>% 
	inner_join(zymo_orthos) %>% 
	filter(genome == "Zymomonas_mobilis_ATCC10988") %>% 
	select(locus_tag) 

fwrite(ZMATCC10988_uniques, "ZMATCC10988_uniques.tsv", sep = "\t")
