library(ggallin)
library(hrbrthemes)

doc_theme <- theme_ipsum(base_family = "Arial", caption_margin = 12, axis_title_size = 12, axis_col = "black")


volcano_plot <- function(median_melted_results, conditions, max_left = 5, max_right = 5, max_top = 5, max_sig = 1e-75){
	median_melted_results %>% 
		mutate(gene_name = case_when(
			gene_name == "." | is.na(gene_name)  ~ locus_tag,
			gene_name != "." ~ gene_name)) %>%
		filter(type == "perfect") %>%
		mutate(FDR = case_when(FDR < 1e-75 ~ FDR == min(FDR[FDR!= 0]), TRUE ~ FDR)) %>%
		mutate(Significance = case_when(
			FDR < 0.05 & abs(medLFC) >= 1 ~ "Significant",
			TRUE ~ "Not Significant")) %>%
		filter(condition %in% conditions) %>%
		group_by(condition) %>% 
		arrange(medLFC) %>% mutate(medLFC_ix = row_number()) %>%
		arrange(desc(medLFC)) %>% mutate(medLFC_desc_ix = row_number()) %>%
		arrange(FDR) %>% mutate(FDR_ix = row_number()) %>%
		mutate(FDR = case_when(FDR != 0 ~ FDR, FDR == 0 ~ min(FDR[FDR != 0]))) %>%
		ggplot(
			aes(x = medLFC,
					y = FDR)) +
		geom_point(aes(alpha = Significance), shape = 21) +
		geom_hline(yintercept = 0.05,
							 linetype = "dashed",
							 color = "#5A5A5A",
							 lwd = 0.5) +
		geom_vline(xintercept = -1,
							 linetype = "dashed",
							 color = "#5A5A5A", 
							 lwd = 0.5) +
		geom_vline(xintercept = 1,
							 linetype = "dashed",
							 color = "#5A5A5A", 
							 lwd = 0.5) +
		geom_vline(xintercept = 0,
							 linetype = "solid",
							 color = "#5A5A5A",
							 lwd = 0.5) +
		doc_theme +
		scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
		geom_label_repel(
			fill = alpha(c("white"),0.75),
			max.iter = 1000000000,
			data = . %>% filter(
				FDR < 0.05 & 
					(medLFC_ix <= max_left & medLFC < -0.5 | 
					 	(medLFC_desc_ix <= max_right & medLFC > 0.5 ) | 
					 	FDR_ix <= max_top)),
			aes(label = gene_name),
			force = 7.5,
			segment.size = 0.15,
			min.segment.length = 0,
			box.padding = 2,
			point.padding = .25,
			size = 2,
			max.overlaps = Inf,
			colour = "black") +
		facet_wrap(~ condition, scales = "free") +
		scale_alpha_manual(
			values = c(
				"Significant" = 0.75,
				"Not Significant" = 0.15))
	# +
	# 	scale_fill_manual(
	# 		values = c(
	# 			"Other" = "grey",
	# 			"Ribosome" = "#E31A1C",
	# 			"Ribosome+Other" = "#FB9A99",
	# 			"Ox Phos" = "#6A3D9A",
	# 			"Ox Phos+Other" = "#CAB2D6",
	# 			"LOS" = "#33A02C",
	# 			"LOS+Other" = "#B2DF8A",
	# 			"Cell Wall/PG" = "#FF7F00",
	# 			"Cell Wall/PG+Other" = "#FDBF6F",
	# 			"tRNA Ligase" = "#1F78B4",
	# 			"tRNA Ligase+Other" = "#A6CEE3")) +
	# 	guides(fill = guide_legend(
	# 		override.aes = list(shape = 21, size = 5))) + 
	# 	guides(alpha = guide_legend(
	# 		override.aes = list(shape = 21, size = 5, fill = "black"))) +
	# 	scale_alpha_manual(
	# 		values = c(
	# 			"Significant" = 0.75,
	# 			"Not Significant" = 0.15)) +
	# 	scale_size_area(breaks = c(1, 5, 10, 13))
}

