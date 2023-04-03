
y_pred <- fread("y_pred_recalc2.tsv")

zmo_key2 <- zmo_key %>% left_join(y_pred %>% rename(spacer = variant)) %>% select(spacer, original, locus_tag, y_pred_calc, type)

zmo_key2[type == "perfect", original:= spacer]

melted_results2 <- zmo_key2 %>% inner_join(melted_results) 

results_mismatch <- melted_results2 %>% filter(type == "mismatch")
results_perfect <- melted_results2 %>% filter(type == "perfect")


results_mismatch <- results_mismatch %>% select(-type) %>% rename(LFC_mismatch = LFC, FDR_mismatch = FDR)
results_perfect <- results_perfect %>% select(-spacer, -y_pred_calc, -type) %>% rename(LFC_perfect = LFC, FDR_perfect = FDR)

# Randomly select a locus_tag
random_locus_tag <- sample(unique(results_mismatch$locus_tag), 1)

# Filter results_mismatch based on the random_locus_tag and plot using ggplot2
results_mismatch %>% 
	filter(locus_tag == random_locus_tag) %>% 
	ggplot(aes(x = y_pred_calc, y = LFC_mismatch)) + 
	facet_wrap(~condition) + 
	geom_line() + ggtitle(random_locus_tag) + geom_point()

# Draw galaxy plots
results_perfect %>%
	inner_join(results_mismatch, by = join_by(original, locus_tag, gene_name, condition)) %>%
	ggplot(aes(x = LFC_perfect, y = LFC_mismatch, color = y_pred_calc)) +
	geom_point(size = 0.5, alpha = 0.5) +
	scale_color_viridis(option = "C") + facet_wrap(~condition)