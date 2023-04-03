
y_pred <- fread("y_pred_recalc2.tsv")

zmo_key2 <- zmo_key %>% left_join(y_pred %>% rename(spacer = variant)) %>% select(spacer, original, locus_tag, y_pred_calc, type)

zmo_key2[type == "perfect", original:= spacer]

melted_results2 <- zmo_key2 %>% inner_join(melted_results) 

results_mismatch <- melted_results2 %>% filter(type == "mismatch")
results_perfect <- melted_results2 %>% filter(type == "perfect")


results_mismatch <- results_mismatch %>% select(-type) %>% rename(LFC_mismatch = LFC, FDR_mismatch = FDR)
results_perfect <- results_perfect %>% select(-spacer, -y_pred_calc, -type) %>% rename(LFC_perfect = LFC, FDR_perfect = FDR)

