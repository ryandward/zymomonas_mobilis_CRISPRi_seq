require('pacman');
p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)

zmo_bed <- fread("CP023715.1.bed",
								 col.names = c("chromosome",
								 							"left",
								 							"right",
								 							"locus_tag",
								 							"gene_name",
								 							"strand",
								 							"coding",
								 							"completeness"))

zmo1 <- fread("ZM1_run1_summary.tsv.gz",
							col.names = c("spacer",
														"count",
														"condition"))

zmo2 <- fread("ZM1_run2_summary.tsv.gz",
							col.names = c("spacer",
														"count",
														"condition"))

zmo3 <- fread("ZM2_run_summary.tsv.gz",
							col.names = c("spacer",
														"count",
														"condition"))

zmo <- rbind(zmo1, zmo2, zmo3)[, .(count = sum(count)), 
															 by = .(spacer, condition)]

zmo_key <- fread("zmo_key.tsv.gz")

zmo_design <- fread("ZMO1_experimental_design_comprehensive.tsv.gz", 
										na.strings = c("#N/A"))

# define the experimental design space to only take into consideration "plates"
zmo_design <- zmo_design[experiment == "plate" & verbose != "glycerol" ]

# keep only the counts that are in the experimental design space
zmo <- zmo[condition %in% zmo_design$condition]

# convert single column into a table 
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html
zmo_grid <- data.table::dcast(
	zmo, 
	spacer ~ factor(condition, levels = unique(condition)),
	value.var = "count", 
	fill = 0)

zmo_grid_matrix <- data.matrix(zmo_grid[, -c("spacer")])
row.names(zmo_grid_matrix) <- zmo_grid$spacer

zmo_group <- factor(zmo_design[,  paste(verbose, timing, sep = "_")])

zmo_permut <-  model.matrix(~ 0 + zmo_group)
colnames(zmo_permut) <- levels(zmo_group)

zmo_y <- DGEList(counts = zmo_grid_matrix, 
								 group = zmo_group, 
								 genes = row.names(zmo_grid_matrix))

zmo_keep <- filterByExpr(
	y = zmo_y, 
	design = zmo_permut,
	group = zmo_group)

zmo_y <- zmo_y[zmo_keep, , keep.lib.sizes = FALSE]
zmo_y <- calcNormFactors(zmo_y)
zmo_y <- estimateDisp(zmo_y, zmo_permut)

plotBCV(zmo_y)

zmo_fit <- glmQLFit(zmo_y, zmo_permut, robust = TRUE)

plotQLDisp(zmo_fit)

zmo_CPM <- cpm(zmo_y, prior.count = 1)

colnames(zmo_CPM) <- factor(zmo_design[,  paste(verbose, timing, sep = "_")])

colors <- rep(
	c(alpha("black", 0.25),
		alpha("red", 0.5)), 2)

pch <- c(
	rep(19, 2), 
	rep(21, 2))

plotMDS(zmo_y, 
				col = colors[zmo_group], 
				pch = pch[zmo_group], 
				cex = 3)

legend("bottomright", 
			 legend = levels(zmo_group), 
			 pch = pch, 
			 col = colors, 
			 ncol = 2)

title("MDS: Guide Count by Aerobicity (Plates)")

########################

zmo_genome <- zmo_bed[zmo_key[, .(spacer, type, locus_tag)], on = .(locus_tag)]

########################

contrast_levels <- c("aerobic_T1 - aerobic_T0",
										 "anaerobic_T1 - anaerobic_T0")

zmo_contrast <- makeContrasts(contrasts = contrast_levels, 
															levels = zmo_permut)

########################

results_FDR <- zmo_key[, .(genes = unique(spacer))]
results_LFC <- zmo_key[, .(genes = unique(spacer))]

########################

for (i in 1:ncol(zmo_contrast)){
	
	results <- glmQLFTest(zmo_fit, contrast = zmo_contrast[,i])
	results <- topTags(results, n = Inf)
	results <- data.table(results$table)
	
	print(paste("Processing results for", contrast_levels[i], "..."))
	
	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
	setnames(results_FDR, "FDR", contrast_levels[i])
	
	results_LFC <- results[, .(genes, logFC)][results_LFC, on = .(genes)]
	setnames(results_LFC, "logFC", contrast_levels[i])
	
}

########################

results_FDR <- zmo_genome[results_FDR, on = .(spacer == genes)]
results_LFC <- zmo_genome[results_LFC, on = .(spacer == genes)]

########################

melted_results_FDR <- 
	data.table::melt(
		results_FDR, 
		id.vars = c(
			"locus_tag",
			"gene_name",
			"type",
			"spacer"),
		variable.name = "condition", 
		value.name = "FDR",
		measure.vars = contrast_levels)

melted_results_LFC <- 
	data.table::melt(
		results_LFC, 
		id.vars = c(
			"locus_tag",
			"gene_name",
			"type",
			"spacer"),
		variable.name = "condition", 
		value.name = "LFC",
		measure.vars = contrast_levels)

melted_results_LFC[
	, LFC := melted_results_LFC[
		i  = type == "control",
		j  = .(ctrl_medLFC = median(LFC, na.rm = TRUE)),
		by = .(condition)][
			i  = .SD, # In this case, .SD is multiple in nature -- it refers to each of these sub-data.tables, one-at-a-time
			on = .(condition),
			j  = .(adj_medLFC = LFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

melted_results <- 
	melted_results_LFC[
		melted_results_FDR, 
		on = .(locus_tag, gene_name, type, spacer, condition)]

melted_results <- melted_results[!is.na(FDR)]

median_melted_results <- 
	melted_results[
		, 
		.(medLFC = median(LFC), 
			FDR = stouffer(FDR)$p), 
		by = .(locus_tag, gene_name, type, condition)]

median_melted_results[
	, medLFC := median_melted_results[
		i = type == "control",
		j = .(ctrl_medLFC = median(medLFC)),
		by = .(condition)][
			i  = .SD, # In this case, .SD is multiple in nature -- it refers to each of these sub-data.tables, one-at-a-time
			on = .(condition),
			j  = .(adj_medLFC = medLFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

plate_results_LFC <- data.table::dcast(melted_results, locus_tag + gene_name + type + spacer ~ condition , value.var = "LFC")
plate_results_FDR <- data.table::dcast(melted_results, locus_tag + gene_name + type + spacer ~ condition , value.var = "FDR")

fwrite(plate_results_LFC, "plate_results_LFC.tsv", sep = "\t")
fwrite(plate_results_FDR, "plate_results_FDR.tsv", sep = "\t")
