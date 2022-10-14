require('pacman');
p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)


RowVar <- function(x, ...) {
	
	rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
	
}

if (!exists("zmo_max_outlier")) {zmo_max_outlier = 0}
if (!exists("zmo_min_outlier")) {zmo_min_outlier = 0}


if (!exists("zmo_bed")) {
	
	zmo_bed <- fread("CP023715.1.bed",
									 col.names = c("chromosome",
									 							"left",
									 							"right",
									 							"locus_tag",
									 							"gene_name",
									 							"strand",
									 							"coding",
									 							"completeness"))
}

if (!exists("zmo")) {
	
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
	
	zmo <- rbind(zmo1, zmo2, zmo3)[, .(count = sum(count)), by = .(spacer, condition)]
	
	# setorderv(zmo, condition, spacer)
	
	rm(zmo1, zmo2, zmo3)
	
}

# if (!exists("zmo_dictionary")) {
#
# 	zmo_dictionary <- fread("/home/ryandward/sgrna_analysis/Orthos/zmou_dictionary.tsv.gz")
#
# }

if (!exists("zmo_key")) {
	
	zmo_key <- fread("Z._mobilis_CRISPRi_library_annotations.txt")
	
	zmo_key <- unique(zmo_key)
	
	zmo_key[, spacer := target]
	
	zmo_key[, target := NULL]
	
	zmo_key[, ZMO := gsub("ZMO1_", "", locus_tag)]
	
}

zmo_design <- fread("ZMO1_experimental_design.tsv.gz", na.strings = c("#N/A"))

zmo <- zmo[condition %in% zmo_design$condition]

zmo_grid <- dcast(zmo, spacer ~ factor(condition, levels = unique(condition)), value.var = "count", fill = 1)

zmo_grid_matrix <- data.matrix(zmo_grid[, -c("spacer")])

row.names(zmo_grid_matrix) <- zmo_grid$spacer

# zmo_group <- factor(zmo_design[,  paste(timing, verbose, bio_rep, sep = "_")])
zmo_group <- factor(zmo_design[,  paste(verbose, timing, sep = "_")])

zmo_bio_reps <- factor(zmo_design$bio_rep)

zmo_y <- DGEList(counts = zmo_grid_matrix, group = zmo_group, genes = row.names(zmo_grid_matrix))

zmo_permut <-  model.matrix(~ 0 + zmo_group)
zmo_permut_br <-  model.matrix(~ 0 + zmo_bio_reps)

# zmo_permut <-  cbind(zmo_permut, zmo_permut_br)

colnames(zmo_permut) <- levels(zmo_group)
# colnames(zmo_permut) <- c(levels(zmo_group), paste("bio_rep", levels(zmo_bio_reps)[-1], sep = "_"))

zmo_keep <- filterByExpr(zmo_y, zmo_permut)

zmo_y <- zmo_y[zmo_keep, , keep.lib.sizes = FALSE]

zmo_y <- calcNormFactors(zmo_y)

zmo_y <- estimateDisp(zmo_y, zmo_permut)

plotBCV(zmo_y)

zmo_fit <- glmQLFit(zmo_y, zmo_permut, robust = TRUE)

plotQLDisp(zmo_fit)

zmo_CPM <- cpm(zmo_y, prior.count = 1)

colnames(zmo_CPM) <- factor(zmo_design[,  paste(verbose, timing, sep = "_")])

colors <- rep(c(alpha("black",0.25), alpha("green",0.5), alpha("magenta",0.5), alpha("red",0.5)), 2)

pch <- c(rep(19,4), rep(21,4))

plotMDS(zmo_y, col = colors[zmo_group], pch = pch[zmo_group], cex = 3)

legend("topleft", legend = levels(zmo_group), pch = pch, col = colors, ncol = 2)

title("Multidimensional Scaling: Guide Count by Aerobicity")

########################

zmo_genome <- zmo_bed[zmo_key[, .(spacer, type, locus_tag)], on = .(locus_tag)]

########################

contrast_levels <- c("aerobic_T1 - aerobic_T0",
										 "aerobic_T2 - aerobic_T0",
										 "aerobic_T3 - aerobic_T0",
										 "anaerobic_T1 - anaerobic_T0",
										 "anaerobic_T2 - anaerobic_T0",
										 "anaerobic_T3 - anaerobic_T0")

zmo_contrast <- makeContrasts(contrasts = contrast_levels, 
															levels = zmo_permut)

results <- glmQLFTest(zmo_fit, contrast = zmo_contrast)

results <- topTags(results, n = Inf)
results <- data.table(results$table)
########################
results <- zmo_genome[results, on = .(spacer == genes)]
#######################
contrast_levels <- gsub(" - ", "...", contrast_levels)
contrast_levels <- gsub("^", "logFC.", contrast_levels)

######################
melted_results <- 
	data.table::melt(
		results, 
		id.vars = c(
			"locus_tag",
			"gene_name",
			"type",
			"spacer",
			"logCPM",
			"F",
			"FDR",
			"PValue"),
		variable.name = "condition", 
		value.name = "logFC",
		measure.vars = contrast_levels)


melted_results[
	, logFC := melted_results[
		i  = type == "control",
		j  = .(ctrl_medLFC = median(logFC)),
		by = .(condition)][
			i  = .SD, # In this case, .SD is multiple in nature -- it refers to each of these sub-data.tables, one-at-a-time
			on = .(condition),
			j  = .(adj_medLFC = logFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

gene_level_results <- 
	melted_results[
		, 
		.(medLFC = median(logFC), 
			FDR = stouffer(FDR)$p), 
		by = .(locus_tag, gene_name, type, condition)]

gene_level_results[
	, medLFC := gene_level_results[
		i = type == "control",
		j = .(ctrl_medLFC = median(medLFC)),
		by = .(condition)][
			i  = .SD, # In this case, .SD is multiple in nature -- it refers to each of these sub-data.tables, one-at-a-time
			on = .(condition),
			j  = .(adj_medLFC = medLFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

liquid_results <- data.table::dcast(gene_level_results, locus_tag + gene_name + type + FDR ~ condition , value.var = "medLFC")

fwrite(liquid_results, "liquid_results.tsv.gz", sep = "\t")