require('pacman');
p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)


RowVar <- function(x, ...) {
	
	rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
	
}

if (!exists("zmo_max_outlier")) {zmo_max_outlier = 0}
if (!exists("zmo_min_outlier")) {zmo_min_outlier = 0}


if (!exists("zmo_bed")) {
	
	zmo_bed <- fread("/home/ryandward/zmo_sgrna_analysis/CP023715.1.bed",
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
	
	zmo1 <- fread("/home/ryandward/zmo_sgrna_analysis/ZM1_run1_summary.tsv",
								col.names = c("spacer",
															"count",
															"condition"))
	
	zmo2 <- fread("/home/ryandward/zmo_sgrna_analysis/ZM1_run2_summary.tsv",
								col.names = c("spacer",
															"count",
															"condition"))
	
	zmo3 <- fread("/home/ryandward/zmo_sgrna_analysis/ZM2_run_summary.tsv",
								col.names = c("spacer",
															"count",
															"condition"))
	
	zmo <- rbind(zmo1, zmo2, zmo3)[, .(count = sum(count)), by = .(spacer, condition)]
	
	# setorderv(zmo, condition, spacer)
	
	rm(zmo1, zmo2, zmo3)
	
}

# if (!exists("zmo_dictionary")) {
#
# 	aba_dictionary <- fread("/home/ryandward/sgrna_analysis/Orthos/abau_dictionary.tsv")
#
# }

if (!exists("zmo_key")) {
	
	zmo_key <- fread("/home/ryandward/zmo_sgrna_analysis/Z._mobilis_CRISPRi_library_annotations.txt")
	
	zmo_key <- unique(zmo_key)
	
	zmo_key[, spacer := target]
	
	zmo_key[, target := NULL]
	
	zmo_key[, ZMO := gsub("ZMO1_", "", locus_tag)]
	
}

if (!exists("zmo_design")) {
	
	zmo_design <- fread("/home/ryandward/zmo_sgrna_analysis/ZMO1_experimental_design.tsv", na.strings = c("#N/A"))

	zmo_design <- zmo_design[!(condition %in% c("ZM1-23", "ZM1-34", "ZM1-35", "ZM2-2", "ZM2-9", "ZM2-10", "ZM2-1"))]
		
}

zmo <- zmo[condition %in% zmo_design$condition]

zmo_grid <- dcast(zmo, spacer ~ factor(condition, levels = unique(condition)), value.var = "count", fill = 1)

zmo_grid_matrix <- data.matrix(zmo_grid[, -c("spacer")])

row.names(zmo_grid_matrix) <- zmo_grid$spacer

# zmo_group <- factor(zmo_design[,  paste(timing, verbose, bio_rep, sep = "_")])
zmo_group <- factor(zmo_design[,  paste(verbose, timing, sep = "_")])

levels(zmo_group) <- c("T0", levels(zmo_group)[-1])

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

zmo_contrast <- makeContrasts(contrasts = c("aerobic_T1 - T0",
																						"aerobic_T2 - T0",
																						"aerobic_T3 - T0",
																						"anaerobic_T1 - T0",
																						"anaerobic_T2 - T0",
																						"anaerobic_T3 - T0"), 
															levels = zmo_permut)

results <- glmQLFTest(zmo_fit, contrast = zmo_contrast)

results <- topTags(results, n = Inf)
results <- data.table(results$table)
########################
results <- zmo_genome[results, on = .(spacer == genes)]

########################

plot(density(-log10(results[, FDR])), main = "Distribution of FDR, aerobic vs non-aerobic growth", xlab = "-log10 FD")


#####################

clipboard <- function(
	x , 
	sep = "\t" , 
	row.names = FALSE , 
	col.names = TRUE 
) {
	con <- pipe(
		"xclip -selection clipboard -i" , 
		open = "w"
	)
	
	write.table(
		x , 
		con , 
		sep = sep , 
		row.names = row.names , 
		col.names = col.names
	)
	close(con)}

RowVar <- function(x , ...) {
	rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

z <- estimateGLMCommonDisp(zmo_y, zmo_permut)
z <- estimateGLMTrendedDisp(z, zmo_permut, method="power")
z <- estimateGLMTagwiseDisp(z, zmo_permut)
plotBCV(z)
