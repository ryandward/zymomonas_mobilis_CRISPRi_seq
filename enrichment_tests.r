require("pacman")

p_load(
  data.table,
  scales,
  edgeR,
  statmod,
  poolr,
  pheatmap,
  svglite,
  ggplot2,
  ggrepel,
  Rtsne,
  pracma,
  colourpicker,
  RColorBrewer,
  vegan,
  tidyverse,
  magrittr,
  ggtext,
  ggforce
)

p_load_current_gh("shabbychef/ggallin")

# aggregate the reads from different runs
zmo1 <- fread(
  "ZM1_run1_summary.tsv.gz",
  col.names = c(
    "spacer",
    "count",
    "condition"
  )
)

zmo2 <- fread(
  "ZM1_run2_summary.tsv.gz",
  col.names = c(
    "spacer",
    "count",
    "condition"
  )
)

zmo3 <- fread(
  "ZM2_run_summary.tsv.gz",
  col.names = c(
    "spacer",
    "count",
    "condition"
  )
)

zmo <- rbind(zmo1, zmo2, zmo3)[, .(count = sum(count)), by = .(spacer, condition)]
targets <- fread("zmo_targets.tsv")
targets <- targets %>% mutate(locus_tag = str_replace(locus_tag, "ZMO1_", ""))

duplicated_guides <- targets %>%
  filter(note %like% "sites") %>%
  select(spacer)

zmo <- zmo %>% filter(!spacer %in% duplicated_guides$spacer)

exp_verbose <- fread(
  "ZMO1_experimental_design_comprehensive.tsv.gz",
  na.strings = c("#N/A")
)

# define the experimental design space to only take into consideration "liquids"
exp_verbose <-
  exp_verbose[
    experiment == "liquid" & verbose != "glycerol"
  ]

# keep only the counts that are in the experimental design space
zmo <- zmo[
  condition %in% exp_verbose$condition
]

# Explicitly reorder 'zmo' based on the 'condition' order in 'exp_verbose'
zmo <- exp_verbose %>%
  select(condition) %>%
  left_join(zmo, by = "condition")


count_grid <- data.table::dcast(
  zmo,
  spacer ~ factor(condition, levels = unique(condition)),
  value.var = "count",
  fill = 0
)

count_grid_matrix <- data.matrix(
  count_grid[
    ,
    -c("spacer")
  ]
)

row.names(count_grid_matrix) <- count_grid$spacer

exp_group <- factor(
  exp_verbose[, paste(verbose, timing, sep = "_")],
  levels = unique(
    exp_verbose[, paste(verbose, timing, sep = "_")]
  )
)

data_design <- model.matrix(~ 0 + exp_group)
colnames(data_design) <- levels(exp_group)
rownames(data_design) <- exp_verbose$condition

dge <- DGEList(
  counts = count_grid_matrix,
  group = exp_group,
  genes = row.names(count_grid_matrix),
  design = data_design
)

dge <- calcNormFactors(dge)
dge <- estimateGLMRobustDisp(y = dge, design = data_design)

contrast_levels <- c(
  "aerobic_T1 - aerobic_T0",
  "aerobic_T2 - aerobic_T0",
  "aerobic_T3 - aerobic_T0",
  "anaerobic_T1 - anaerobic_T0",
  "anaerobic_T2 - anaerobic_T0",
  "anaerobic_T3 - anaerobic_T0",
  "aerobic_T1 - anaerobic_T1",
  "aerobic_T2 - anaerobic_T2",
  "aerobic_T3 - anaerobic_T3"
)

zmo_contrast <- makeContrasts(
  contrasts = contrast_levels,
  levels = data_design
)

dge$design <- data_design

colnames(dge$design) <- gsub("[^[:alnum:]]", "_", colnames(dge$design))
colnames(dge$design) <- gsub("___", " - ", colnames(dge$design))

dge$samples$group <- gsub("[^[:alnum:]]", "_", dge$samples$group)
dge$samples$group <- gsub("___", " - ", dge$samples$group)

contrast_levels <- gsub("[^[:alnum:]]", "_", contrast_levels)
contrast_levels <- gsub("___", " - ", contrast_levels)

contrast_list <- setNames(as.list(contrast_levels), contrast_levels)
contrast_list$levels <- dge$design

# Then use these contrast levels in the makeContrasts function
contrasts <- do.call(makeContrasts, contrast_list)



all_string <- fread("~/R/zymomonas_mobilis_CRISPRi_seq/264203.protein.enrichment.terms.v12.0.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()

targets <- targets


gene_groups <- all_string %>%
  # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
  group_by(category, term, description) %>%
  summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
  mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))


term_stats <- gene_groups %>%
  unnest(locus_tag) %>%
  inner_join(targets %>%
    select(locus_tag) %>% unique()) %>%
  group_by(term, gene_count, description) %>%
  summarize(genes_targeted = n())

complete_terms <- term_stats %>%
  filter(gene_count == genes_targeted)

# only perform enrichments where all genes are available
# gene_groups <- complete_terms %>% inner_join(gene_groups)

repeated_gene_groups <- gene_groups %>%
  group_by(locus_tag) %>%
  mutate(times_listed = n()) %>%
  arrange(locus_tag) %>%
  ungroup()


# pick the best annotation for each locus_tag_group, i.e., highest in term, and the lowest in the category_rank
ranked_annotations <- repeated_gene_groups %>%
  group_by(locus_tag_group, category) %>%
  arrange(versionsort::ver_sort(term)) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(category_rank = case_when(
    category == "Biological Process (Gene Ontology)" ~ 1,
    category == "Molecular Function (Gene Ontology)" ~ 2,
    category == "Cellular Component (Gene Ontology)" ~ 3,
    category == "Protein Domains and Features (InterPro)" ~ 4,
    category == "Protein Domains (SMART)" ~ 5,
    category == "Protein Domains (Pfam)" ~ 6,
    category == "Annotated Keywords (UniProt)" ~ 7,
    category == "Reactome Pathways" ~ 8,
    category == "Subcellular localization (COMPARTMENTS)" ~ 9,
    category == "Local Network Cluster (STRING)" ~ 10,
    TRUE ~ NA_integer_
  )) %>%
  group_by(locus_tag_group) %>%
  filter(category_rank == min(category_rank))

enrichments <- ranked_annotations %>%
  ungroup() %>%
  distinct(locus_tag_group, .keep_all = TRUE) %>%
  select(-locus_tag_group) %>%
  unnest(locus_tag) %>%
  inner_join(term_stats)


# Get the unique terms
unique_terms <- unique(enrichments$term)

target_spacers_for_terms <- term_stats %>%
  inner_join(enrichments, relationship = "many-to-many") %>%
  inner_join(targets, relationship = "many-to-many")


#########################################################################################

# Split the spacer column by term
sets_to_locus_tags <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
sets_to_locus_tags_indices <- lapply(sets_to_locus_tags, function(locus_tags) which(rownames(dge) %in% locus_tags))

v <- voomWithQualityWeights(dge, dge$design, plot = TRUE)

# filter out guides that do not have bona-fide targets, but keep non-targeting guides
# i.e. guides that have no target, or guides that have a target in the same (wrong) direction
# overlap of spacer to gene should be maximal, i.e. 20

edge_targets <- dge$counts %>% data.table(keep.rownames = "spacer")

v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  filter(spacer %in% edge_targets$spacer) %>%
  select(spacer) %>%
  left_join(
    targets %>% select(spacer, y_pred, target, mismatches) %>% unique()
  )

# Assign the weight to the guides based on y_pred to be between 1 and 100
v_targets[y_pred == "None", y_pred := NA_integer_]

v_targets$y_pred <- as.numeric(v_targets$y_pred)

v_targets[is.na(target) | target == "None", weight := min(v_targets$y_pred, na.rm = TRUE)]

v_targets[spacer == target, weight := max(v_targets$y_pred, na.rm = TRUE)]

v_targets[mismatches >= 1, weight := y_pred]

v_targets$weight <- rescale(as.numeric(v_targets$weight), to = c(1, 100))

v_targets[is.na(weight), weight := 1]

# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = sets_to_locus_tags_indices,
    design = dge$design,
    weights = v_targets$weight,
    # inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)

all_sets <- all_sets %>%
  inner_join(enrichments) %>%
  inner_join(targets %>% inner_join(v_targets %>% select(-y_pred))) %>%
  group_by(contrast, term, description) %>%
  nest(locus_tags = locus_tag) %>%
  group_by(locus_tags, contrast) %>%
  mutate(missing_genes = gene_count - genes_targeted) %>%
  arrange(FDR, missing_genes) %>%
  ungroup() %>%
  rename(guide_count = NGenes) %>%
  select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count) %>%
  unique() %>%
  data.table()

contrast_assignments <- contrasts %>%
  data.table(keep.rownames = "group") %>%
  melt(
    id.vars = "group",
    variable.name = "contrast",
    value.name = "assignment"
  ) %>%
  filter(assignment != 0)

group_assignments <- dge$design %>%
  data.table(keep.rownames = "sample") %>%
  melt(id.vars = "sample", variable.name = "group") %>%
  filter(value != 0) %>%
  select(-value)

original_data <- dge$counts %>%
  data.table(keep.rownames = "spacer") %>%
  melt(
    value.name = "count",
    id.vars = "spacer",
    variable.name = "sample"
  )

annotated_data <- dge$samples %>%
  data.table(keep.rownames = "sample") %>%
  inner_join(original_data) %>%
  group_by(sample) %>%
  mutate(cpm = 1e6 * count / sum(count))


### you could create a function out of this

library(grDevices)

# create a color palette that gradually changes from the start color to the end color

# Define the start and end colors for each condition
aerobic_start_color <- rgb(0.7, 0.7, 1)  # light blue
aerobic_end_color <- rgb(0, 0, 1)  # dark blue
anaerobic_start_color <- rgb(1, 0.7, 0.7)  # light red
anaerobic_end_color <- rgb(1, 0, 0)  # dark red

# Create color palettes that gradually change from the start color to the end color
aerobic_palette <- colorRampPalette(c(aerobic_start_color, aerobic_end_color))
anaerobic_palette <- colorRampPalette(c(anaerobic_start_color, anaerobic_end_color))

# Define the conditions
aerobic_conditions <- c("aerobic_T0", "aerobic_T1", "aerobic_T2", "aerobic_T3")
anaerobic_conditions <- c("anaerobic_T0", "anaerobic_T1", "anaerobic_T2", "anaerobic_T3")

# Assign a color to each condition
aerobic_colors <- setNames(aerobic_palette(length(aerobic_conditions)), aerobic_conditions)
anaerobic_colors <- setNames(anaerobic_palette(length(anaerobic_conditions)), anaerobic_conditions)

# Combine the color assignments for both conditions
colors <- c(aerobic_colors, anaerobic_colors)

####

this_term <- "GO:0046219"

title <- term_stats %>%
  filter(term == this_term) %>%
  pull(description)

title <- paste(title, " (", this_term, ")", sep = "")

enrichment_plot <- contrast_assignments %>%
  inner_join(
    group_assignments,
    relationship = "many-to-many"
  ) %>%
  inner_join(
    all_sets %>%
      filter(term == this_term) %>%
      inner_join(contrast_assignments) %>%
      mutate(contrast = factor(contrast, levels = unique(contrast)))
  ) %>%
  mutate(label = paste(contrast, paste("FDR:", signif(FDR, 3)), paste(Direction, paste(paste(genes_targeted, gene_count, sep = "/"), " genes present", sep = "")), sep = "\n")) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  inner_join(enrichments) %>%
  inner_join(targets) %>%
  inner_join(annotated_data) %>%
  inner_join(v_targets %>% select(-y_pred)) %>%
  arrange(assignment) %>%
  mutate(group = factor(group, levels = unique(group))) %>%
  ggplot(aes(x = as.character(assignment), y = cpm)) +
  geom_sina(aes(weight = as.numeric(weight), size = weight, color = group), scale = "width") +
  geom_violin(aes(weight = as.numeric(weight)), alpha = 0.35, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width", lwd = 1.25) +
  geom_tile(aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))), width = Inf, height = Inf, fill = "light grey") +
  scale_alpha_manual(values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE) +
  scale_size(range = c(0.1, 3)) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(10^(0:5)),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  facet_wrap(~label) +
  ggtitle(title) +
  scale_color_manual(values = colors) +

  theme_minimal()


plot(enrichment_plot)













