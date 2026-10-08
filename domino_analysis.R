#!/usr/bin/env Rscript
# ==============================================================================
# DOMINO - Microbial continuity across an olive farm-to-food continuum
#
# Reproducibility script for the statistical analyses and figures of:
#   Traina C., Antonielli L., Holcik L., et al. "Microbial continuity across an
#   olive farm-to-food continuum reveals stable reservoirs, ecological filtering
#   and microbiome transmission".
#
# It contains only the analyses that appear in the paper and reproduces the
# published figures and tables.
#
# ------------------------------------------------------------------------------
# HOW TO RUN
#   From the project root (the folder containing `bac/`, `fun/` and the
#   metadata file):
#       Rscript -e 'renv::restore()'      # once, installs the pinned packages
#       Rscript domino_analysis.R
#   All outputs are written to `output/figures/` and `output/tables/`.
#
# INPUTS (relative to the project root)
#   bac_meta_full.tsv                     sample metadata (Area, Tree,
#                                         Compartment, Year, Time)
#   bac/03.dada2_output/asv_table.txt     16S ASV table (ASVs x samples + taxonomy)
#   bac/04.stats/02.filtered_reads.txt    16S read counts per processing step
#   fun/03.dada2_output/asv_table.txt     ITS ASV table
#   fun/04.stats/02.filtered_reads.txt    ITS read counts per processing step
#
# PAPER ITEM -> SCRIPT SECTION
#   Figure 2   (A, B)  Alluvial ASV flow 2022 -> 2023 -> fermentation   Sec. 4
#   Figure 3A          Fermentation traceability (2022 / 2023 origins)  Sec. 5.1
#   Figure 3B          Traceability over fermentation time, per tree    Sec. 5.2
#   Figure 3C          Taxonomy of origin ∩ BrinedOlive ASVs            Sec. 5.3
#   Figure 4A          Alpha diversity 2022 vs 2023                     Sec. 6.1
#   Figure 4B          PCoA (Bray-Curtis), tree compartments            Sec. 6.2
#   Figure 4C          Cross-season ASV persistence                     Sec. 6.3
#   Figure 4D          Taxonomy of 2022 ∩ 2023 ASVs                     Sec. 6.4
#   Figure 5A          Richness / Simpson over fermentation             Sec. 7.1
#   Figure 5B          ASV turnover between consecutive timepoints     Sec. 7.2
#   Figure 5C          Dominant-taxon trajectories (start -> end)       Sec. 7.3
#   Supp. Table 1      Area comparison 2022 (alpha + PERMANOVA)         Sec. 8
#   Supp. Table 2      Cross-season alpha (A) + PERMANOVA (B)           Sec. 6.1, 6.2
#   Supp. Table 3      Cross-season persistence binomial GLMM           Sec. 6.3
#   Supp. Table 4      Bacteria vs fungi alpha during fermentation      Sec. 7.1
#   (Figure 1 is a map; Supplementary Figure 1 is weather data - not produced here.)
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(ggalluvial)
  library(vegan)
  library(rtk)
  library(decontam)
  library(emmeans)
  library(lme4)
})
# multcomp is only used through multcomp::cld() and is deliberately not
# attached, because it attaches MASS, which masks dplyr::select().
# cld() on emmeans objects also needs multcompView, which emmeans only suggests;
# naming it here makes renv record it in renv.lock.
stopifnot(requireNamespace("multcompView", quietly = TRUE))

# --- Parameters ---------------------------------------------------------------
# `seed` is passed to every rtk() call: with rtk's default (seed = 0) each call
# draws a new random seed, so rarefied results would change between runs.
seed        <- 42
threads_num <- max(1, min(10, parallel::detectCores() - 1))

min_reads_sample   <- 100   # samples with <= this many reads are removed
alpha_rare_depth   <- 500   # rarefaction depth for alpha diversity
alpha_rare_repeats <- 100   # rarefaction iterations for alpha diversity (Area A)
area_rare_repeats  <- 20    # rarefaction iterations for the 2022 area comparison
beta_min_reads     <- 1000  # samples below this are dropped before beta diversity
beta_rare_repeats  <- 100   # rarefaction iterations for beta diversity
pa_min_reads       <- 500   # min reads for presence/absence analyses (Fig 3A, 4C, 5B)
pa_rare_repeats    <- 20    # rarefaction iterations for presence/absence analyses
other_threshold    <- 2     # taxa below this % are lumped into "Other"
n_perm             <- 999   # PERMANOVA permutations

tree_compartments <- c("Bulksoil", "Rhizosphere", "Bark", "Leaves", "Olive")
compartments      <- c(tree_compartments, "BrinedOlive")
kingdom_names     <- c(bac = "Bacteria", fun = "Fungi")

# --- Colour palettes (as used in the paper) -----------------------------------
color_palette_compartment <- c(
  BrinedOlive = "#E69F00", Olive = "#2E7D32", Leaves = "#0097A7",
  Bark = "#8C510A", Rhizosphere = "#D98C8C", Bulksoil = "#660033"
)
color_palette_tree    <- c("1" = "#E69F00", "2" = "#009E73", "3" = "#CC79A7")
color_palette_season  <- c("2022" = "#5B0E7D", "2023" = "#E69F00")
color_palette_kingdom <- c(Bacteria = "#1f78b4", Fungi = "#e31a1c")
color_palette_change  <- c(Gained = "#4daf4a", Lost = "#e41a1c",
                           Maintained = "#377eb8", Persistent = "#377eb8")

# --- Output -------------------------------------------------------------------
fig_dir   <- file.path("output", "figures")
table_dir <- file.path("output", "tables")
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".svg")), plot, device = grDevices::svg,
         width = width, height = height, units = "in")
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300)
  message("Saved figure: ", name)
}

save_table <- function(df, name) {
  write.csv(df, file.path(table_dir, paste0(name, ".csv")), row.names = FALSE)
  message("Saved table:  ", name)
}

# Common ggplot theme of the panel figures
theme_panel <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text       = element_text(face = "bold", size = 18),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 16),
      axis.text.y      = element_text(size = 16),
      axis.title       = element_text(size = 18),
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 16),
      panel.grid.minor = element_blank()
    )
}

set.seed(seed)


# ==============================================================================
# 2. DATA IMPORT AND QUALITY CONTROL
# ==============================================================================

read_tsv_table <- function(path) {
  read.table(path, header = TRUE, row.names = 1, sep = "\t",
             comment.char = "", check.names = FALSE)
}

# --- Metadata (Area A, trees 1-3) ---------------------------------------------
# 2023 olives at T0 are relabelled "Olive" (they are the pre-brine fruits);
# 2022 samples get Time = 0.
meta_full <- read_tsv_table("bac_meta_full.tsv")

meta <- meta_full %>%
  filter(Compartment != "OliveTE", Area == "A", Tree %in% c("1", "2", "3")) %>%
  mutate(
    Compartment = ifelse(Time == 0 & Compartment == "BrinedOlive", "Olive", Compartment),
    Time        = ifelse(is.na(Time) & Year == "2022", 0, Time),
    Compartment = ifelse(Compartment == "BulkSoil", "Bulksoil", Compartment),
    Area        = factor(Area),
    Tree        = factor(Tree),
    Year        = factor(Year),
    Compartment = factor(Compartment, levels = compartments)
  )

# --- ASV tables and read-count statistics -------------------------------------
asv_raw <- list(
  bac = read_tsv_table("bac/03.dada2_output/asv_table.txt"),
  fun = read_tsv_table("fun/03.dada2_output/asv_table.txt")
)
read_stats <- list(
  bac = read_tsv_table("bac/04.stats/02.filtered_reads.txt"),
  fun = read_tsv_table("fun/04.stats/02.filtered_reads.txt")
)

# Keep the study samples plus negative controls; drop ASVs that become empty
subset_biom <- function(biom, meta) {
  keep_samples <- c(intersect(rownames(meta), colnames(biom)),
                    grep("neg", colnames(biom), ignore.case = TRUE, value = TRUE))
  biom_sub <- biom[, intersect(c(keep_samples, "taxonomy"), colnames(biom)), drop = FALSE]
  biom_sub[rowSums(biom_sub[, colnames(biom_sub) != "taxonomy"]) > 0, , drop = FALSE]
}

# Remove non-target sequences
#   Bacteria: keep Bacteria/Archaea, drop chloroplasts and mitochondria
#   Fungi:    keep Fungi, strip "NA" ranks from the taxonomy string
select_target_taxa <- function(biom, kingdom) {
  if (kingdom == "bac") {
    biom <- biom[grepl("k__Archaea|k__Bacteria", biom$taxonomy, ignore.case = TRUE), ]
    biom <- biom[!grepl("chloroplast|mitochondria", biom$taxonomy, ignore.case = TRUE), ]
  } else {
    biom <- biom[grepl("k__Fungi", biom$taxonomy, ignore.case = TRUE), ]
    biom$taxonomy <- gsub(";NA|; NA", "", biom$taxonomy)
  }
  biom
}

samples_x_asvs <- function(biom) as.data.frame(t(biom[, colnames(biom) != "taxonomy"]))

# Decontamination (decontam, prevalence method, threshold 0.5, negative
# extraction controls) followed by filtering of controls, ASVs with <= 1 read
# and samples with <= `min_reads_sample` reads.
decontaminate_and_filter <- function(read_counts, biom, tab, kingdom_label) {
  # Samples with 0 reads at any processing step are treated as failed
  decon_map <- read_counts %>%
    rownames_to_column("sample") %>%
    mutate(across(where(is.numeric), ~ na_if(., 0))) %>%
    drop_na() %>%
    mutate(is_neg = grepl("NEG", sample))

  biom <- biom[, c(decon_map$sample, "taxonomy")]
  tab  <- tab[decon_map$sample, ]

  contam <- isContaminant(as.matrix(tab), method = "prevalence",
                          neg = decon_map$is_neg, threshold = 0.5, normalize = TRUE)
  message(kingdom_label, ": ", sum(contam$contaminant), " contaminant ASVs removed")
  biom <- biom[!contam$contaminant, ]
  tab  <- tab[, !contam$contaminant]

  tab <- tab[!grepl("NEG|POS", rownames(tab)), ]
  tab <- tab[, colSums(tab) > 1]
  tab <- tab[rowSums(tab) > min_reads_sample, ]

  list(tab = tab, biom = biom[colnames(tab), c(rownames(tab), "taxonomy")])
}

biom_sel  <- list()   # target taxa, all study samples (used for taxonomy lookups)
filt_tab  <- list()   # filtered count tables (samples x ASVs)
filt_biom <- list()   # filtered ASV tables with taxonomy (ASVs x samples)
filt_meta <- list()   # metadata matching filt_tab

for (k in names(kingdom_names)) {
  biom_sub <- subset_biom(asv_raw[[k]], meta)
  stats    <- read_stats[[k]][colnames(biom_sub)[colnames(biom_sub) != "taxonomy"], ]

  biom_sel[[k]] <- select_target_taxa(biom_sub, k)
  tab_sel       <- samples_x_asvs(biom_sel[[k]])

  # Supplementary read-count summary per sample
  read_counts <- cbind(stats, read_counts = rowSums(tab_sel))
  write.table(read_counts, file.path(table_dir, paste0(kingdom_names[[k]], "_read_counts.tsv")),
              sep = "\t", col.names = NA)

  clean <- decontaminate_and_filter(read_counts, biom_sel[[k]], tab_sel, kingdom_names[[k]])
  filt_tab[[k]]  <- clean$tab
  filt_biom[[k]] <- clean$biom
  filt_meta[[k]] <- meta[rownames(clean$tab), ]
}


# ==============================================================================
# 3. SHARED HELPERS
# ==============================================================================

# --- Taxonomy -----------------------------------------------------------------
uninformative <- c("", "NA", "uncultured", "unidentified", "Incertae_sedis", "Incertae sedis")

# Rank label from a "k__;p__;...;g__" string; NA if absent or uninformative
get_tax_label <- function(taxonomy_string, level) {
  if (is.na(taxonomy_string)) return(NA_character_)
  prefix <- switch(level, genus = "g__", family = "f__", order = "o__",
                   class = "c__", phylum = "p__", stop("Invalid taxonomic level"))
  parts <- trimws(unlist(strsplit(taxonomy_string, ";")))
  val   <- gsub(prefix, "", grep(paste0("^", prefix), parts, value = TRUE))
  if (length(val) == 0 || val %in% uninformative) return(NA_character_)
  val
}

asv_tax_labels <- function(asvs, biom, level) {
  tax <- setNames(biom$taxonomy, rownames(biom))
  sapply(tax[asvs], get_tax_label, level = level)
}

# Label with fallback to a higher rank (used for Fig. 5C)
#   Fungi:    genus > family > order;   Bacteria: family > order > class
get_fallback_label <- function(tax_strings, kingdom) {
  sapply(tax_strings, function(x) {
    if (is.na(x)) return("Unclassified")
    parts <- unlist(strsplit(x, ";"))
    rank  <- function(prefix) gsub(prefix, "", grep(prefix, parts, value = TRUE))
    if (kingdom == "Fungi") {
      g <- rank("g__"); f <- rank("f__"); o <- rank("o__")
      if (length(g) > 0 && !g %in% c("unidentified", "Incertae_sedis", "")) return(g)
      if (length(f) > 0 && !f %in% c("unidentified", ""))                   return(paste0(f, " (Fam)"))
      if (length(o) > 0 && o != "unidentified")                             return(paste0(o, " (Ord)"))
      return("Unclassified Fungi")
    } else {
      f <- rank("f__"); o <- rank("o__"); cl <- rank("c__")
      if (length(f) > 0 && !f %in% c("uncultured", ""))  return(f)
      if (length(o) > 0 && !o %in% c("uncultured", ""))  return(paste0(o, " (Order)"))
      if (length(cl) > 0)                                return(paste0(cl, " (Class)"))
      return("Unclassified")
    }
  })
}

# Categorical palette with "Other" in grey
build_palette <- function(taxa, base_cols) {
  named <- sort(unique(taxa[taxa != "Other"]))
  pal   <- if (length(named) > length(base_cols)) colorRampPalette(base_cols)(length(named))
           else base_cols[seq_along(named)]
  c(setNames(pal, named), Other = "#CCCCCC")
}
base_cols_set1 <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999",
  "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#1B9E77",
  "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#003C30","#01665E",
  "#35978F","#80CDC1","#C7EAE5","#543005","#8C510A","#BF812D")
base_cols_tab20 <- c(
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2","#7F7F7F",
  "#BCBD22","#17BECF","#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5","#F7B6D2",
  "#C49C94","#DBDB8D","#9EDAE5","#393B79")

# --- Presence / rarefaction ---------------------------------------------------
# ASVs with > 0 reads in at least one of `samples`
get_present_asvs <- function(tab, samples) {
  samples <- intersect(samples, rownames(tab))
  if (length(samples) == 0) return(character(0))
  colnames(tab)[colSums(tab[samples, , drop = FALSE]) > 0]
}

# Drop samples below `min_reads`, then rarefy all remaining samples to the
# smallest library size (`repeats` iterations). Returns the list of rarefied
# matrices and the metadata of the kept samples.
rarefy_presence <- function(tab, meta, min_reads, repeats) {
  keep <- rownames(tab)[rowSums(tab) >= min_reads]
  dropped <- setdiff(rownames(tab), keep)
  if (length(dropped) > 0)
    message("  Dropping ", length(dropped), " samples < ", min_reads, " reads: ",
            paste(dropped, collapse = ", "))
  tab   <- tab[keep, ]
  tab   <- tab[, colSums(tab) > 0, drop = FALSE]
  depth <- min(rowSums(tab))
  message("  Rarefying ", nrow(tab), " samples to ", depth, " reads (", repeats, " iterations)")
  rar <- rtk(tab, depth = depth, repeats = repeats, ReturnMatrix = repeats,
             margin = 1, threads = threads_num, verbose = FALSE, seed = seed)
  list(mats = rar$raremat, meta = meta[keep, ], depth = depth)
}

# Tidy an adonis2 table
tidy_adonis2 <- function(adon, kingdom, analysis, model) {
  df <- as.data.frame(adon)
  data.frame(Kingdom = kingdom, Analysis = analysis, Model = model, Term = rownames(df),
             df[, c("Df", "SumOfSqs", "R2", "F", "Pr(>F)")],
             check.names = FALSE, row.names = NULL)
}


# ==============================================================================
# 4. FIGURE 2 - ALLUVIAL ASV FLOW: 2022 -> 2023 -> FERMENTATION
# ==============================================================================
# Per compartment (pooled across trees, no rarefaction), ASVs are classified as
#   Persisted (2022 & 2023) / Lost (2022 only) / New (2023 only),
# and persisted + new ASVs are intersected with the pooled BrinedOlive (T > 0)
# ASV set. Ribbon weights:
#   lost -> 2022 reads; persisted but not fermented -> mean 2022/2023 reads;
#   new but not fermented -> 2023 reads; reaching fermentation -> fermentation reads.
# Families < 2 % in every compartment are grouped as "Other".

compute_3stage_flow <- function(tab, meta, biom, kingdom) {
  asv_family <- asv_tax_labels(rownames(biom), biom, "family")

  ferm_samples <- rownames(meta)[meta$Compartment == "BrinedOlive" & meta$Time > 0]
  counts_ferm  <- colSums(tab[ferm_samples, , drop = FALSE])
  asvs_ferm    <- names(counts_ferm)[counts_ferm > 0]

  family_counts <- function(asvs, weights, fate, comp) {
    fam   <- asv_family[asvs]
    valid <- !is.na(fam)
    if (sum(valid) == 0) return(NULL)
    data.frame(Family = fam[valid], Count = weights[asvs[valid]]) %>%
      group_by(Family) %>%
      summarise(Count = sum(Count), .groups = "drop") %>%
      mutate(Fate = fate, Compartment = comp)
  }

  flows <- list()
  for (comp in tree_compartments) {
    s22 <- rownames(meta)[meta$Compartment == comp & meta$Year == "2022"]
    s23 <- rownames(meta)[meta$Compartment == comp & meta$Year == "2023"]
    if (length(s22) == 0 || length(s23) == 0) next

    counts_22 <- colSums(tab[s22, , drop = FALSE])
    counts_23 <- colSums(tab[s23, , drop = FALSE])
    asvs_22   <- names(counts_22)[counts_22 > 0]
    asvs_23   <- names(counts_23)[counts_23 > 0]
    persisted <- intersect(asvs_22, asvs_23)
    new       <- setdiff(asvs_23, asvs_22)

    flows[[comp]] <- bind_rows(
      family_counts(intersect(persisted, asvs_ferm), counts_ferm,                 "Persisted_to_Ferm",  comp),
      family_counts(setdiff(persisted, asvs_ferm),   (counts_22 + counts_23) / 2, "Persisted_not_Ferm", comp),
      family_counts(setdiff(asvs_22, asvs_23),       counts_22,                   "Lost",               comp),
      family_counts(intersect(new, asvs_ferm),       counts_ferm,                 "Gained_to_Ferm",     comp),
      family_counts(setdiff(new, asvs_ferm),         counts_23,                   "Gained_not_Ferm",    comp)
    )
  }
  bind_rows(flows) %>% mutate(Kingdom = kingdom)
}

fam_flow <- bind_rows(lapply(names(kingdom_names), function(k)
  compute_3stage_flow(filt_tab[[k]], filt_meta[[k]], filt_biom[[k]], kingdom_names[[k]])))

# Lump families that never reach `other_threshold` % within any compartment
major_families <- fam_flow %>%
  group_by(Kingdom, Compartment, Family) %>%
  summarise(Reads = sum(Count), .groups = "drop") %>%
  group_by(Kingdom, Compartment) %>%
  mutate(Pct = Reads / sum(Reads) * 100) %>%
  group_by(Kingdom, Family) %>%
  summarise(MaxPct = max(Pct), .groups = "drop") %>%
  filter(MaxPct >= other_threshold)

fam_flow <- fam_flow %>%
  left_join(major_families %>% transmute(Kingdom, Family, Major = TRUE), by = c("Kingdom", "Family")) %>%
  mutate(Family = ifelse(is.na(Major), "Other", Family)) %>%
  group_by(Kingdom, Compartment, Fate, Family) %>%
  summarise(Count = sum(Count), .groups = "drop")

# Long format: one row per flow and axis (2022, 2023, Fermented)
x_levels <- c("2022", "2023", "Fermented")
flow_strata <- function(fate, comp) {
  switch(fate,
         Persisted_to_Ferm  = c(comp,  comp,   "Fermented"),
         Persisted_not_Ferm = c(comp,  comp,   "Lost"),
         Lost               = c(comp,  "Lost", "Lost"),
         Gained_to_Ferm     = c("New", comp,   "Fermented"),
         Gained_not_Ferm    = c("New", comp,   "Lost"))
}

alluv_long <- bind_rows(lapply(seq_len(nrow(fam_flow)), function(i) {
  r <- fam_flow[i, ]
  data.frame(Kingdom = r$Kingdom, FlowID = i, Family = r$Family,
             x = factor(x_levels, levels = x_levels),
             stratum = flow_strata(r$Fate, r$Compartment), Count = r$Count)
})) %>%
  mutate(FlowID  = factor(FlowID),
         stratum = factor(stratum, levels = c(tree_compartments, "Fermented", "Lost", "New")))

for (king in kingdom_names) {
  d <- alluv_long %>% filter(Kingdom == king)
  fam_colors <- build_palette(as.character(d$Family), base_cols_set1)
  d$Family   <- factor(d$Family, levels = c(sort(setdiff(unique(d$Family), "Other")), "Other"))

  p <- ggplot(d, aes(x = x, stratum = stratum, alluvium = FlowID, y = Count, fill = Family)) +
    geom_stratum(width = 0.30, color = "gray30", linewidth = 0.3, fill = "white") +
    geom_flow(width = 0.30, alpha = 0.50) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5, fontface = "bold") +
    scale_fill_manual(values = fam_colors, name = "Family") +
    guides(fill = guide_legend(ncol = 1)) +
    labs(title = paste0(king, " — ASV Flow: Environment → Fermentation"),
         x = NULL, y = "Weighted Read abundance") +
    theme_bw(base_size = 12) +
    theme(legend.position  = "right",
          legend.text      = element_text(size = 12, face = "italic"),
          legend.key.size  = unit(0.4, "cm"),
          panel.grid.minor = element_blank())

  panel <- if (king == "Bacteria") "Fig2A" else "Fig2B"
  save_fig(p, paste0(panel, "_alluvial_", tolower(king)), width = 16, height = 9)
}

# Numbers quoted in Results 3.1.1
save_table(fam_flow %>%
             group_by(Kingdom, Compartment, Fate) %>%
             summarise(Abundance = sum(Count), .groups = "drop") %>%
             group_by(Kingdom, Compartment) %>%
             mutate(Percent = round(Abundance / sum(Abundance) * 100, 1)) %>%
             arrange(Kingdom, Compartment, desc(Percent)),
           "Fig2_alluvial_fate_summary")
save_table(fam_flow %>%
             group_by(Kingdom, Compartment, Fate) %>%
             mutate(Percent = round(Count / sum(Count) * 100, 1)) %>%
             arrange(Kingdom, Compartment, Fate, desc(Percent)) %>%
             rename(Abundance = Count),
           "Fig2_alluvial_taxonomic_summary")


# ==============================================================================
# 5. FIGURE 3 - ORCHARD-TO-FERMENTATION TRACEABILITY
# ==============================================================================

origins <- c("Olive", "Leaves", "Bark", "Rhizosphere", "Bulksoil")

# ------------------------------------------------------------------------------
# 5.1 Figure 3A - fraction of each BrinedOlive replicate's ASVs found in the
#     pooled origin compartment of the same tree (rarefied, 20 iterations)
# ------------------------------------------------------------------------------
calc_fermentation_traceability <- function(tab, meta, kingdom, source_year) {
  message(kingdom, " fermentation traceability (origins ", source_year, ")")
  trees <- sort(unique(as.character(meta$Tree)))

  is_origin <- meta$Compartment %in% origins & meta$Time == 0 & meta$Year == source_year
  is_ferm   <- meta$Compartment == "BrinedOlive" & meta$Time > 0
  idx       <- which(is_origin | is_ferm)
  rar       <- rarefy_presence(tab[rownames(meta)[idx], , drop = FALSE], meta[idx, ],
                               pa_min_reads, pa_rare_repeats)
  m <- rar$meta

  res <- list()
  for (it in seq_along(rar$mats)) {
    mat <- rar$mats[[it]]
    for (tree in trees) {
      ferm_samples <- intersect(rownames(m)[m$Tree == tree & m$Compartment == "BrinedOlive" & m$Time > 0],
                                rownames(mat))
      for (samp in ferm_samples) {
        asvs_ferm <- colnames(mat)[mat[samp, ] > 0]
        if (length(asvs_ferm) == 0) next
        for (origin in origins) {
          origin_samples <- rownames(m)[m$Tree == tree & m$Compartment == origin &
                                          m$Time == 0 & m$Year == source_year]
          pool <- get_present_asvs(mat, origin_samples)
          if (length(pool) == 0) next
          n_shared <- length(intersect(asvs_ferm, pool))
          res[[length(res) + 1]] <- data.frame(
            Kingdom = kingdom, Source_Year = source_year, Iteration = it, Tree = tree,
            Origin = origin, Time = as.numeric(as.character(m[samp, "Time"])), Sample = samp,
            N_Pool_Origin = length(pool), N_Ferm_Rep = length(asvs_ferm),
            N_Shared = n_shared, N_Not_In_Origin = length(asvs_ferm) - n_shared,
            Fraction = n_shared / length(asvs_ferm))
        }
      }
    }
  }

  bind_rows(res) %>%
    group_by(Kingdom, Source_Year, Tree, Origin, Time, Sample) %>%
    summarise(across(c(N_Pool_Origin, N_Ferm_Rep, N_Shared, N_Not_In_Origin, Fraction), mean),
              .groups = "drop") %>%
    mutate(Tree = factor(Tree), Origin = factor(Origin, levels = origins))
}

ferm_trace <- list()
for (yr in c("2022", "2023")) {
  ferm_trace[[yr]] <- bind_rows(lapply(names(kingdom_names), function(k)
    calc_fermentation_traceability(filt_tab[[k]], filt_meta[[k]], kingdom_names[[k]], yr)))
  save_table(ferm_trace[[yr]], paste0("Fig3A_traceability_", yr, "origins"))
  save_table(ferm_trace[[yr]] %>%
               group_by(Kingdom, Origin) %>%
               summarise(Mean_Pct = mean(Fraction) * 100, SD_Pct = sd(Fraction) * 100,
                         .groups = "drop") %>%
               arrange(Kingdom, desc(Mean_Pct)),
             paste0("Fig3A_traceability_", yr, "origins_summary"))

  p <- ggplot(ferm_trace[[yr]], aes(x = Origin, y = Fraction * 100, color = Tree)) +
    geom_jitter(size = 3, alpha = 0.85, width = 0.15) +
    stat_summary(fun = mean, geom = "crossbar", aes(group = 1),
                 width = 0.5, color = "black", linewidth = 0.6) +
    scale_color_manual(values = color_palette_tree, name = "Tree") +
    scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1)) +
    facet_wrap(~ Kingdom, ncol = 2) +
    labs(x = paste0("Origin Compartment (", yr, ")"),
         y = "% of fermentation replicate ASVs also present in origin") +
    theme_panel()
  save_fig(p, paste0("Fig3A_traceability_", yr, "origins"), width = 10, height = 5)
}

# ------------------------------------------------------------------------------
# 5.2 Figure 3B - traceability over fermentation time per tree
#     (2023 origins, non-rarefied presence/absence; replicates pooled per
#      tree x compartment x timepoint)
# ------------------------------------------------------------------------------
build_presence_df <- function(tab, meta, kingdom, year = "2023") {
  m <- meta[meta$Year == year, ]
  groups <- m %>% distinct(Tree, Compartment, Time)
  bind_rows(lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    samples <- rownames(m)[m$Tree == g$Tree & m$Compartment == g$Compartment & m$Time == g$Time]
    asvs <- get_present_asvs(tab, samples)
    if (length(asvs) == 0) return(NULL)
    data.frame(Kingdom = kingdom, Tree = as.character(g$Tree),
               Compartment = as.character(g$Compartment), Time = as.numeric(g$Time), ASV = asvs)
  }))
}

calc_traceability_by_time <- function(presence) {
  ferm_times <- sort(unique(presence$Time[presence$Compartment == "BrinedOlive" & presence$Time > 0]))
  res <- list()
  for (tree in unique(presence$Tree)) {
    origin_asvs <- lapply(setNames(origins, origins), function(comp)
      unique(presence$ASV[presence$Tree == tree & presence$Compartment == comp & presence$Time == 0]))
    for (tp in ferm_times) {
      ferm_asvs <- unique(presence$ASV[presence$Tree == tree &
                                         presence$Compartment == "BrinedOlive" & presence$Time == tp])
      if (length(ferm_asvs) == 0) next
      for (comp in origins) {
        shared <- length(intersect(ferm_asvs, origin_asvs[[comp]]))
        res[[length(res) + 1]] <- data.frame(
          Kingdom = presence$Kingdom[1], Tree = tree, Time = tp, Origin = comp,
          N_Origin = length(origin_asvs[[comp]]), N_Fermentation = length(ferm_asvs),
          N_Shared = shared, Fraction = shared / length(ferm_asvs))
      }
    }
  }
  bind_rows(res)
}

trace_time <- bind_rows(lapply(names(kingdom_names), function(k)
  calc_traceability_by_time(build_presence_df(filt_tab[[k]], filt_meta[[k]], kingdom_names[[k]])))) %>%
  mutate(Origin = factor(Origin, levels = origins))
save_table(trace_time, "Fig3B_traceability_by_tree")

color_palette_origin <- color_palette_compartment[origins]

p_trace_time <- trace_time %>%
  filter(!(Origin == "Olive" & Tree == "3")) %>%   # Tree 3 has no 2023 olive samples after QC
  ggplot(aes(x = Time, y = Fraction * 100, color = Origin, group = Origin)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  facet_grid(Tree ~ Kingdom, labeller = labeller(Tree = function(x) paste("Tree", x))) +
  scale_color_manual(values = color_palette_origin, name = "Origin\nCompartment") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Fermentation Time (Days)", y = "Fermentation ASVs traceable to origin (%)") +
  theme_panel()
save_fig(p_trace_time, "Fig3B_traceability_by_tree", width = 10, height = 10)

# ------------------------------------------------------------------------------
# 5.3 Figure 3C (and 4D) - taxonomic composition of ASV intersections
#     Non-rarefied, pooled across trees. Bacteria at family, fungi at genus
#     level; unclassified ASVs excluded; taxa < 2 % of a bar -> "Other".
# ------------------------------------------------------------------------------
tax_level_bar <- c(Bacteria = "family", Fungi = "genus")

make_bar_df <- function(asvs, weights, biom, kingdom, group_label) {
  if (length(asvs) == 0) return(NULL)
  labels <- asv_tax_labels(asvs, biom, tax_level_bar[[kingdom]])
  valid  <- !is.na(labels)
  if (sum(valid) == 0) return(NULL)
  data.frame(Taxon = labels[valid], Count = weights[asvs[valid]]) %>%
    group_by(Taxon) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    mutate(Pct = Count / sum(Count) * 100,
           Taxon = if_else(Pct < other_threshold, "Other", Taxon)) %>%
    group_by(Taxon) %>%
    summarise(Count = sum(Count), Pct = sum(Pct), .groups = "drop") %>%
    mutate(Kingdom = kingdom, Group = group_label, .before = 1)
}

# Cross-season: 2022 ∩ 2023 per compartment, weighted by mean 2022/2023 reads
calc_cs_intersection_bar <- function(tab, meta, biom, kingdom) {
  bind_rows(lapply(tree_compartments, function(comp) {
    s22 <- intersect(rownames(meta)[meta$Compartment == comp & meta$Year == "2022"], rownames(tab))
    s23 <- intersect(rownames(meta)[meta$Compartment == comp & meta$Year == "2023"], rownames(tab))
    shared <- intersect(get_present_asvs(tab, s22), get_present_asvs(tab, s23))
    if (length(shared) == 0) return(NULL)
    w <- (colSums(tab[s22, shared, drop = FALSE]) + colSums(tab[s23, shared, drop = FALSE])) / 2
    make_bar_df(shared, w, biom, kingdom, comp)
  }))
}

# Fermentation: origin (T0 of `source_year`) ∩ pooled BrinedOlive (T > 0),
# weighted by fermentation reads
calc_ferm_intersection_bar <- function(tab, meta, biom, kingdom, source_year) {
  ferm_samples <- intersect(rownames(meta)[meta$Compartment == "BrinedOlive" & meta$Time > 0], rownames(tab))
  counts_ferm  <- colSums(tab[ferm_samples, , drop = FALSE])
  asvs_ferm    <- names(counts_ferm)[counts_ferm > 0]
  bind_rows(lapply(tree_compartments, function(origin) {
    s <- intersect(rownames(meta)[meta$Compartment == origin & meta$Year == source_year & meta$Time == 0],
                   rownames(tab))
    shared <- intersect(asvs_ferm, get_present_asvs(tab, s))
    make_bar_df(shared, counts_ferm[shared], biom, kingdom, origin)
  })) %>% mutate(Source_Year = source_year)
}

cs_bar <- bind_rows(lapply(names(kingdom_names), function(k)
  calc_cs_intersection_bar(filt_tab[[k]], filt_meta[[k]], filt_biom[[k]], kingdom_names[[k]])))
ferm_bar <- bind_rows(lapply(names(kingdom_names), function(k) bind_rows(
  calc_ferm_intersection_bar(filt_tab[[k]], filt_meta[[k]], filt_biom[[k]], kingdom_names[[k]], "2022"),
  calc_ferm_intersection_bar(filt_tab[[k]], filt_meta[[k]], filt_biom[[k]], kingdom_names[[k]], "2023"))))

save_table(cs_bar,   "Fig4D_intersection_taxonomy_cross_season")
save_table(ferm_bar, "Fig3C_intersection_taxonomy_fermentation")

make_bar_panel <- function(df, palette, title, y_label = TRUE, tax_label = "Family") {
  df <- df %>%
    group_by(Group, Taxon) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(Pct = Count / sum(Count) * 100) %>%
    ungroup() %>%
    mutate(Group = factor(Group, levels = tree_compartments),
           Taxon = factor(Taxon, levels = rev(c(sort(setdiff(names(palette), "Other")), "Other"))))

  legend_face <- setNames(rep("italic", length(palette)), names(palette))
  legend_face["Other"] <- "plain"

  p <- ggplot(df, aes(x = Group, y = Pct, fill = Taxon)) +
    geom_col(width = 0.72, color = "gray30", linewidth = 0.2) +
    scale_fill_manual(values = palette, name = tax_label,
                      guide = guide_legend(reverse = TRUE, ncol = 2)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = NULL, y = if (y_label) "% of intersection reads (classified)" else NULL) +
    theme_bw(base_size = 16) +
    theme(plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x        = element_text(angle = 35, hjust = 1, size = 14),
          axis.text.y        = element_text(size = 13),
          axis.title.y       = element_text(size = 15, margin = margin(r = 12)),
          legend.text        = element_text(size = 9.5, face = legend_face),
          legend.title       = element_text(size = 12, face = "bold"),
          legend.key.size    = unit(0.4, "cm"),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank())
  if (!y_label) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

# Figure 3C: one palette per kingdom covering both cross-season and fermentation taxa
ferm_bar_row <- function(kingdom) {
  tax_label <- if (kingdom == "Bacteria") "Family" else "Genus"
  pal <- build_palette(c(cs_bar$Taxon[cs_bar$Kingdom == kingdom],
                         ferm_bar$Taxon[ferm_bar$Kingdom == kingdom]), base_cols_tab20)
  d   <- ferm_bar %>% filter(Kingdom == kingdom)
  p22 <- make_bar_panel(filter(d, Source_Year == "2022"), pal,
                        paste(kingdom, "-\n2022 origins ∩ BrinedOlive"), TRUE, tax_label)
  p23 <- make_bar_panel(filter(d, Source_Year == "2023"), pal,
                        paste(kingdom, "-\n2023 origins ∩ BrinedOlive"), FALSE, tax_label)
  # Legend from both years combined so that every plotted taxon is listed
  legend <- get_legend(make_bar_panel(d, pal, "", TRUE, tax_label))
  panels <- plot_grid(p22 + theme(legend.position = "none"), p23 + theme(legend.position = "none"),
                      ncol = 2, rel_widths = c(1.05, 1), align = "h", axis = "tb")
  plot_grid(panels, legend, ncol = 2, rel_widths = c(2.2, 1.1))
}

p_ferm_bar <- plot_grid(ferm_bar_row("Bacteria"), ferm_bar_row("Fungi"),
                        ncol = 1, labels = c("A", "B"), label_size = 22, label_x = 0.01)
save_fig(p_ferm_bar, "Fig3C_intersection_taxonomy_fermentation", width = 16.5, height = 13)


# ==============================================================================
# 6. FIGURE 4 - CROSS-SEASON DYNAMICS (AREA A, 2022 vs 2023)
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 Figure 4A + Supplementary Table 2A - alpha diversity
#     Simpson's index and richness averaged over 100 rarefactions to 500 reads.
#     (div_alpha is reused for Figure 5A.)
# ------------------------------------------------------------------------------
div_alpha <- list()
for (k in names(kingdom_names)) {
  rar <- rtk(filt_tab[[k]], repeats = alpha_rare_repeats, depth = alpha_rare_depth,
             margin = 1, threads = threads_num, verbose = FALSE, seed = seed)
  m <- filt_meta[[k]]
  stopifnot(length(rar$divvs) == nrow(filt_tab[[k]]))
  div_alpha[[k]] <- data.frame(
    Sample      = rownames(filt_tab[[k]]),
    Area        = m$Area,
    Tree        = m$Tree,
    Compartment = m$Compartment,
    Time        = factor(m$Time),
    Year        = m$Year,
    Richness    = sapply(rar$divvs, function(x) mean(x$richness)),
    Diversity   = sapply(rar$divvs, function(x) mean(x$simpson))
  )
}

# Model selection: fixed (Tree + Compartment * Year) vs Tree-nested structure;
# the nested model is kept if it fits significantly better (p < 0.05).
# Post-hoc: 2022 vs 2023 within each compartment (FDR, compact letter display).
alpha_season_anova <- list()
alpha_season_cld   <- list()
for (k in names(kingdom_names)) {
  d <- div_alpha[[k]] %>%
    filter(Compartment %in% tree_compartments) %>%
    mutate(Compartment = factor(Compartment, levels = tree_compartments))

  lm_fix    <- lm(Diversity ~ Tree + Compartment * Year, data = d)
  lm_nested <- lm(Diversity ~ Tree + Tree/Compartment + Compartment * Year, data = d)
  p_nested  <- anova(lm_fix, lm_nested)$`Pr(>F)`[2]
  use_nested <- !is.na(p_nested) && p_nested < 0.05
  fit <- if (use_nested) lm_nested else lm_fix

  alpha_season_anova[[k]] <- as.data.frame(anova(fit)) %>%
    rownames_to_column("Term") %>%
    mutate(Kingdom = kingdom_names[[k]],
           Model = if (use_nested) "Tree + Tree/Compartment + Compartment*Year"
                   else "Tree + Compartment*Year",
           .before = 1)

  emm <- emmeans(fit, ~ Year | Compartment, data = d)
  alpha_season_cld[[k]] <- as.data.frame(
    multcomp::cld(emm, alpha = 0.05, Letters = letters, sort = TRUE, adjust = "fdr")) %>%
    mutate(Kingdom = kingdom_names[[k]], Metric = "Simpson's Diversity")
}
save_table(bind_rows(alpha_season_anova), "SuppTable2A_alpha_season_anova")
save_table(bind_rows(alpha_season_cld),   "SuppTable2A_alpha_season_cld")

alpha_plots <- lapply(names(kingdom_names), function(k) {
  div_alpha[[k]] %>%
    filter(Compartment %in% tree_compartments) %>%
    mutate(Compartment = factor(Compartment, levels = tree_compartments)) %>%
    ggplot(aes(x = Compartment, y = Diversity, fill = Year, color = Year)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.65,
                 position = position_dodge(width = 0.75)) +
    geom_jitter(aes(group = Year), alpha = 0.7, size = 1.5,
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) +
    scale_fill_manual(values = color_palette_season, name = "Season") +
    scale_color_manual(values = color_palette_season, name = "Season") +
    labs(title = kingdom_names[[k]], x = "Compartment", y = "Simpson's Index") +
    theme_bw(base_size = 18) +
    theme(axis.text.x      = element_text(angle = 35, hjust = 1, size = 18),
          axis.text.y      = element_text(size = 18),
          axis.title       = element_text(size = 20),
          legend.text      = element_text(size = 16),
          legend.title     = element_text(size = 18),
          plot.title       = element_text(size = 18),
          panel.grid.minor = element_blank())
})
save_fig(plot_grid(plotlist = alpha_plots, ncol = 2),
         "Fig4A_alpha_diversity_season", width = 14, height = 5)

# ------------------------------------------------------------------------------
# 6.2 Figure 4B + Supplementary Table 2B - beta diversity
#     Samples < beta_min_reads dropped; 100 rarefactions to the smallest
#     library, averaged; Bray-Curtis; PCoA of both seasons together.
# ------------------------------------------------------------------------------
rarefy_average <- function(tab, repeats) {
  depth <- min(rowSums(tab))
  message("  Rarefying to ", depth, " reads (", repeats, " iterations, averaged)")
  rar <- rtk(tab, repeats = repeats, depth = depth, ReturnMatrix = repeats,
             margin = 1, threads = threads_num, verbose = FALSE, seed = seed)
  avg <- as.data.frame((Reduce(`+`, rar$raremat) / repeats)[rownames(tab), ])
  avg[, colSums(avg) > 0]
}

beta_tab  <- list()
beta_meta <- list()
for (k in names(kingdom_names)) {
  tab <- filt_tab[[k]]
  tab <- tab[rowSums(tab) >= beta_min_reads, ]
  beta_tab[[k]]  <- rarefy_average(tab, beta_rare_repeats)
  beta_meta[[k]] <- filt_meta[[k]][rownames(beta_tab[[k]]), ]
}

# Restrict to tree compartments
tree_subset <- function(k) {
  idx <- beta_meta[[k]]$Compartment %in% tree_compartments
  tab <- beta_tab[[k]][idx, ]
  list(tab  = tab[, colSums(tab) > 0],
       meta = beta_meta[[k]][idx, ] %>%
         mutate(Compartment = factor(Compartment, levels = tree_compartments)))
}

pcoa_plots <- lapply(names(kingdom_names), function(k) {
  s   <- tree_subset(k)
  ord <- capscale(s$tab ~ 1, distance = "bray")
  var_exp <- round(ord$CA$eig / sum(ord$CA$eig) * 100, 1)
  scrs <- as.data.frame(scores(ord, display = "sites", choices = 1:2)) %>% bind_cols(s$meta)

  ggplot(scrs, aes(x = MDS1, y = MDS2, fill = Compartment, shape = Year)) +
    geom_point(size = 4, stroke = 0.6, alpha = 0.85) +
    scale_fill_manual(values = color_palette_compartment) +
    scale_shape_manual(values = c("2022" = 21, "2023" = 24), name = "Season") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4))) +
    labs(title = kingdom_names[[k]],
         x = paste0("PCoA 1 (", var_exp[1], "%)"), y = paste0("PCoA 2 (", var_exp[2], "%)")) +
    theme_bw(base_size = 18) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
          legend.text = element_text(size = 14), legend.title = element_text(size = 16),
          plot.title = element_text(size = 18), panel.grid.minor = element_blank())
})
save_fig(plot_grid(plotlist = pcoa_plots, ncol = 2), "Fig4B_pcoa", width = 16, height = 6)

# PERMANOVA (sequential terms): dist ~ Year + Tree + Compartment, and per
# compartment dist ~ Year + Tree (Year effect FDR-adjusted across compartments)
permanova <- list()
permanova_year <- list()
for (k in names(kingdom_names)) {
  s <- tree_subset(k)
  set.seed(seed)
  perm <- adonis2(vegdist(s$tab, method = "bray") ~ Year + Tree + Compartment,
                  data = s$meta, permutations = n_perm, by = "terms")
  permanova[[k]] <- tidy_adonis2(perm, kingdom_names[[k]], "Tree compartments",
                                 "dist ~ Year + Tree + Compartment")

  rows <- list()
  for (comp in tree_compartments) {
    idx <- which(s$meta$Compartment == comp)
    if (length(idx) < 3 || length(unique(s$meta$Year[idx])) < 2) next
    tab_c <- s$tab[idx, ]
    tab_c <- tab_c[, colSums(tab_c) > 0]
    set.seed(seed)
    perm_c <- adonis2(vegdist(tab_c, method = "bray") ~ Year + Tree,
                      data = s$meta[idx, ], permutations = n_perm, by = "terms")
    yr <- as.data.frame(perm_c)["Year", ]
    rows[[comp]] <- data.frame(Kingdom = kingdom_names[[k]], Compartment = comp, Df = yr$Df,
                               SumOfSqs = yr$SumOfSqs, R2 = yr$R2, F = yr$F, p_value = yr$`Pr(>F)`)
  }
  permanova_year[[k]] <- bind_rows(rows) %>% mutate(p_adj = p.adjust(p_value, method = "fdr"))
}
save_table(bind_rows(permanova),      "SuppTable2B_permanova")
save_table(bind_rows(permanova_year), "SuppTable2B_permanova_year_within_compartment")

# ------------------------------------------------------------------------------
# 6.3 Figure 4C + Supplementary Table 3 - cross-season ASV persistence
#     Per 2023 replicate: fraction of its ASVs also detected in the pooled 2022
#     samples of the same tree x compartment (rarefied, 20 iterations).
# ------------------------------------------------------------------------------
cs_compartments <- c("Bark", "Leaves", "Rhizosphere", "Bulksoil", "Olive")

calc_cross_season_persistence <- function(tab, meta, kingdom) {
  message(kingdom, " cross-season persistence")
  trees <- sort(unique(as.character(meta$Tree)))
  idx <- which(meta$Compartment %in% cs_compartments & meta$Year %in% c("2022", "2023"))
  rar <- rarefy_presence(tab[rownames(meta)[idx], , drop = FALSE], meta[idx, ],
                         pa_min_reads, pa_rare_repeats)
  m <- rar$meta

  res <- list()
  for (it in seq_along(rar$mats)) {
    mat <- rar$mats[[it]]
    for (tree in trees) {
      for (comp in cs_compartments) {
        in_group <- m$Tree == tree & m$Compartment == comp
        pool_2022 <- get_present_asvs(mat, rownames(m)[in_group & m$Year == "2022"])
        if (length(pool_2022) == 0) next
        for (samp in intersect(rownames(m)[in_group & m$Year == "2023"], rownames(mat))) {
          asvs_2023 <- colnames(mat)[mat[samp, ] > 0]
          if (length(asvs_2023) == 0) next
          n_shared <- length(intersect(asvs_2023, pool_2022))
          res[[length(res) + 1]] <- data.frame(
            Kingdom = kingdom, Iteration = it, Tree = tree, Compartment = comp, Sample = samp,
            N_Pool_2022 = length(pool_2022), N_2023_Rep = length(asvs_2023),
            N_Shared = n_shared, N_Not_In_2022 = length(asvs_2023) - n_shared,
            Fraction = n_shared / length(asvs_2023))
        }
      }
    }
  }

  bind_rows(res) %>%
    group_by(Kingdom, Tree, Compartment, Sample) %>%
    summarise(across(c(N_Pool_2022, N_2023_Rep, N_Shared, N_Not_In_2022, Fraction), mean),
              .groups = "drop") %>%
    mutate(Tree = factor(Tree), Compartment = factor(Compartment, levels = cs_compartments))
}

cs_trace <- bind_rows(lapply(names(kingdom_names), function(k)
  calc_cross_season_persistence(filt_tab[[k]], filt_meta[[k]], kingdom_names[[k]])))
save_table(cs_trace, "Fig4C_cross_season_persistence")
save_table(cs_trace %>%
             group_by(Kingdom, Compartment) %>%
             summarise(N_Replicates = n(),
                       Mean_Rep_Richness = mean(N_2023_Rep),
                       Mean_Pct_From_2022 = mean(Fraction) * 100,
                       SD_Pct_From_2022 = sd(Fraction) * 100,
                       .groups = "drop"),
           "Fig4C_cross_season_persistence_summary")

p_cs <- ggplot(cs_trace, aes(x = Compartment, y = Fraction * 100, color = Tree)) +
  geom_jitter(size = 3, alpha = 0.85, width = 0.15) +
  stat_summary(fun = mean, geom = "crossbar", aes(group = 1),
               width = 0.5, color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_palette_tree, name = "Tree") +
  scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ Kingdom, ncol = 2) +
  labs(x = "Compartment", y = "% of 2023 ASVs also present in 2022") +
  theme_panel()
save_fig(p_cs, "Fig4C_cross_season_persistence", width = 10, height = 5)

# Binomial GLMM: cbind(N_Shared, N_Lost) ~ Compartment + (1|Tree) + (1|obs_id)
# (obs_id = observation-level random effect for overdispersion). Compartment
# tested by LRT against the null model; post-hoc contrasts with FDR.
glmm_lrt <- list(); glmm_cld <- list(); glmm_posthoc <- list()
for (king in kingdom_names) {
  d <- cs_trace %>%
    filter(Kingdom == king) %>%
    mutate(N_Shared_int = round(N_Shared), N_Lost_int = round(N_Not_In_2022)) %>%
    filter(N_Shared_int + N_Lost_int > 0) %>%
    mutate(obs_id = row_number())

  fit  <- glmer(cbind(N_Shared_int, N_Lost_int) ~ Compartment + (1 | Tree) + (1 | obs_id),
                data = d, family = binomial)
  null <- glmer(cbind(N_Shared_int, N_Lost_int) ~ 1 + (1 | Tree) + (1 | obs_id),
                data = d, family = binomial)

  glmm_lrt[[king]] <- as.data.frame(anova(null, fit)) %>%
    rownames_to_column("Model") %>% mutate(Kingdom = king, .before = 1)

  emm <- emmeans(fit, ~ Compartment, type = "response")
  glmm_cld[[king]] <- as.data.frame(
    multcomp::cld(emm, alpha = 0.05, Letters = letters, sort = TRUE, adjust = "fdr")) %>%
    mutate(Kingdom = king)
  glmm_posthoc[[king]] <- bind_rows(
    as.data.frame(emm) %>% mutate(Type = "Marginal estimate (response scale)"),
    as.data.frame(summary(pairs(emm, adjust = "fdr"))) %>% mutate(Type = "Pairwise contrast (FDR-adjusted)")
  ) %>% mutate(Kingdom = king)
}
save_table(bind_rows(glmm_lrt),     "SuppTable3_persistence_glmm_lrt")
save_table(bind_rows(glmm_posthoc), "SuppTable3_persistence_glmm_posthoc")
save_table(bind_rows(glmm_cld),     "SuppTable3_persistence_glmm_cld")

# ------------------------------------------------------------------------------
# 6.4 Figure 4D - taxonomy of ASVs shared between 2022 and 2023
#     (data computed in Section 5.3; palette built from cross-season taxa only)
# ------------------------------------------------------------------------------
cs_bar_row <- function(kingdom) {
  d   <- cs_bar %>% filter(Kingdom == kingdom)
  pal <- build_palette(d$Taxon, base_cols_tab20)
  p   <- make_bar_panel(d, pal, paste(kingdom, "- Cross-season (2022 ∩ 2023)"), TRUE,
                        if (kingdom == "Bacteria") "Family" else "Genus")
  plot_grid(p + theme(legend.position = "none"), get_legend(p), ncol = 2, rel_widths = c(1, 1.2))
}
p_cs_bar <- plot_grid(cs_bar_row("Bacteria"), cs_bar_row("Fungi"),
                      ncol = 1, labels = c("A", "B"), label_size = 22, label_x = 0.01)
save_fig(p_cs_bar, "Fig4D_intersection_taxonomy_cross_season", width = 13.5, height = 13)


# ==============================================================================
# 7. FIGURE 5 - COMMUNITY DYNAMICS DURING FERMENTATION (2023)
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Figure 5A + Supplementary Table 4 - richness and Simpson over time
#     Olive T0 + BrinedOlive T3-T145, rarefied values from Section 6.1.
#     Model: metric ~ kingdom * Time + Tree:kingdom; bacteria vs fungi
#     contrasts at each timepoint (FDR).
# ------------------------------------------------------------------------------
ts_alpha <- bind_rows(div_alpha, .id = "kingdom") %>%
  filter(Compartment %in% c("Olive", "BrinedOlive"), Year == "2023") %>%
  mutate(Time          = as.numeric(as.character(Time)),
         kingdom_label = unname(kingdom_names[kingdom]),
         Time_factor   = factor(Time),
         Tree_kingdom  = interaction(Tree, kingdom_label, drop = TRUE))

kingdom_anova <- list(); kingdom_contrasts <- list()
for (metric in c("Richness", "Diversity")) {
  p <- ggplot(ts_alpha, aes(x = Time, y = .data[[metric]])) +
    geom_point(aes(color = kingdom_label, shape = Tree), size = 3) +
    geom_smooth(aes(group = kingdom_label, color = kingdom_label, fill = kingdom_label),
                method = "loess", alpha = 0.2) +
    scale_color_manual(values = color_palette_kingdom, name = "Kingdom") +
    scale_fill_manual(values = color_palette_kingdom, name = "Kingdom") +
    scale_shape_manual(values = c(16, 17, 15), name = "Tree") +
    labs(title = paste("Combined", metric, "Dynamics during Fermentation"),
         x = "Time (Days)", y = metric) +
    theme_bw()
  save_fig(p, paste0("Fig5A_", tolower(metric), "_vs_time"), width = 8, height = 6)

  fit <- lm(as.formula(paste(metric, "~ kingdom_label * Time_factor + Tree_kingdom")), data = ts_alpha)
  kingdom_anova[[metric]] <- as.data.frame(anova(fit)) %>%
    rownames_to_column("Term") %>% mutate(Metric = metric, .before = 1)
  kingdom_contrasts[[metric]] <- as.data.frame(summary(
    pairs(emmeans(fit, ~ kingdom_label | Time_factor), adjust = "fdr"))) %>%
    mutate(Metric = metric)
}
save_table(bind_rows(kingdom_anova),     "SuppTable4_ferm_alpha_kingdom_anova")
save_table(bind_rows(kingdom_contrasts), "SuppTable4_ferm_alpha_kingdom_contrasts")

# ------------------------------------------------------------------------------
# 7.2 Figure 5B - ASV turnover between consecutive timepoints, per tree
#     Olive + BrinedOlive samples >= 500 reads, rarefied (20 iterations) to the
#     smallest library. Replicates are pooled per tree x timepoint; ASVs at T_i
#     vs T_(i-1) are classified maintained / gained / lost.
#     NOTE: no year filter - the T0 pool contains 2022 and 2023 olives.
# ------------------------------------------------------------------------------
for (k in names(kingdom_names)) {
  title_prefix <- if (k == "bac") "Bacterial" else "Fungal"

  m <- div_alpha[[k]] %>%
    filter(Compartment %in% c("Olive", "BrinedOlive")) %>%
    mutate(Time = as.numeric(as.character(Time)))
  rownames(m) <- m$Sample
  rar <- rarefy_presence(filt_tab[[k]][m$Sample, ], m, pa_min_reads, pa_rare_repeats)
  m <- rar$meta

  res <- list()
  for (it in seq_along(rar$mats)) {
    mat <- rar$mats[[it]]
    for (tr in sort(unique(m$Tree))) {
      mt <- m[m$Tree == tr, ]
      timepoints <- sort(unique(mt$Time))
      asvs_at <- function(t) get_present_asvs(mat, mt$Sample[mt$Time == t])
      for (i in seq_along(timepoints)) {
        now <- asvs_at(timepoints[i])
        if (i == 1) {
          res[[length(res) + 1]] <- data.frame(
            Iteration = it, Tree = tr, Time = timepoints[i],
            Category = if (timepoints[i] == 0) "Baseline (T0)" else paste0("Start (T", timepoints[i], ")"),
            Count = length(now))
        } else {
          prev <- asvs_at(timepoints[i - 1])
          res[[length(res) + 1]] <- data.frame(
            Iteration = it, Tree = tr, Time = timepoints[i],
            Category = c("Maintained", "Gained", "Lost"),
            Count = c(length(intersect(now, prev)), length(setdiff(now, prev)),
                      -length(setdiff(prev, now))))
        }
      }
    }
  }

  turnover <- bind_rows(res) %>%
    group_by(Tree, Time, Category) %>%
    summarise(Mean_Count = mean(Count), SD_Count = sd(Count), .groups = "drop")
  save_table(turnover, paste0("Fig5B_turnover_", tolower(title_prefix)))

  p <- ggplot(filter(turnover, Time > 0), aes(x = factor(Time), y = Mean_Count, fill = Category)) +
    geom_col() +
    geom_hline(yintercept = 0, color = "black") +
    scale_fill_manual(values = color_palette_change[c("Gained", "Lost", "Maintained")]) +
    facet_wrap(~ Tree) +
    labs(title = paste0(title_prefix, " ASV Turnover by Tree"),
         x = "Time (Days)", y = "Number of ASVs") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold"),
          legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  save_fig(p, paste0("Fig5B_turnover_", tolower(title_prefix)), width = 10, height = 6)
}

# ------------------------------------------------------------------------------
# 7.3 Figure 5C - dominant-taxon trajectories between start and end
#     Non-rarefied relative abundances, mean over replicates at the first and
#     last timepoint of each tree. ASVs > 0.01 % at both ends = persistent, only
#     at start = lost, only at end = gained. Aggregated to family (bacteria) /
#     genus (fungi) with fallback; taxa > 5 % at either end are shown.
#     NOTE: no year filter - the start pool contains 2022 and 2023 olives.
# ------------------------------------------------------------------------------
rel_abund_threshold <- 0.01   # %
traj_min_pct        <- 5      # %

calc_trajectories <- function(tab, meta, biom, kingdom) {
  samples <- intersect(rownames(meta)[meta$Compartment %in% c("Olive", "BrinedOlive")], rownames(tab))
  rel <- prop.table(as.matrix(tab[samples, ]), margin = 1) * 100
  m   <- meta[samples, ]
  m$TimeNum <- as.numeric(as.character(m$Time))

  per_tree <- lapply(unique(m$Tree), function(tree) {
    mt <- m[m$Tree == tree, ]
    if (length(unique(mt$TimeNum)) < 2) return(NULL)
    start <- colMeans(rel[rownames(mt)[mt$TimeNum == min(mt$TimeNum)], , drop = FALSE])
    end   <- colMeans(rel[rownames(mt)[mt$TimeNum == max(mt$TimeNum)], , drop = FALSE])
    data.frame(Tree = tree, ASV = colnames(rel), Start_Abund = start, End_Abund = end) %>%
      mutate(Category = case_when(
        Start_Abund >  rel_abund_threshold & End_Abund >  rel_abund_threshold ~ "Persistent",
        Start_Abund >  rel_abund_threshold & End_Abund <= rel_abund_threshold ~ "Lost",
        Start_Abund <= rel_abund_threshold & End_Abund >  rel_abund_threshold ~ "Gained")) %>%
      filter(!is.na(Category))
  })

  bind_rows(per_tree) %>%
    mutate(Tree = factor(Tree),
           Taxon = get_fallback_label(biom[ASV, "taxonomy"], kingdom),
           Dataset = kingdom) %>%
    group_by(Tree, Taxon, Category, Dataset) %>%
    summarise(Start_Abund = sum(Start_Abund), End_Abund = sum(End_Abund), .groups = "drop")
}

traj <- bind_rows(lapply(names(kingdom_names), function(k)
  calc_trajectories(filt_tab[[k]], filt_meta[[k]], biom_sel[[k]], kingdom_names[[k]])))

# One line per taxon; colour = category of the ASVs contributing most at the end
traj_plot_data <- traj %>%
  group_by(Dataset, Tree, Taxon) %>%
  summarise(Total_Start = sum(Start_Abund), Total_End = sum(End_Abund),
            Dom_Cat_End = Category[which.max(End_Abund)], .groups = "drop") %>%
  mutate(Final_Category = ifelse(Total_End < 1, "Lost", Dom_Cat_End)) %>%
  filter(Total_Start > traj_min_pct | Total_End > traj_min_pct)
save_table(traj,           "Fig5C_trajectories_all_taxa")
save_table(traj_plot_data, "Fig5C_trajectories_plotted")

for (king in kingdom_names) {
  p <- ggplot(filter(traj_plot_data, Dataset == king)) +
    geom_segment(aes(x = 0, xend = 1, y = Total_Start, yend = Total_End, color = Final_Category),
                 linewidth = 1, alpha = 0.7) +
    geom_point(aes(x = 0, y = Total_Start), size = 3) +
    geom_point(aes(x = 1, y = Total_End), size = 3) +
    geom_text_repel(aes(x = 0, y = Total_Start, label = Taxon), hjust = 1.2, direction = "y", size = 4) +
    geom_text_repel(aes(x = 1, y = Total_End, label = paste0(round(Total_End, 1), "%")),
                    hjust = -0.2, direction = "y", size = 4) +
    facet_wrap(~ Tree) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End"), limits = c(-0.5, 1.5)) +
    scale_color_manual(values = color_palette_change[c("Gained", "Persistent", "Lost")]) +
    labs(title = paste0(king, " Trajectories (Family/Genus Aggregated)"),
         x = "", y = "Relative Abundance (%)") +
    theme_bw() +
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  save_fig(p, paste0("Fig5C_trajectories_", tolower(king)), width = 10, height = 6)
}


# ==============================================================================
# 8. SUPPLEMENTARY TABLE 1 - AREA COMPARISON (2022, AREAS A/B/C)
# ==============================================================================
# Dataset: 2022 tree-compartment samples of all areas, target taxa only,
# empty ASVs and samples with <= 100 reads removed (no decontamination step).

meta_area <- meta_full %>%
  filter(Year == "2022", Compartment %in% c("BulkSoil", "Bulksoil", tree_compartments[-1]),
         Area %in% c("A", "B", "C")) %>%
  mutate(Compartment = factor(ifelse(Compartment == "BulkSoil", "Bulksoil", Compartment),
                              levels = tree_compartments),
         Area = factor(Area, levels = c("A", "B", "C")),
         Tree = factor(Tree))

area_tab  <- list()
area_meta <- list()
for (k in names(kingdom_names)) {
  biom <- asv_raw[[k]][, c(intersect(rownames(meta_area), colnames(asv_raw[[k]])), "taxonomy")]
  biom <- if (k == "bac") {
    biom[grepl("k__Archaea|k__Bacteria", biom$taxonomy, ignore.case = TRUE) &
           !grepl("chloroplast|mitochondria", biom$taxonomy, ignore.case = TRUE), ]
  } else {
    biom[grepl("k__Fungi", biom$taxonomy, ignore.case = TRUE), ]
  }
  tab <- samples_x_asvs(biom)
  tab <- tab[, colSums(tab) > 0]
  area_tab[[k]]  <- tab[rowSums(tab) > min_reads_sample, ]
  area_meta[[k]] <- meta_area[rownames(area_tab[[k]]), ]
}

# --- 8.1 Alpha diversity: Diversity ~ Area / Tree per compartment ------------
# Area tested against the Tree(Area) mean square (quasi-F); pairwise area
# contrasts from emmeans (FDR).
area_anova <- list(); area_posthoc <- list(); area_cld <- list()
for (k in names(kingdom_names)) {
  rar <- rtk(area_tab[[k]], repeats = area_rare_repeats, depth = alpha_rare_depth,
             margin = 1, threads = threads_num, verbose = FALSE, seed = seed)
  stopifnot(length(rar$divvs) == nrow(area_tab[[k]]))
  alpha <- area_meta[[k]] %>%
    transmute(Area, Tree, Compartment,
              Diversity = sapply(rar$divvs, function(x) mean(x$simpson)))

  for (comp in tree_compartments) {
    d <- filter(alpha, Compartment == comp)
    key <- paste(k, comp)
    if (n_distinct(d$Area) < 2) next
    if (nrow(d) / n_distinct(paste(d$Area, d$Tree)) <= 1) {
      area_anova[[key]] <- data.frame(Kingdom = kingdom_names[[k]], Compartment = comp,
                                      Note = "No residual df (1 rep/tree)")
      next
    }
    fit <- lm(Diversity ~ Area / Tree, data = d)
    aov_tab <- anova(fit)
    F_area  <- aov_tab["Area", "Mean Sq"] / aov_tab["Area:Tree", "Mean Sq"]
    area_anova[[key]] <- data.frame(
      Kingdom = kingdom_names[[k]], Compartment = comp, Term = "Area",
      Df_Area = aov_tab["Area", "Df"], Df_TreeArea = aov_tab["Area:Tree", "Df"],
      MS_Area = aov_tab["Area", "Mean Sq"], MS_TreeArea = aov_tab["Area:Tree", "Mean Sq"],
      F_value = F_area,
      p_value = pf(F_area, aov_tab["Area", "Df"], aov_tab["Area:Tree", "Df"], lower.tail = FALSE))

    emm <- emmeans(fit, ~ Area, data = d)
    area_posthoc[[key]] <- as.data.frame(summary(pairs(emm, adjust = "fdr"))) %>%
      mutate(Kingdom = kingdom_names[[k]], Compartment = comp)
    area_cld[[key]] <- as.data.frame(
      multcomp::cld(emm, alpha = 0.05, Letters = letters, sort = TRUE, adjust = "fdr")) %>%
      mutate(Kingdom = kingdom_names[[k]], Compartment = comp)
  }
}
save_table(bind_rows(area_anova),   "SuppTable1A_alpha_area_anova")
save_table(bind_rows(area_posthoc), "SuppTable1A_alpha_area_posthoc")
save_table(bind_rows(area_cld),     "SuppTable1A_alpha_area_cld")

# --- 8.2 Beta diversity: PERMANOVA -------------------------------------------
# 20 rarefactions to the smallest library, averaged and rounded; Bray-Curtis.
# Global: dist ~ Area + Compartment + Area:Compartment
# Per compartment: dist ~ Area / Tree
area_permanova <- list()
for (k in names(kingdom_names)) {
  tab <- round(rarefy_average(area_tab[[k]], area_rare_repeats))
  tab <- tab[rowSums(tab) > 0, ]
  tab <- tab[, colSums(tab) > 0]
  m   <- area_meta[[k]][rownames(tab), ]
  king <- if (k == "bac") "Bacterial" else "Fungal"

  set.seed(seed)
  perm <- adonis2(vegdist(tab, method = "bray") ~ Area + Compartment + Area:Compartment,
                  data = m, permutations = n_perm, by = "terms")
  area_permanova[[paste(k, "global")]] <- tidy_adonis2(
    perm, king, "All compartments (2022)", "dist ~ Area + Compartment + Area:Compartment")

  for (comp in tree_compartments) {
    idx <- which(m$Compartment == comp)
    if (length(idx) < 3 || length(unique(m$Area[idx])) < 2) next
    tab_c <- tab[idx, ]
    tab_c <- tab_c[, colSums(tab_c) > 0]
    set.seed(seed)
    perm_c <- adonis2(vegdist(tab_c, method = "bray") ~ Area / Tree,
                      data = m[idx, ], permutations = n_perm, by = "terms")
    area_permanova[[paste(k, comp)]] <- tidy_adonis2(perm_c, king, paste0(comp, " (2022)"),
                                                     "dist ~ Area / Tree")
  }
}
save_table(bind_rows(area_permanova), "SuppTable1B_beta_area_permanova")


# ==============================================================================
# 9. SESSION INFO
# ==============================================================================
writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))
message("Done.")
