#===================================================#
#  Single clone analysis
#===================================================#
#
# This script analyses OTX2 expression in single clones of edited cells,
# measured by flow cytometry, and relates it to RNA-seq and gene-expression
# scores.
#
# Compared to the original script, repeated logic (sampling, normalising,
# summarising by replicate, plotting) has been extracted into small
# reusable functions. The two flow cytometry batches ("13/05/26" and
# "01/09/26") are processed with the same functions instead of duplicating
# the same block of code twice.

# ---------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------

packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel")
for (p in packages) {
  library(p, character.only = TRUE)
  cat(p, "version:", as.character(packageVersion(p)), "\n")
}

# ---------------------------------------------------------------------
# 1. Backup (run from a shell, not from R)
# ---------------------------------------------------------------------

# cd /lustre/scratch126/gengen/projects_v2/escramble
# git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/005_single_clone_analysis_v2_0.R
# git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/005_single_clone_analysis/flow_data
# git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/005_single_clone_analysis/rplots
# git commit -m "Last update committed on $(date +'%D %r')"
# git push

# ---------------------------------------------------------------------
# 2. Working directory
# ---------------------------------------------------------------------

# setwd("~/Desktop/mount_escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes")
# or
# setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes")

# ---------------------------------------------------------------------
# 3. Clone nomenclature (unchanged)
# ---------------------------------------------------------------------
#
# The clones are named after the deletion sites according to the
# following nomenclature:
#
#                                                                     1        2   3   4
#                                                                     |        |   |   |
# ------[--OTX2--]------------------------------------------------------------------------------
#                           |           |           |           |           |             |
#                           A           B           C           D           E             F
#                               R1           R2          R3          R4           R5

# ---------------------------------------------------------------------
# 4. Reusable functions
# ---------------------------------------------------------------------

#' Read a raw flow cytometry CSV and down-/up-sample every column to `n`
#' observations, so that distributions are comparable across clones.
load_and_resample <- function(path, n = 1500, seed = 42) {
  df <- read.csv(path, check.names = FALSE)
  set.seed(seed)
  as.data.frame(
    lapply(df, function(x) sample(x[!is.na(x)], size = n, replace = TRUE)),
    check.names = FALSE
  )
}

#' Reshape a wide clone-expression table to long format and normalise each
#' single-cell value to 0-100 % using the WT and HAP1 medians as anchors.
normalise_single_cells <- function(df, wt_pattern = "WT", hap1_pattern = "HAP1") {
  long <- df %>% pivot_longer(cols = everything(), names_to = "clone", values_to = "exp")
  
  wt_median   <- long %>% filter(grepl(wt_pattern,   clone)) %>% pull(exp) %>% median(na.rm = TRUE)
  hap1_median <- long %>% filter(grepl(hap1_pattern, clone)) %>% pull(exp) %>% median(na.rm = TRUE)
  
  long <- long %>% mutate(exp.norm = ((exp - hap1_median) / (wt_median - hap1_median)) * 100)
  
  # order clones by median normalised expression, for plotting
  clone_order <- long %>%
    group_by(clone) %>%
    summarise(median_exp.norm = median(exp.norm, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_exp.norm) %>%
    pull(clone)
  
  long %>% mutate(clone = factor(clone, levels = clone_order))
}

#' Split "clone_replicate" columns into clone + replicate, then summarise
#' per-replicate and per-clone normalised expression (mean, SD, SEM),
#' anchored on WT / HAP1 replicate means.
summarise_replicates <- function(long_df, wt_pattern = "WT", hap1_pattern = "HAP1") {
  
  rep_means <- long_df %>%
    group_by(clone) %>%
    summarise(rep.mean = mean(exp.norm, na.rm = TRUE), .groups = "drop") %>%
    separate(col = clone, into = c("clone", "rep"), sep = "_", fill = "right")
  
  wt_anchor   <- rep_means %>% filter(grepl(wt_pattern,   clone)) %>% pull(rep.mean) %>% mean(na.rm = TRUE)
  hap1_anchor <- rep_means %>% filter(grepl(hap1_pattern, clone)) %>% pull(rep.mean) %>% mean(na.rm = TRUE)
  
  rep_means <- rep_means %>%
    mutate(rep.mean.norm = ((rep.mean - hap1_anchor) / (wt_anchor - hap1_anchor)) * 100)
  
  clone_summary <- rep_means %>%
    group_by(clone) %>%
    summarise(
      exp.norm.mean = mean(rep.mean.norm, na.rm = TRUE),
      exp.norm.sd   = sd(rep.mean.norm,   na.rm = TRUE),
      exp.norm.sem  = sd(rep.mean.norm,   na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(clone = factor(clone, levels = arrange(., exp.norm.mean) %>% pull(clone)))
  
  # individual replicate points, ordered using the same clone levels
  rep_points <- rep_means %>%
    mutate(clone = factor(clone, levels = levels(clone_summary$clone)))
  
  list(clone_summary = clone_summary, rep_points = rep_points)
}

#' Boxplot of single-cell normalised expression per clone.
plot_single_cells <- function(long_df, highlight_pattern, title, subtitle) {
  cols <- setNames(
    ifelse(grepl(highlight_pattern, levels(long_df$clone)), "#ED1C24", "#E6E7E8"),
    levels(long_df$clone)
  )
  
  ggplot(long_df, aes(x = exp.norm, y = clone, fill = clone)) +
    geom_vline(xintercept = c(0, 100), linetype = "dashed", linewidth = 0.4, colour = "grey40") +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, outlier.colour = "grey50",
                 linewidth = 0.4, width = 0.7, colour = "grey20") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(limits = c(0, 300), breaks = seq(0, 300, 50), expand = c(0, 0)) +
    labs(x = "Relative OTX2 expression (% of WT)", y = NULL, title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 3.3,
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          axis.text = element_text(size = 9),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey30"))
}

#' Point + error-bar plot of replicate-mean normalised expression per clone.
plot_replicate_means <- function(clone_summary, rep_points, highlight_pattern, title, subtitle, x_max = 225) {
  cols <- setNames(
    ifelse(grepl(highlight_pattern, levels(clone_summary$clone)), "#ED1C24", "#E6E7E8"),
    levels(clone_summary$clone)
  )
  
  ggplot(clone_summary, aes(x = exp.norm.mean, y = clone)) +
    geom_vline(xintercept = c(0, 100), linetype = "dashed", linewidth = 0.4, colour = "grey40") +
    geom_errorbarh(aes(xmin = exp.norm.mean - exp.norm.sd, xmax = exp.norm.mean + exp.norm.sd),
                   height = 0.35, linewidth = 0.6, na.rm = TRUE, colour = "grey30") +
    geom_point(data = rep_points, aes(x = rep.mean.norm, fill = clone),
               shape = 21, size = 2.5, colour = "grey20", stroke = 0.4) +
    geom_point(aes(fill = clone), shape = 124, size = 7, colour = "black") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, 300, 50), expand = c(0, 0)) +
    labs(x = "Relative OTX2 expression (% of WT)", y = NULL, title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1.4,
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          axis.text = element_text(size = 9),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey30"))
}

#' Scatter plot with linear regression, error bars and text labels -
#' shared by the flow-vs-RNAseq correlation and the calibration plot.
plot_correlation <- function(df, x, y, xerr, yerr, label_col,
                             x_limits, y_limits, x_lab, y_lab, title, subtitle,
                             eb_width = 2, eb_height = 2) {
  x <- rlang::enquo(x); y <- rlang::enquo(y)
  xerr <- rlang::enquo(xerr); yerr <- rlang::enquo(yerr); label_col <- rlang::enquo(label_col)
  
  lm.fit <- lm(reformulate(rlang::as_label(x), rlang::as_label(y)), data = df)
  lm.lab <- sprintf("R\u00b2 = %.3f\np = %.3e", summary(lm.fit)$r.squared, summary(lm.fit)$coefficients[2, 4])
  
  ggplot(df, aes(x = !!x, y = !!y)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, colour = "grey30", fill = "grey80", alpha = 0.4) +
    geom_errorbar(aes(ymin = !!y - !!yerr, ymax = !!y + !!yerr), width = eb_width, linewidth = 0.4, colour = "grey40") +
    geom_errorbarh(aes(xmin = !!x - !!xerr, xmax = !!x + !!xerr), height = eb_height, linewidth = 0.4, colour = "grey40") +
    geom_point(shape = 21, size = 3, fill = "#ED1C24", colour = "grey20", stroke = 0.4) +
    geom_text_repel(aes(label = !!label_col), size = 3, colour = "grey20", box.padding = 0.4, max.overlaps = Inf) +
    annotate("text", x = -Inf, y = Inf, label = lm.lab, hjust = -0.1, vjust = 1.3, size = 3.2, colour = "grey20", fontface = "italic") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1,
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          axis.text = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey30"))
}

# ---------------------------------------------------------------------
# 5. Batch 1 - data from 13/05/26
# ---------------------------------------------------------------------

flow.df <- load_and_resample("./005_single_clone_analysis/flow_data/260513_OTX2_clones_single_cell_clean.csv")

flow.long.df.1 <- normalise_single_cells(flow.df)
rep.1 <- summarise_replicates(flow.long.df.1)

p1 <- plot_single_cells(flow.long.df.1, "WT",
                        "OTX2 expression across clones",
                        "Single-cell distributions; normalised to WT (100 %) and HAP1 (0 %)")

p2 <- plot_replicate_means(rep.1$clone_summary, rep.1$rep_points, "WT",
                           "OTX2 expression across clones",
                           "Replicate means \u00b1 SD; normalised to WT (100 %) and HAP1 (0 %)",
                           x_max = 225)

pdf("./005_single_clone_analysis/rplots/260513_OTX2_flow_expression_single_cells.pdf", width = 5, height = 7, useDingbats = FALSE); print(p1); dev.off()
pdf("./005_single_clone_analysis/rplots/260513_OTX2_flow_expression_clones.pdf",       width = 3, height = 4, useDingbats = FALSE); print(p2); dev.off()

# ---------------------------------------------------------------------
# 6. Flow cytometry vs RNA-seq correlation (batch 1 only)
# ---------------------------------------------------------------------
# Requires `rnaseq.rep.df` from:
# /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/003_RNA-seq_v2_0.R

corr.df <- rep.1$clone_summary %>%
  select(clone, flow.mean = exp.norm.mean, flow.sd = exp.norm.sd) %>%
  inner_join(
    rnaseq.rep.df %>% select(clone, rnaseq.mean = exp.norm.mean, rnaseq.sd = exp.norm.sd),
    by = "clone"
  )

corr.plot <- plot_correlation(
  corr.df, x = flow.mean, y = rnaseq.mean, xerr = flow.sd, yerr = rnaseq.sd, label_col = clone,
  x_limits = c(-25, 225), y_limits = c(-25, 225),
  x_lab = "Relative OTX2 expression from flow cytometry (% of WT)",
  y_lab = "Relative OTX2 expression from RNA-seq (% of WT)",
  title = "Flow cytometry vs RNA-seq OTX2 expression",
  subtitle = "Means \u00b1 SD; both axes normalised to WT (100 %)"
)
print(corr.plot)

pdf("./005_single_clone_analysis/rplots/260513_OTX2_expression_flow_vs_rna_seq.pdf", width = 3, height = 4, useDingbats = FALSE)
print(corr.plot)
dev.off()

# ---------------------------------------------------------------------
# 7. Calibration curve: flow cytometry vs gene expression score
# ---------------------------------------------------------------------

gene.exp.df <- read.csv("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/134_architectures_scores.csv")

# Clones of interest and their equivalent SAVANA deletion notation
#
# clone savana_notation
# 12    DEL_2-4
# 13    DEL_2-5
# 1E    DEL_2-3
# 1F    DEL_2-7
# 23    DEL_4-5
# D2    DEL_1-4
# DE    DEL_1-3
# DF    DEL_1-7
# E2    DEL_3-4
# EF    DEL_3-7
# WT    WT

clone.lup.df <- data.frame(
  clone = c("12", "13", "1E", "1F", "23", "D2", "DE", "DF", "E2", "EF", "WT"),
  savana_notation_grammar = c("DEL_2-4", "DEL_2-5", "DEL_2-3", "DEL_2-7", "DEL_4-5",
                              "DEL_1-4", "DEL_1-3", "DEL_1-7", "DEL_3-4", "DEL_3-7", "WT")
)

gene.exp.clone.df <- gene.exp.df %>%
  filter(savana_notation_grammar %in% clone.lup.df$savana_notation_grammar) %>%
  left_join(clone.lup.df, by = "savana_notation_grammar") %>%
  transmute(
    clone, savana_notation_grammar,
    gene_exp_score_mean = (mean_expression_score_1_scaled + mean_expression_score_2_scaled) / 2,
    gene_exp_score_sd   = abs(mean_expression_score_1_scaled - mean_expression_score_2_scaled) / 2
  )

calib.df <- rep.1$clone_summary %>%
  left_join(gene.exp.clone.df, by = "clone") %>%
  filter(!is.na(gene_exp_score_mean)) %>%
  select(clone, savana_notation_grammar, gene_exp_score_mean, gene_exp_score_sd, exp.norm.mean, exp.norm.sd)

calib.plot <- plot_correlation(
  calib.df, x = gene_exp_score_mean, y = exp.norm.mean, xerr = gene_exp_score_sd, yerr = exp.norm.sd,
  label_col = savana_notation_grammar,
  x_limits = c(-0.2, 1.2), y_limits = c(-25, 125),
  x_lab = "Gene expression score",
  y_lab = "Relative OTX2 expression\nfrom flow cytometry(% of WT)",
  title = "loxP7 screen calibration",
  subtitle = "Means \u00b1 SD or range, n = 10 clones",
  eb_width = 0.01
)

pdf("./005_single_clone_analysis/rplots/260513_OTX2_flow_expression_calibration.pdf", width = 3, height = 4, useDingbats = FALSE)
print(calib.plot)
dev.off()

# ---------------------------------------------------------------------
# 8. Batch 2 - data from 01/09/26
# ---------------------------------------------------------------------

flow.2.df <- load_and_resample("./005_single_clone_analysis/flow_data/260901_OTX2_clones_single_cell_clean.csv")

# Add WT / DF columns from batch 1, as a control for the impact of attB
# landing-pad insertion.
flow.2.df <- flow.2.df %>% cbind(flow.df %>% select(WT_1, WT_2, DF_1, DF_2))

flow.long.df.2 <- normalise_single_cells(flow.2.df)
rep.2 <- summarise_replicates(flow.long.df.2)

p3 <- plot_single_cells(flow.long.df.2, "WT",
                        "OTX2 expression across clones",
                        "Single-cell distributions; normalised to WT (100 %) and HAP1 (0 %)")

p4 <- plot_replicate_means(rep.2$clone_summary, rep.2$rep_points, "HAP1",
                           "OTX2 expression across clones",
                           "Replicate means \u00b1 SD; normalised to WT (100 %) and HAP1 (0 %)",
                           x_max = 130)

pdf("./005_single_clone_analysis/rplots/260901_OTX2_flow_expression_single_cells.pdf", width = 5, height = 7, useDingbats = FALSE); print(p3); dev.off()
pdf("./005_single_clone_analysis/rplots/260901_OTX2_flow_expression_clones.pdf",       width = 3, height = 4, useDingbats = FALSE); print(p4); dev.off()

# ---------------------------------------------------------------------
# 9. Combine the two batches and export
# ---------------------------------------------------------------------
#
# Two tables are combined and exported here:
#   - `combined.clone.summary`   : per-clone summary (mean/SD/SEM) of
#                                  normalised expression.
#   - `combined.single.cell.df`  : the underlying single-cell normalised
#                                  values (`exp`, `exp.norm`) that the
#                                  summary was computed from.
#
# For each, before combining, we:
#   1. Confirm the two tables have the same columns, in the same order
#      (so bind_rows()/rbind() lines things up correctly).
#   2. Add a `batch` column so rows stay traceable to their source, since
#      some clones (e.g. WT, DF) appear in both batches as controls.

# -- 9a. Per-clone summary table (mean / SD / SEM) --------------------

stopifnot(identical(names(rep.1$clone_summary), names(rep.2$clone_summary)))

combined.clone.summary <- bind_rows(
  rep.1$clone_summary %>% mutate(clone = as.character(clone), batch = "260513"),
  rep.2$clone_summary %>% mutate(clone = as.character(clone), batch = "260901")
) %>%
  relocate(batch, .after = clone)

# Report (rather than silently drop) any clone that was measured in both
# batches, since duplicated clone/batch combinations would indicate a
# data problem.
dup.rows <- combined.clone.summary %>% count(clone, batch) %>% filter(n > 1)
if (nrow(dup.rows) > 0) {
  warning("Duplicate clone/batch combinations found:\n",
          paste(capture.output(print(dup.rows)), collapse = "\n"))
}

write.csv(
  combined.clone.summary,
  "./005_single_clone_analysis/OTX2_clone_expression_summary.csv",
  row.names = FALSE
)


# -- 9b. Underlying single-cell values, one column per condition -------
#
# `flow.long.df.1` / `flow.long.df.2` hold one row per single-cell
# observation, in long format (`clone`, `exp`, `exp.norm`). Here we pivot
# each back to wide format - one column per clone/replicate condition,
# holding its normalised (`exp.norm`) values - and combine the two
# batches side by side. Since every condition was resampled to the same
# number of observations (`n` in `load_and_resample()`), rows don't
# represent matched cells across conditions; `obs_id` is just a 1..n
# position within each condition, needed to pivot to wide format.
#
# Some conditions (e.g. "WT_1", "DF_1") are shared between batches (see
# section 8), so columns are suffixed with the batch to keep every
# column name unique before combining.

widen_single_cell <- function(long_df, batch_label) {
  long_df %>%
    mutate(clone = as.character(clone)) %>%
    group_by(clone) %>%
    mutate(obs_id = row_number()) %>%
    ungroup() %>%
    dplyr::select(obs_id, clone, exp.norm) %>%
    pivot_wider(names_from = clone, values_from = exp.norm) %>%
    rename_with(~ paste0(.x, "_", batch_label), -obs_id)
}

single.cell.wide.1 <- widen_single_cell(flow.long.df.1, "260513")
single.cell.wide.2 <- widen_single_cell(flow.long.df.2, "260901")

combined.single.cell.wide.df <- full_join(single.cell.wide.1, single.cell.wide.2, by = "obs_id")

# Confirm every condition column ended up with a unique name.
stopifnot(!any(duplicated(names(combined.single.cell.wide.df))))

write.csv(
  combined.single.cell.wide.df,
  "./005_single_clone_analysis/OTX2_single_cell_expression.csv",
  row.names = FALSE
)
