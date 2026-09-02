
#==================================================#
# R2C1. The rationale behind the weights applied is not clear and seems somewhat arbitrary but has a big effect on how the data are interpreted.
#==================================================#

#==============================≠≠#
# Set up environement

setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2")

packages = c('dplyr', 'tidyr', 'ggplot2')
for(p in packages){library(p, character.only = T)}

#==============================≠≠#
# Load loxp7 high confidence data and compute expression scores with different weights

loxP7.score.df <- read.csv("./61_high_conf_architectures_scores.tsv", sep = "\t")
nrow(loxP7.score.df)
# 61 architectures

# Define function to compute score with different weights

# Original weights of (-1, -0.5, 0.5, 1)

compute_scores <- function(df, weights){
  w1 <- weights[1]
  w2 <- weights[2]
  w3 <- weights[3]
  w4 <- weights[4]
  df <- df %>%
    mutate(score.1 = (read_count_very_dim_1*w1 + read_count_dim_1*w2 + read_count_bright_1*w3 + read_count_very_bright_1*w4) / (read_count_very_dim_1 + read_count_dim_1 + read_count_bright_1 + read_count_very_bright_1),
           score.2 = (read_count_very_dim_2*w1 + read_count_dim_2*w2 + read_count_bright_2*w3 + read_count_very_bright_2*w4) / (read_count_very_dim_2 + read_count_dim_2 + read_count_bright_2 + read_count_very_bright_2)) %>%
    mutate(score.mean=(score.1+score.2)/2, score.max = pmax(score.1, score.2, na.rm = T), diff=score.max-score.mean)
  return(df)
}

# Compute scores with different weights

# We evaluated ten weight configurations spanning scale, sparsity, nonlinearity, dominance, and random sampling.

# Define weights
weights.list <- list(
  original          = c(-1, -0.5, 0.5, 1),
  strong_scaled     = c(-2, -1, 1, 2),
  medium            = c(-0.75, -0.25, 0.25, 0.75),
  weak              = c(-0.5, -0.1, 0.1, 0.5),
  zero_centered     = c(-1, 0, 0, 1),
  rank_only         = c(-1, -1, 1, 1),
  nonlinear         = c(-1, -0.2, 0.2, 1),
  dominant_left     = c(-5, -0.1, 0.1, 1),
  dominant_right    = c(-1, -0.1, 0.1, 5),
  compressed_middle = c(-1, -0.9, 0.9, 1)
)

# Compute normalized scores for all weights
scores.results <- lapply(weights.list, function(w) {
  df <- compute_scores(loxP7.score.df, w) %>% 
    dplyr::select(genotype = savana_notation_grammar, score.1, score.2)
  ref_del <- mean(df$score.1[df$genotype=="DEL_1-7"], na.rm=TRUE)
  ref_wt  <- mean(df$score.1[df$genotype=="WT"], na.rm=TRUE)
  df %>% mutate(across(score.1:score.2, ~ (. - ref_del)/(ref_wt-ref_del)*100))
})

# Prepare reference scores
ref <- scores.results$original %>% dplyr::select(genotype, score.1) %>% rename(score_ref = score.1)

# Long data for plotting
plot.df <- bind_rows(lapply(names(scores.results), function(w) {
  scores.results[[w]] %>%
    dplyr::select(genotype, score.1) %>% rename(score_alt=score.1) %>%
    left_join(ref, by="genotype") %>% mutate(weights=w)
})) %>% filter(weights!="original")

# Spearman correlations
cors <- plot.df %>% group_by(weights) %>%
  summarise(spearman=round(cor(score_ref, score_alt, method="spearman"),3), .groups="drop")

# Labels with weights
weights.labels <- tibble(weights=names(weights.list),
                         label=paste0(weights,"\n(", sapply(weights.list, paste, collapse=", "), ")"))
plot.df <- plot.df %>% left_join(weights.labels, by="weights")

# Summary for crosshair error bars (both axes)
summary.df <- bind_rows(lapply(names(scores.results), function(w) {
  df <- scores.results[[w]] %>% mutate(weights=w) %>%
    rowwise() %>%
    mutate(
      mean_alt = mean(c(score.1, score.2), na.rm=TRUE),
      diff_alt = abs(score.1 - score.2),
      mean_ref = mean(c(
        scores.results$original$score.1[scores.results$original$genotype==genotype],
        scores.results$original$score.2[scores.results$original$genotype==genotype]
      ), na.rm=TRUE),
      diff_ref = abs(
        diff(c(
          scores.results$original$score.1[scores.results$original$genotype==genotype],
          scores.results$original$score.2[scores.results$original$genotype==genotype]
        ))
      )
    ) %>%
    ungroup() %>%
    dplyr::select(genotype, weights, mean_alt, diff_alt, mean_ref, diff_ref)
  df
})) %>% filter(weights!="original") %>% left_join(weights.labels, by="weights")

# Final ggplot
ggplot(plot.df, aes(score_ref, score_alt)) +
  geom_point(size=1.5, alpha=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=0.8, color="red") +
  # geom_errorbar(data=summary.df, aes(x=mean_ref, ymin=mean_alt-diff_alt/2, ymax=mean_alt+diff_alt/2),
  #               inherit.aes=FALSE, color="black", linewidth=0.5) +
  # geom_errorbarh(data=summary.df, aes(y=mean_alt, xmin=mean_ref-diff_ref/2, xmax=mean_ref+diff_ref/2),
  #                inherit.aes=FALSE, color="black", linewidth=0.5) +
  # geom_point(data=summary.df, aes(x=mean_ref, y=mean_alt),
  #            inherit.aes=FALSE, color="black", size=2.5, shape=18) +
  geom_text(data=cors %>% left_join(weights.labels, by="weights"),
            aes(x=-Inf, y=Inf, label=paste0("ρ=", spearman)),
            inherit.aes=FALSE, hjust=-0.1, vjust=1.1, size=3.5) +
  facet_wrap(~ label, ncol=3) +
  labs(x="Mean score (original (-1, -0.5, 0.5, 1) weights)", y="Mean score (alternative weights)") +
  xlim(-30,130) + ylim(-30,130) +
  theme_bw() +
  theme(aspect.ratio=1,
        strip.background=element_rect(fill="grey95"),
        strip.text=element_text(size=9),
        panel.grid.minor=element_blank())

# We assessed the robustness of gene expression scores to different weighting schemes by computing normalized scores for 10 alternative weight sets,
# comparing each to the original weights. Spearman correlation coefficients were calculated for each weight set to quantify the concordance with the
# original scores. The correlations ranged from 0.897 (dominant_right) to 1 (strong_scaled), with most alternative weights showing very high concordance (>0.95).
# Specifically, medium and nonlinear weights yielded correlations above 0.98, while rank_only and compressed_middle were slightly lower but still highly correlated.

# The results demonstrate that the overall ranking and relative gene expression scores are highly robust to variations in the weighting scheme.
# Even extreme or asymmetric weight configurations (dominant_left, dominant_right, compressed_middle) maintain strong correlations with the original scores.
# This indicates that the choice of specific weights has minimal impact on the relative assessment of genotypes, supporting the reliability of the scoring method.
# In practice, this suggests that conclusions drawn from the original weighting scheme are stable and not sensitive to reasonable variations in weight selection.










