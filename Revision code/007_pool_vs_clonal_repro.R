
#===================================================#
#  Screen reproducibility between clonal and pool-based approaches
#===================================================#

# Load packages

packages = c('readxl', 'dplyr', 'tidyr', 'ggplot2', 'patchwork', 'ggrepel')
for(p in packages) {
  library(p, character.only = TRUE)
  cat(p, "version:", as.character(packageVersion(p)), "\n")
}

# Backup

cd /lustre/scratch126/gengen/projects_v2/escramble
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/007_pool_vs_clonal_repro.R
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/007_pool_vs_clonal_repro/rplots
git commit -m "Last update committed on $(date +'%D %r')"
git push

# Set working directories

setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes")

# Load data

escramble.data.path <- "./Supplementary_table.xlsx"
escramble.data <- excel_sheets(escramble.data.path)

loxP6.clonal.df <- read_excel(escramble.data.path, sheet = "loxPsym_6_clonal_cas9")
loxP6.pool.df <- read_excel(escramble.data.path, sheet = "loxPsym_6_pool_cas9")

# Compute gene expression score for all variants

loxP6.clonal.score.df <- loxP6.clonal.df %>%
  group_by(name, replicate) %>%
  mutate(scale = read_percentage/sum(read_percentage)) %>% 
  dplyr::select(name, gate, replicate, scale) %>% ungroup() %>% 
  pivot_wider(names_from = gate, values_from = scale) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(score = (-1*dark)+(-0.5*dim)+(0.5*medium)+(1*bright)) %>% 
  ungroup() %>% 
  group_by(name) %>% summarise(score.clonal = mean(score), range.clonal = max(score)-min(score)) %>% 
  mutate(name = case_when(name == "no_sv_evidence" ~ "WT", TRUE ~ name))

loxP6.pool.score.df <- loxP6.pool.df %>%
  group_by(name, replicate) %>%
  mutate(scale = read_percentage/sum(read_percentage)) %>% 
  dplyr::select(name, gate, replicate, scale) %>% ungroup() %>% 
  pivot_wider(names_from = gate, values_from = scale) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(score = (-1*dark)+(-0.5*dim)+(0.5*medium)+(1*bright)) %>% 
  ungroup() %>% 
  group_by(name) %>% summarise(score.pool = mean(score), range.pool = max(score)-min(score)) %>% 
  mutate(name = case_when(name == "no sv evidence" ~ "WT", TRUE ~ name))

# Combine

loxP6.score.df <- loxP6.clonal.score.df %>% full_join(loxP6.pool.score.df, by = "name")

# Select consistent architectures

loxP6.score.df <- loxP6.score.df %>% filter(name %in% c("WT", "∆R123", "∆R1234", "iR1234", "∆R12345", "∆R34"))

# Normalise scores

wt.clonal <- loxP6.score.df$score.clonal[loxP6.score.df$name == "WT"]
del.clonal <- loxP6.score.df$score.clonal[loxP6.score.df$name == "∆R12345"]

wt.pool <- loxP6.score.df$score.pool[loxP6.score.df$name == "WT"]
del.pool <- loxP6.score.df$score.pool[loxP6.score.df$name == "∆R12345"]

clonal.factor <- 100 / (wt.clonal - del.clonal)
pool.factor   <- 100 / (wt.pool - del.pool)

loxP6.score.df <- loxP6.score.df %>%
  mutate(
    score.clonal.norm = (score.clonal - del.clonal) * clonal.factor,
    range.clonal.norm = range.clonal * clonal.factor,
    score.pool.norm   = (score.pool - del.pool) * pool.factor,
    range.pool.norm   = range.pool * pool.factor
  )

# Fit linear regression
lm.fit <- lm(score.pool.norm ~ score.clonal.norm, data = loxP6.score.df)
lm.r2  <- summary(lm.fit)$r.squared
lm.p   <- summary(lm.fit)$coefficients[2, 4]
lm.lab <- sprintf("R² = %.3f\np = %.3e", lm.r2, lm.p)

# Plot

loxP6.score.plot <- loxP6.score.df %>%
  ggplot(aes(x = score.pool.norm, y = score.clonal.norm)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
              colour = "grey30", fill = "grey80", alpha = 0.4) +
  geom_errorbar(aes(ymin = score.clonal.norm - range.clonal.norm, ymax = score.clonal.norm + range.clonal.norm),
                width = 2, linewidth = 0.4, colour = "grey40") +
  geom_errorbarh(aes(xmin = score.pool.norm - range.pool.norm, xmax = score.pool.norm + range.pool.norm),
                 height = 2, linewidth = 0.4, colour = "grey40") +
  geom_point(shape = 21, size = 3, fill = "#ED1C24", colour = "grey20", stroke = 0.4) +
  geom_text_repel(aes(label = name), size = 3, colour = "grey20",
                  box.padding = 0.4, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, label = lm.lab,
           hjust = -0.1, vjust = 1.3, size = 3.2, colour = "grey20", fontface = "italic") +
  scale_x_continuous(limits = c(-50, 175), breaks = seq(0, 150, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-100, 350), breaks = seq(0, 300, 50), expand = c(0, 0)) +
  labs(x = "Relative gene expression\n[pooled transfection experiment]",
       y = "Relative gene expression\n[clonal experiment]",
       subtitle = "Means ± ranges") +
  theme_bw(base_size = 11) +
  theme(aspect.ratio     = 1,
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        axis.text        = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey30"))
loxP6.score.plot

# Save plot

pdf("./Revision_R_codes/007_pool_vs_clonal_repro/rplots/pool.vs.clonal.repro.pdf", width=3, height=4, useDingbats=FALSE)
loxP6.score.plot
dev.off()














