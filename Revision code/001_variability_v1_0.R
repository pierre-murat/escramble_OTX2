
#==================================================#
# R1C3. The gene expression scores in Figure 5 show strong deviations between experiments.
# Could the authors comment on the origin of these differences between experiments?
#==================================================#

#==============================≠≠#
# Set up environement

setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2")

packages = c('dplyr', 'tidyr', 'ggplot2')
for(p in packages){library(p, character.only = T)}

#==============================≠≠#
# Load loxp7 data and analyse expression scores

loxP7.score.df <- read.csv("./134_architectures_scores.csv")

# Reformat table

loxP7.score.conf.df <- loxP7.score.df %>%
  mutate(read.count.1 = read_count_very_dim_1+read_count_dim_1+read_count_bright_1+read_count_very_bright_1,
         read.count.2 = read_count_very_dim_2+read_count_dim_2+read_count_bright_2+read_count_very_bright_2) %>% 
  dplyr::select(genotype = savana_notation_grammar, read.count.1, read.count.2, score.1 = mean_expression_score_1_scaled, score.2 = mean_expression_score_2_scaled)

# Compute error associated with the gene expression scores
loxP7.score.conf.df <- loxP7.score.conf.df %>% mutate(read.count=read.count.1+read.count.2, score.mean=(score.1+score.2)/2, score.max = pmax(score.1, score.2, na.rm = T), diff=score.max-score.mean)

# Reproduce Supllementary Figure 12a

# Compute values for error bars

loxP7.score.conf.df <- loxP7.score.conf.df %>% 
  mutate(read.count.min = log10(pmin(read.count.1, read.count.2, na.rm=T)),
         read.count.max = log10(pmax(read.count.1, read.count.2, na.rm=T)),
         read.count.mean = (read.count.min+read.count.max)/2)

# Plot

loxP7.score.conf.plot <- loxP7.score.conf.df %>% ggplot(aes(x = diff*100, y = read.count.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = read.count.min, ymax = read.count.max), width = 0.05) +
  ylab("Read count [log10]") + xlab("Gene expression score difference") +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  # geom_vline(xintercept=25, linetype="dashed") +
  # geom_hline(yintercept=log10(15), linetype="dashed") +
  theme_bw() + theme(aspect.ratio=1, panel.grid.minor = element_blank())
loxP7.score.conf.plot

# Compute corellation between read count and score difference (for the 134 detected architectures)

cor.test(loxP7.score.conf.df$read.count.mean, loxP7.score.conf.df$diff*100, method = "spearman")
# rho = -0.457, p-value = 1.718e-06
# There is a significant negative correlation between read count and score difference between experiments









