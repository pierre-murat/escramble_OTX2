
#==================================================#
# OTX2 expression from RNA-seq data in selected clones
#==================================================#

# Data are generated using the Plasmidsaurus RNA-seq service

#================================#
# Set up environement

setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2")

packages = c('dplyr', 'tidyr', 'ggplot2', 'tibble', 'DESeq2', 'ggrepel', 'scales', 'patchwork', 'GGally')
for(p in packages) {
  library(p, character.only = TRUE, verbose = FALSE, quietly = TRUE)
  cat(p, "version:", as.character(packageVersion(p)), "\n")
}

#================================#
# Backup

cd /lustre/scratch126/gengen/projects_v2/escramble
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/003_RNA-seq_v2_0.R
git commit -m "Last update committed on $(date +'%D %r')"
git push

#==================================================#
# Preliminary analyses to assess consistency between replicates and batches
#==================================================#

#================================#
# Load data

# 260209

L4KG9K.count.df <- read.table("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260209_ep10_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/L4KG9K_results/L4KG9K-expression-matrix.tsv", header = TRUE, row.names = 1, sep = "\t")  %>% 
  rownames_to_column("gene_id")
L4KG9K.summary.df <- read.csv("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260209_ep10_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/L4KG9K_results/L4KG9K-gene-biotype-5plus_reads-summary.csv") %>%
  dplyr::rename(sample = X) %>% separate(sample, into = c("sample", "sample.name"), sep = "_") %>% separate(sample.name, into = c("sample.name", "rep"), sep = "-rep")
# Rename columns
colnames(L4KG9K.count.df) <- c("gene_id", "gene_name", "gene_biotype",
                               "WT.batch1_rep1_cpm", "WT.batch1_rep1_count", "WT.batch1_rep2_cpm", "WT.batch1_rep2_count", "WT.batch1_rep3_cpm", "WT.batch1_rep3_count",
                               "del_R5_2_A1_rep1_cpm", "del_R5_2_A1_rep1_count", "del_R5_2_A1_rep2_cpm", "del_R5_2_A1_rep2_count", "del_R5_2_A1_rep3_cpm", "del_R5_2_A1_rep3_count")

# 260224
# Resubmit WT to assess variability and batch effect

J76TNK.count.df <- read.table("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260224_ep10_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/J76TNK_results/J76TNK-expression-matrix.tsv", header = TRUE, row.names = 1, sep = "\t")  %>% 
  rownames_to_column("gene_id")
J76TNK.summary.df <- read.csv("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260224_ep10_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/J76TNK_results/J76TNK-gene-biotype-5plus_reads-summary.csv") %>%
  dplyr::rename(sample = X) %>% separate(sample, into = c("sample.name", "rep"), sep = "__") #%>% separate(sample.name, into = c("sample.name", "rep"), sep = "-rep")
# Rename columns
colnames(J76TNK.count.df) <- c("gene_id", "gene_name", "gene_biotype",
                               "WT.batch2_rep1_cpm", "WT.batch2_rep1_count", "WT.batch2_rep2_cpm", "WT.batch2_rep2_count", "WT.batch2_rep3_cpm", "WT.batch2_rep3_count",
                               "del_R5_2_C5_rep1_cpm", "del_R5_2_C5_rep1_count", "del_R5_2_C5_rep2_cpm", "del_R5_2_C5_rep2_count", "del_R5_2_C5_rep3_cpm", "del_R5_2_C5_rep3_count",
                               "del_R5_2_D6_rep1_cpm", "del_R5_2_D6_rep1_count", "del_R5_2_D6_rep2_cpm", "del_R5_2_D6_rep2_count", "del_R5_2_D6_rep3_cpm", "del_R5_2_D6_rep3_count")

# Aggregate all data
count.df <- L4KG9K.count.df %>% 
  full_join(J76TNK.count.df, by = c("gene_id", "gene_name", "gene_biotype"))

#================================#
# Compare expression profiles of WT cells submitted independently in triplicates

WT.count.df <- count.df %>% 
  dplyr::select(gene_id,
                WT.batch1_rep1_cpm, WT.batch1_rep2_cpm, WT.batch1_rep3_cpm,
                WT.batch2_rep1_cpm, WT.batch2_rep2_cpm, WT.batch2_rep3_cpm)

# Compute and plot a pearson correlation matrix

mat <- WT.count.df %>% dplyr::select(-gene_id) %>% as.matrix()
# log1p transformation
keep <- rowMeans(mat) > 1
mat.log <- log2(mat[keep, ] + 1)

cor.mat <- cor(mat.log, method = "pearson", use = "pairwise.complete.obs")

cor.df <- as.data.frame(cor.mat) %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") %>% 
  mutate(Var1 = gsub("_cpm", "", Var1), Var2 = gsub("_cpm", "", Var2))

cor.df %>% ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", correlation)), size = 3, col = "white") +
  scale_fill_gradient2(
    low = "#2C75FF",
    mid = "white",
    high = "#ED1C24",
    midpoint = 0.98,   # adjust depending on your correlation range
    limits = c(min(cor.df$correlation), 1)
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  coord_fixed() +
  labs(fill = "Pearson r", x = NULL, y = NULL)

# Compute parowise corelations

ggpairs(as.data.frame(mat.log)) +
  theme_minimal() +
  ggtitle("Pairwise comparison (log1p CPM)")

#================================#
# Compare OTX2 expression between WT and del_R5_2 clones

# Unnormalised

otx2.df <- count.df %>%
  filter(gene_name == "OTX2") %>%
  dplyr::select(gene_name, contains("_cpm")) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample",
    values_to = "CPM"
  ) %>% mutate(sample = gsub("_cpm", "", sample)) %>% separate(sample, into = c("sample", "rep"), sep = "_rep") %>% 
  mutate(CPM = as.numeric(CPM), sample = factor(sample, levels=c("WT.batch1", "WT.batch2", "del_R5_2_A1", "del_R5_2_C5", "del_R5_2_D6")))

a <- otx2.df %>% ggplot(aes(x=sample, y=CPM)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean, width = 0.6, geom = "crossbar", fill = "#ED1C24", color = "#ED1C24") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.1, color = "#ED1C24") +
  geom_jitter(width = 0.1, size = 3) +
  ylim(0,120) + 
  ylab("OTX2 expression (CPM)") + ggtitle("Unnormalised") +
  theme_bw() + theme(aspect.ratio=1.5, panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

# Normalised to WT and per batch

# Compute baseline per sequencing batches
WT.mean.batch1 <- otx2.df %>% filter(grepl("WT.batch1", sample)) %>% pull(CPM) %>% mean()
WT.mean.batch2 <- otx2.df %>% filter(grepl("WT.batch2", sample)) %>% pull(CPM) %>% mean()

otx2.norm.df <- otx2.df %>% 
  mutate(CPM = case_when(sample == "WT.batch1" ~ CPM / WT.mean.batch1 * 100,
                         sample == "WT.batch2" ~ CPM / WT.mean.batch2 * 100,
                         sample == "del_R5_2_A1" ~ CPM / WT.mean.batch1 * 100,
                         sample == "del_R5_2_C5" ~ CPM / WT.mean.batch2 * 100,
                         sample == "del_R5_2_D6" ~ CPM / WT.mean.batch2 * 100))

b <- otx2.norm.df %>% ggplot(aes(x=sample, y=CPM)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean, width = 0.6, geom = "crossbar", fill = "#ED1C24", color = "#ED1C24") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.1, color = "#ED1C24") +
  geom_jitter(width = 0.1, size = 3) +
  ylim(0,120) + 
  ylab("Relative OTX2 expression") + ggtitle("Normalised to WT (per batch)") +
  theme_bw() + theme(aspect.ratio=1.5, panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

a + b

# Compare to GAPDH (and other house-keeping genes)

# GAPDH.df <- count.df %>% filter(gene_name == "GAPDH") %>% dplyr::select(gene_name, contains("_cpm")) %>% pivot_longer(
#   cols = -gene_name, names_to = "sample", values_to = "CPM")
# GAPDH.mean.batch1 <- GAPDH.df %>% filter(grepl("WT.batch1", sample)) %>% pull(CPM) %>% mean()
# GAPDH.mean.batch2 <- GAPDH.df %>% filter(grepl("WT.batch2", sample)) %>% pull(CPM) %>% mean()
# 
# otx2.gapdh.norm.df <- otx2.df %>% 
#   mutate(CPM = case_when(sample == "WT.batch1" ~ CPM / GAPDH.mean.batch1,
#                          sample == "WT.batch2" ~ CPM / GAPDH.mean.batch2,
#                          sample == "del_R5_2_A1" ~ CPM / GAPDH.mean.batch1 ,
#                          sample == "del_R5_2_C5" ~ CPM / GAPDH.mean.batch2,
#                          sample == "del_R5_2_D6" ~ CPM / GAPDH.mean.batch2))
# 
# otx2.gapdh.norm.df %>% ggplot(aes(x=sample, y=CPM)) +
#   #geom_hline(yintercept = 100, linetype = "dashed", color = "black") +
#   #geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   stat_summary(fun = mean, width = 0.6, geom = "crossbar", fill = "#ED1C24", color = "#ED1C24") +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.1, color = "#ED1C24") +
#   geom_jitter(width = 0.1, size = 3) +
#   #ylim(0,120) + 
#   ylab("Relative OTX2 expression") + ggtitle("Normalised to GAPDH (per batch)") +
#   theme_bw() + theme(aspect.ratio=1.5, panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

GAPDH.df <- count.df %>%
  filter(gene_name == "GAPDH") %>%
  dplyr::select(gene_name, contains("_cpm")) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample",
    values_to = "CPM"
  ) %>% mutate(sample = gsub("_cpm", "", sample)) %>% separate(sample, into = c("sample", "rep"), sep = "_rep") %>% 
  mutate(CPM = as.numeric(CPM), sample = factor(sample, levels=c("WT.batch1", "WT.batch2", "del_R5_2_A1", "del_R5_2_C5", "del_R5_2_D6")))

GAPDH.plot <- GAPDH.df %>% ggplot(aes(x=sample, y=CPM)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = mean, width = 0.6, geom = "crossbar", fill = "#ED1C24", color = "#ED1C24") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.1, color = "#ED1C24") +
  geom_jitter(width = 0.1, size = 3) +
  #ylim(0,120) + 
  ylab("GAPDH expression (CPM)") + #ggtitle("Unnormalised") +
  theme_bw() + theme(aspect.ratio=1.5, panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

GAPDH.plot + ACTB.plot + tub1a.plot

# Assess which genes are greatly different in WT from the different batches

WT.count.df <- count.df %>%
  rowwise() %>%
  mutate(
    WT.batch1.cpm = mean(c(WT.batch1_rep1_cpm,
                           WT.batch1_rep2_cpm,
                           WT.batch1_rep3_cpm)),
    
    WT.batch2.cpm = mean(c(WT.batch2_rep1_cpm,
                           WT.batch2_rep2_cpm,
                           WT.batch2_rep3_cpm)),
    
    WT.batch.pval = t.test(
      c(WT.batch1_rep1_cpm,
        WT.batch1_rep2_cpm,
        WT.batch1_rep3_cpm),
      c(WT.batch2_rep1_cpm,
        WT.batch2_rep2_cpm,
        WT.batch2_rep3_cpm)
    )$p.value
  ) %>%
  ungroup() %>% dplyr::select(gene_id, gene_name, WT.batch1.cpm, WT.batch2.cpm, WT.batch.pval) %>% 
  mutate(WT.LFC = log2(WT.batch2.cpm/WT.batch1.cpm)) %>% mutate(WT.cpm = (WT.batch1.cpm+WT.batch2.cpm)/2) %>%
  filter(!is.na(WT.LFC) & is.finite(WT.LFC))


plot(WT.count.df$WT.LFC, -log10(WT.count.df$WT.batch.pval))
plot(log10(WT.count.df$WT.cpm), WT.count.df$WT.LFC)
temp <- WT.count.df %>% filter(WT.LFC < 2.5 & log10(WT.cpm) > 1)

gene.to.highlight <- c("OTX2", "GAPDH", "CD74", "IGHM", "RPL29P11")

WT.count.df %>%
  mutate(highlight = case_when(gene_name == "OTX2" ~ "OTX2", gene_name %in% gene.to.highlight ~ "highlight_other", TRUE ~ "other")) %>%
  ggplot(aes(x = WT.cpm, y = WT.LFC)) +
  geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
  geom_text_repel(data = subset(batch.plot.df, highlight != "other"),
                  aes(label = gene_name, color = highlight), size = 3, fontface = "bold") +
  scale_color_manual(values = c("OTX2" = "#ED1C24", "highlight_other" = "#1F78B4", "other" = "#A7A9AC")) +
  scale_size_manual(values = c("OTX2" = 2.5, "highlight_other" = 2, "other" = 1), guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1, "highlight_other" = 1, "other" = 0.6), guide = "none") +
  scale_x_log10() + ylim(-10, 10) +
  xlab("Gene expression (mean, CPM)") + ylab("L2FC") + ggtitle("WT batch2 vs WT batch1") +
  theme_bw() + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "none")

# Differential expression analysis with DESeq2

# Extract count columns
count.df <- count.df %>% dplyr::select(gene_id, contains("_count"))

# Convert to matrix
count.mat <- as.matrix(count.df[,-1])

# Set gene IDs as rownames
rownames(count.mat) <- count.df$gene_id

# Round counts (DESeq2 requires integers)
count.mat <- round(count.mat)

# Create sample information dataframe

# sample <- factor(c(
#   "WT.batch1", "WT.batch1", "WT.batch1",
#   "del_R5_2_A1", "del_R5_2_A1", "del_R5_2_A1",
#   "WT.batch2", "WT.batch2", "WT.batch2",
#   "del_R5_2_C5", "del_R5_2_C5", "del_R5_2_C5",
#   "del_R5_2_D6", "del_R5_2_D6", "del_R5_2_D6"
# ))
# 
# coldata <- data.frame(
#   row.names = colnames(count.mat),
#   sample = sample,
#   condition = c(rep("WT", 3), rep("del_R5_2", 3), rep("WT", 3), rep("del_R5_2", 6)),
#   batch = c(rep("A", 6), rep("B", 9))
# )

condition <- factor(c(
  "WT.batch1", "WT.batch1", "WT.batch1",
  "del_R5_2_D6", "del_R5_2_D6", "del_R5_2_D6"
))

count.mat <- count.mat %>% as.data.frame() %>%
  dplyr::select(WT.batch1_rep1_count, WT.batch1_rep2_count, WT.batch1_rep3_count,
                                         del_R5_2_D6_rep1_count, del_R5_2_D6_rep2_count, del_R5_2_D6_rep3_count) %>%
  as.matrix()

coldata <- data.frame(
  row.names = colnames(count.mat),
  condition = condition
)

# Create DESeq dataset

# dds <- DESeqDataSetFromMatrix(
#   countData = count.mat,
#   colData = coldata,
#   design = ~ batch + condition
# )
dds <- DESeqDataSetFromMatrix(
  countData = count.mat,
  colData = coldata,
  design = ~ condition
)
# Filter out lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results (main contrast: del_R5_2 vs WT)
res.condition <- results(dds, contrast = c("condition", "del_R5_2_D6", "WT.batch1"))

# Add gene names to results
res.condition.df <- as.data.frame(res.condition) %>% rownames_to_column("gene_id") %>% left_join(L4KG9K.count.df %>% dplyr::select(gene_id, gene_name), by = "gene_id")
res.condition.df %>% filter(gene_name == "OTX2")

# Base mean expression plot
c <- res.condition.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
  scale_color_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC")) +
  scale_size_manual(values = c("OTX2" = 2, "other" = 1), guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1, "other" = 0.7), guide = "none") +
  geom_text_repel(data = subset(res.df, gene_name == "OTX2"), aes(label = gene_name),
                  color = "#ED1C24", size = 3, fontface = "bold", nudge_x = 1, hjust = 1, direction = "y", segment.size = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  ylim(-10, 10) +
  xlab("Gene expression (baseMean, CPM)") + ylab("L2FC") + ggtitle("del_R5_2 all vs WT") +
  theme_bw() + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "none")

# Volcano plot
d <- res.condition.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
  scale_color_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC")) +
  scale_size_manual(values = c("OTX2" = 2, "other" = 1), guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1, "other" = 0.7), guide = "none") +
  xlim(-10, 10) +
  geom_text_repel(data = subset(res.condition.df, gene_name == "OTX2"), aes(label = gene_name),
    color = "#ED1C24", size = 3, fontface = "bold", nudge_x = -1, hjust = 1, direction = "y", segment.size = 0) +
  xlab("L2FC") + ylab("P-value (-log10, adjusted)") + ggtitle("del_R5_2_D6 vs WT batch1") +
  theme_bw() + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "none")

# Combine all plots
c + d

# Pull affected genes with very significant adjusted p-values (padj < 1e-40)
affected.genes <- res.condition.df %>% filter(-log10(padj) > 40) %>% pull(gene_name)

# Sample-level Model

dds.sample <-  DESeqDataSetFromMatrix(
  countData = count.mat,
  colData = coldata,
  design = ~ sample
)
# Run DEseq2
dds.sample <- DESeq(dds.sample)

res.WT.df <- results(dds.sample, contrast = c("sample", "WT.batch2", "WT.batch1"))
res.WT.df <- as.data.frame(res.WT.df) %>% rownames_to_column("gene_id") %>% left_join(L4KG9K.count.df %>% dplyr::select(gene_id, gene_name), by = "gene_id")

res.del_R5_2_A1.df <- results(dds.sample, contrast = c("sample", "del_R5_2_A1", "WT.batch1"))
res.del_R5_2_A1.df <- as.data.frame(res.del_R5_2_A1.df) %>% rownames_to_column("gene_id") %>% left_join(L4KG9K.count.df %>% dplyr::select(gene_id, gene_name), by = "gene_id")

res.del_R5_2_C5.df <- results(dds.sample, contrast = c("sample", "del_R5_2_C5", "WT.batch2"))
res.del_R5_2_C5.df <- as.data.frame(res.del_R5_2_C5.df) %>% rownames_to_column("gene_id") %>% left_join(L4KG9K.count.df %>% dplyr::select(gene_id, gene_name), by = "gene_id")

res.del_R5_2_D6.df <- results(dds.sample, contrast = c("sample", "del_R5_2_D6", "WT.batch2"))
res.del_R5_2_D6.df <- as.data.frame(res.del_R5_2_D6.df) %>% rownames_to_column("gene_id") %>% left_join(L4KG9K.count.df %>% dplyr::select(gene_id, gene_name), by = "gene_id")

# Base mean expression plot
e <- res.del_R5_2_A1.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
  scale_color_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC")) +
  scale_size_manual(values = c("OTX2" = 2, "other" = 1), guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1, "other" = 0.7), guide = "none") +
  geom_text_repel(data = subset(res.df, gene_name == "OTX2"), aes(label = gene_name),
                  color = "#ED1C24", size = 3, fontface = "bold", nudge_x = 1, hjust = 1, direction = "y", segment.size = 0) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  ylim(-10, 10) +
  xlab("Gene expression (baseMean, CPM)") + ylab("L2FC") + ggtitle("WT_batch1 vs WT_batch2") +
  theme_bw() + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "none")

# Volcano plot
f <- res.del_R5_2_A1.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = highlight, size = highlight, alpha = highlight)) +
  scale_color_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC")) +
  scale_size_manual(values = c("OTX2" = 2, "other" = 1), guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1, "other" = 0.7), guide = "none") +
  xlim(-10, 10) +
  geom_text_repel(data = subset(res.df, gene_name == "OTX2"), aes(label = gene_name),
                  color = "#ED1C24", size = 3, fontface = "bold", nudge_x = -1, hjust = 1, direction = "y", segment.size = 0) +
  xlab("L2FC") + ylab("P-value (-log10, adjusted)") + ggtitle("WT_batch1 vs WT_batch2") +
  theme_bw() + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "none")

# Combine all plots
e + f

#=============================#
# Extract OTX2 normalised count
# gene id: ENSG00000165588

res.condition["ENSG00000165588",]

vsd <- vst(dds.sample)
plotCounts(dds.sample, gene="ENSG00000165588", intgroup=c("sample"))

plotPCA(vsd, intgroup=c("sample"))

#==============================#
# limma-voom workflow

library(edgeR)
library(limma)

# Create DGEList
dge <- DGEList(counts = count.mat)

# Filter low expression
keep <- filterByExpr(dge, group = coldata$condition)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalise
dge <- calcNormFactors(dge)

# Build design matrix
coldata$batch <- factor(coldata$batch)
coldata$condition <- factor(coldata$condition)
design <- model.matrix(~ batch + condition, data=coldata)
design

# Voom transformation
v <- voom(dge, design, plot=TRUE)

# Fit linear Model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract DEL vs WT
topTable(fit, coef="conditionWT")
fit$coefficients["ENSG00000165588", ]

# Extract batch effect
topTable(fit, coef="batchB")

otx2_expr <- v$E["ENSG00000165588", ]
boxplot(otx2_expr ~ coldata$condition)
boxplot(otx2_expr ~ coldata$sample)

plot.df <- data.frame(
  expression = otx2_expr,
  sample = coldata$sample,
  condition = coldata$condition,
  batch = coldata$batch
)
plot.df$condition <- factor(plot.df$condition, levels = c("WT", "del_R5_2"))

plot.df.norm <- plot.df %>%
  group_by(batch) %>%
  mutate(
    WT_mean_batch = mean(expression[condition == "WT"]),
    expr_norm = expression / WT_mean_batch
  ) %>%
  ungroup()

plot.df %>% ggplot(aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = batch),
              width = 0.15,
              size = 3) +
  scale_fill_manual(values = c("WT" = "#4C72B0", "del_R5_2" = "#DD8452")) +
  theme_classic(base_size = 14) +
  labs(
    title = "OTX2 expression (voom normalized)",
    y = "log2 expression",
    x = ""
  ) + theme_bw() +
  theme(aspect.ratio = 2, legend.position = "right")

plot.df %>% ggplot(aes(x = sample, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = batch),
              width = 0.15,
              size = 3) +
  scale_fill_manual(values = c("WT" = "#4C72B0", "del_R5_2" = "#DD8452")) +
  theme_classic(base_size = 14) +
  labs(
    title = "OTX2 expression (voom normalised)",
    y = "log2 expression",
    x = ""
  ) + theme_bw() +
  theme(aspect.ratio = 2, legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

plot.df.norm %>% ggplot(aes(x = sample, y = expr_norm, fill = condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = batch),
              width = 0.15,
              size = 3) +
  scale_fill_manual(values = c("WT" = "#4C72B0", "del_R5_2" = "#DD8452")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Relative OTX2 expression (voom normalised)",
    y = "log2 expression",
    x = ""
  ) + ylim(0,1.2) + theme_bw() +
  theme(aspect.ratio = 2, legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

#==================================================#
# OTX2 scramble clones
#==================================================#

#================================#
# Load data

# 260514

count.44K5D7.df <- read.table("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260514_gg6_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/44K5D7-expression-matrix.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)  %>% 
  rownames_to_column("gene_id") %>% 
  dplyr::select(1:3, contains("cpm"))
sample.44K5D7.df <- read.csv("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260514_gg6_pm23_OTX2_scramble_plasmidsaurus_RNA_seq/44K5D7-samplesheet.csv") %>%
  mutate(sample = gsub("sample_", "", sample), exp = paste(clone, well, sep = "_")) %>% 
  mutate(sample.id = paste("44K5D7", sample, "cpm", sep = "_")) %>% dplyr::select(exp, sample.id)

# Rename columns

count.44K5D7.df <- count.44K5D7.df %>% 
  rename_with(
    ~ sample.44K5D7.df$exp[
      match(.x, sample.44K5D7.df$sample.id)
    ],
    -c(gene_id, gene_name, gene_biotype)
  )

# Use WT samples from first dataset as the new samples show too much variability
# We also filter out outliers

WT.count.df <- L4KG9K.count.df %>% 
  dplyr::select(gene_id, WT_1 = WT.batch1_rep1_cpm, WT_2 = WT.batch1_rep2_cpm, WT_3 = WT.batch1_rep3_cpm)

count.44K5D7.df <- count.44K5D7.df %>% 
  dplyr::select(-WT_1_1, -WT_1_2, -WT_2_1, -WT_2_2, -E2_H2) %>% 
  left_join(WT.count.df, by = "gene_id")

# Compare OTX2 expression between WT and scrambled clones

otx2.44K5D7.df <- count.44K5D7.df %>%
  filter(gene_name == "OTX2") %>%
  dplyr::select(-gene_id, -gene_biotype) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample",
    values_to = "CPM"
  ) %>% mutate(sample = gsub("_cpm", "", sample))


otx2.44K5D7.df %>% ggplot(aes(x=sample, y=CPM)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black") +
  stat_summary(fun = mean, width = 0.6, geom = "crossbar", fill = "#ED1C24", color = "#ED1C24") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.1, color = "#ED1C24") +
  geom_jitter(width = 0.1, size = 3) +
  ylab("OTX2 expression (CPM)") + ggtitle("Unnormalised") +
  theme_bw() + theme(aspect.ratio=0.5, panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

# Normalised to WT and per batch

# Compute baseline

WT.cpm <- otx2.44K5D7.df %>% filter(grepl("WT", sample)) %>% pull(CPM) %>% median(na.rm = T)

otx2.44K5D7.norm.df <- otx2.44K5D7.df %>% 
  mutate(exp = (CPM / WT.cpm)*100)

# Summarise RNA-seq data by replicate and clone

rnaseq.rep.df <- otx2.44K5D7.norm.df %>%
  separate(col = sample, into = c("clone", "rep"), sep = "_", extra = "drop", fill = "right") %>%
  group_by(clone) %>%
  summarise(exp.norm.mean = mean(exp, na.rm = TRUE),
            exp.norm.sd   = sd(exp,   na.rm = TRUE),
            exp.norm.sem  = sd(exp,   na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  bind_rows(tibble(clone = "HAP1", exp.norm.mean = NA_real_, exp.norm.sd = NA_real_, exp.norm.sem = NA_real_)) %>%
  mutate(clone = factor(clone, levels = levels(flow.long.rep.df$clone)))

# Load flow.long.rep.df from /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/005_single_clone_analysis_v1_0.R
# To follow clone order and colour scheme

rnaseq.ind.df <- otx2.44K5D7.norm.df %>%
  separate(col = sample, into = c("clone", "rep"), sep = "_", extra = "drop", fill = "right") %>%
  bind_rows(tibble(clone = "HAP1", rep = NA_character_, exp = NA_real_)) %>%
  mutate(clone = factor(clone, levels = levels(flow.long.rep.df$clone)))

# Colour palette derived from flow.long.rep.df levels which include HAP1
cols.rnaseq <- setNames(
  ifelse(grepl("WT", levels(flow.long.rep.df$clone)), "#ED1C24", "#E6E7E8"),
  levels(flow.long.rep.df$clone)
)

# Plot

p.rnaseq <- rnaseq.rep.df %>%
  ggplot(aes(x = exp.norm.mean, y = clone)) +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_errorbarh(aes(xmin = exp.norm.mean - exp.norm.sd, xmax = exp.norm.mean + exp.norm.sd),
                 height = 0.35, linewidth = 0.6, na.rm = TRUE, colour = "grey30") +
  geom_point(data = rnaseq.ind.df, aes(x = exp, fill = clone),
             shape = 21, size = 2.5, colour = "grey20", stroke = 0.4) +
  geom_point(aes(fill = clone), shape = 124, size = 7, colour = "black") +
  scale_fill_manual(values = cols.rnaseq, guide = "none") +
  scale_x_continuous(limits = c(0, 225), breaks = seq(0, 225, 50), expand = c(0, 0)) +
  labs(x = "Relative OTX2 expression (% of WT)", y = NULL,
       title = "OTX2 expression across clones (RNA-seq)",
       subtitle = "Replicate means ± SD; normalised to WT median CPM (100 %)") +
  theme_bw(base_size = 11) +
  theme(aspect.ratio = 1.4,
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 9),
        axis.text = element_text(size = 9),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey30"))
print(p.rnaseq)

pdf("./005_single_clone_analysis/rplots/44K5D7_OTX2_rna_seq_expression_clones.pdf", width=3, height=4, useDingbats=FALSE)
p.rnaseq
dev.off()

#============================#   
# Calibration curve

# Load original gene expression scores

gene.exp.df <- read.csv("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/134_architectures_scores.csv")

# Recover architectures of interest

# We focus on 

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
  clone = c("12","13","1E","1F","23","D2","DE","DF","E2","EF","WT"),
  savana_notation_grammar = c("DEL_2-4","DEL_2-5","DEL_2-3","DEL_2-7","DEL_4-5","DEL_1-4","DEL_1-3","DEL_1-7","DEL_3-4","DEL_3-7","WT")
)

gene.exp.clone.df <- gene.exp.df %>% 
  filter(savana_notation_grammar %in% clone.lup.df$savana_notation_grammar) %>% 
  left_join(clone.lup.df, by = c("savana_notation_grammar")) %>%
  dplyr::select(clone, savana_notation_grammar, mean_expression_score_1_scaled, mean_expression_score_2_scaled) %>% 
  mutate(gene_exp_score_mean = (mean_expression_score_1_scaled + mean_expression_score_2_scaled)/2,
         gene_exp_score_sd = abs(mean_expression_score_1_scaled - mean_expression_score_2_scaled)/2) %>% 
  dplyr::select(-mean_expression_score_1_scaled, -mean_expression_score_2_scaled)

# Combine with RNA-seq data
calib.rnaseq.df <- rnaseq.rep.df %>%
  left_join(gene.exp.clone.df, by = "clone") %>%
  filter(!is.na(gene_exp_score_mean)) %>%
  dplyr::select(clone, savana_notation_grammar, gene_exp_score_mean, gene_exp_score_sd, exp.norm.mean, exp.norm.sd)

# Fit linear regression
lm.fit.rnaseq <- lm(exp.norm.mean ~ gene_exp_score_mean, data = calib.rnaseq.df)
lm.r2.rnaseq  <- summary(lm.fit.rnaseq)$r.squared
lm.p.rnaseq   <- summary(lm.fit.rnaseq)$coefficients[2, 4]
lm.lab.rnaseq <- sprintf("R² = %.3f\np = %.3e", lm.r2.rnaseq, lm.p.rnaseq)

# Prepare plot
calib.rnaseq.plot <- calib.rnaseq.df %>%
  ggplot(aes(x = gene_exp_score_mean, y = exp.norm.mean)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
              colour = "grey30", fill = "grey80", alpha = 0.4) +
  geom_errorbar(aes(ymin = exp.norm.mean - exp.norm.sd,
                    ymax = exp.norm.mean + exp.norm.sd),
                width = 0.01, linewidth = 0.4, colour = "grey40") +
  geom_errorbarh(aes(xmin = gene_exp_score_mean - gene_exp_score_sd,
                     xmax = gene_exp_score_mean + gene_exp_score_sd),
                 height = 2, linewidth = 0.4, colour = "grey40") +
  geom_point(shape = 21, size = 3, fill = "#ED1C24", colour = "grey20", stroke = 0.4) +
  geom_text_repel(aes(label = savana_notation_grammar), size = 3, colour = "grey20",
                  box.padding = 0.4, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, label = lm.lab.rnaseq,
           hjust = -0.1, vjust = 1.3, size = 3.2, colour = "grey20", fontface = "italic") +
  scale_x_continuous(limits = c(-0.2, 1.2), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-25, 190), breaks = seq(0, 200, 50), expand = c(0, 0)) +
  labs(x = "Gene expression score",
       y = "Relative OTX2 expression\nfrom RNA-seq (% of WT)",
       title = "loxP7 screen calibration",
       subtitle = "Means ± SD or range, n = 10 clones") +
  theme_bw(base_size = 11) +
  theme(aspect.ratio     = 1,
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(colour = "grey40", size = 9),
        axis.text        = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "grey30"))

pdf("./005_single_clone_analysis/rplots/260513_OTX2_rnaseq_expression_calibration.pdf", width = 3, height = 4, useDingbats = FALSE)
calib.rnaseq.plot
dev.off()

#============================#   
# DE-seq analysis with deltaR5-2

# Prepare count matrix
count.mat <- count.44K5D7.df %>%
  dplyr::select(gene_id, WT_1, WT_2, WT_3, `23_A1`, `23_C5`, `23_D6`) %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round()

# Sample metadata
coldata <- data.frame(
  row.names = colnames(count.mat),
  condition = factor(c("WT", "WT", "WT", "del_23", "del_23", "del_23"),
                     levels = c("WT", "del_23"))
)

# Create DESeq2 dataset, filter & run
dds <- DESeqDataSetFromMatrix(countData = count.mat, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Extract results (del_23 vs WT
res.23.df <- results(dds, contrast = c("condition", "del_23", "WT")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(count.44K5D7.df %>% dplyr::select(gene_id, gene_name, gene_biotype), by = "gene_id")

# Check OTX2
res.23.df %>% filter(gene_name == "OTX2")

# MA plot
p.ma <- res.23.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = highlight, size = highlight, alpha = highlight)) +
  scale_colour_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC"), guide = "none") +
  scale_size_manual(values  = c("OTX2" = 2,          "other" = 0.8),       guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1,          "other" = 0.5),       guide = "none") +
  geom_text_repel(data = filter(res.23.df, gene_name == "OTX2"),
                  aes(label = gene_name), colour = "#ED1C24", size = 3,
                  fontface = "bold", nudge_x = 0.3, segment.size = 0.3) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ylim(-10, 10) +
  labs(x = "Mean normalised expression (baseMean)",
       y = "Log2 fold change",
       title = "deltaR5-2 vs WT; DESeq2 analysis",
       subtitle = "MA plot") +
  theme_bw(base_size = 11) +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 9))

# Volcano plot
p.volcano <- res.23.df %>%
  mutate(highlight = ifelse(gene_name == "OTX2", "OTX2", "other")) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = highlight, size = highlight, alpha = highlight)) +
  scale_colour_manual(values = c("OTX2" = "#ED1C24", "other" = "#A7A9AC"), guide = "none") +
  scale_size_manual(values  = c("OTX2" = 2,          "other" = 0.8),       guide = "none") +
  scale_alpha_manual(values = c("OTX2" = 1,          "other" = 0.5),       guide = "none") +
  geom_text_repel(data = filter(res.23.df, gene_name == "OTX2"),
                  aes(label = gene_name), colour = "#ED1C24", size = 3,
                  fontface = "bold", nudge_x = -1, segment.size = 0.3) +
  xlim(-10, 10) +
  labs(x = "Log2 fold change",
       y = "Adjusted p-value (-log10)",
       subtitle = "Volcano plot") +
  theme_bw(base_size = 11) +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 9))

# Combine & save
p.ma + p.volcano

pdf("./005_single_clone_analysis/rplots/44K5D7_OTX2_deltaR52_clones_DESeq2.pdf", width = 8, height = 4, useDingbats = FALSE)
p.ma + p.volcano
dev.off()














