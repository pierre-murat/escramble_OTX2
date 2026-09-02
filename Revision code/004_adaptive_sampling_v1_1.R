
#=================================================#
# Assess the use of adaptive sampling for 
# genotyping the OTX2 locus
#=================================================#

# Set working directory
setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling")

# Load library
library(tidyverse)
library(rtracklayer)
library(patchwork)

# Backup

cd /lustre/scratch126/gengen/projects_v2/escramble
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling_v1_1.R
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/scripts
git commit -m "Last update committed on $(date +'%D %r')"
git push

#=================================================#
# Prepare files for adaptive sampling

# Information about adaptive sampling for ONT sequencing can be found here:
# https://nanoporetech.com/document/adaptive-sampling

# ONT advises uploading no more than 125 Mb of non-indexed reference (FASTA).
# Above this size, the MinION will experience heavy delays in the adaptive sampling decision.
# So instead, of uploading the full human genome, we use only chromosome 14, which contains the OTX2 locus.
# This chromosome is 107 Mb long, so it is suitable for adaptive sampling.

# This is a good compromise for:
# - Reducing the time for the adaptive sampling decision
# - Keeping a good enough decission accuracy for rejecting reads with seconfary alignments to other chromosomes.

#####
# Prepare the reference FASTA file for adaptive sampling

# bash
# cd /lustre/scratch125/gengen/teams_v2/parts/pm23/genomes/hg38
# module load HGI/softpack/groups/escramble/eSCRAMBLE/9
# samtools faidx hg38.fa.gz chr14 > /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/fasta/chr14.fa
# This is 109 Mb, which follow ONT recommendation

#####
# Prepare the ROI BED file

# We prepare a bed file considering strand directionality without any buffer

# ROI: chr14:56,795,000-56,920,000 (125 kb)

# Strand directionalitiy is defined by defining two indepedent intervals for the forward and reverse strands,
# + strand interval - stop 25 kb before the end of the ROI: chr14:56,795,000-56,895,000
# - strand interval - start 25 kb after the start of the ROI: chr14:56,820,000-56,920,000

# Create bed file

roi.bed <- data.frame(
  chrom = c("chr14", "chr14"),
  start = c(56795000, 56820000),
  end = c(56895000, 56920000),
  gene = c("OTX2_forward", "OTX2_reverse"),
  score = c(0, 0),
  strand = c("+", "-")
)

# Export bed file

write.table(roi.bed, file = "./bed/OTX2_ROI.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# The generated bed file can be checked at https://epi2me.nanoporetech.com/bed-bugs/
# BED file OK

# A mixture of WT, deltaR123, deltaR3 and deltaR45 were prep (adpated ligation kit) and sequenced with adaptive sampling using the above files.
# The resulting POD5 files are /lustre/scratch125/gengen/teams_v2/parts/sequencing/260427_pm23_OTX2_adaptive_sampling for two flow cells.

#=================================================#
# Basecalling

# Generate fastq files from the POD5 files using dorado basecaller.
# Using the fast model

# The shell script for basecalling is in
# /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/scripts/dorado_basecalling.sh

# Check that all bam files were generated correctly (need to run the script iteratively to increase memory allowance)
# 96 jobs
grep -rl "Successfully completed" *.log | wc -l # 96 out of 96 
grep -rl "job killed after reaching LSF memory usage limit" *.log | wc -l # 0 out of 96

#=================================================#
# Alignment

# Align bam files to the hg38 reference genome using dorado/minimap2 with the following script
# /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/scripts/dorado_alignment.sh

#=================================================#
# Merge all aligned bam file

# Merge all aligned bam files using samtools merge and create a bigwig file for visualisation in IGV using deeptools bamCoverage with the following script
# /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling/scripts/samtools_merge.sh

# Check that all bam files were generated correctly (need to run the script iteratively to increase memory allowance)
# 96 jobs
grep -rl "Successfully completed" *.log | wc -l # 96 out of 96 

#=================================================#
# Split bam files for reads that overlap and do not overlap with the ROI (OTX2 locus)

module load HGI/softpack/groups/escramble/eSCRAMBLE/9
echo -e "chr14\t56794999\t56920000" > ROI.bed
samtools view -b -L ROI.bed -U OTX2_as.other.sorted.bam OTX2_as.sorted.bam > OTX2_as.ROI.sorted.bam

# Count number of reads

# OTX2_as.ROI.sorted.bam, 74 reads
# OTX2_as.other.sorted.bam, 1836139 reads
# Total, 1,836,213 reads

#=================================================#
# Extract read statistics

samtools fastq OTX2_as.ROI.sorted.bam | seqkit stats -a > OTX2_as.ROI.read.stats.tsv
samtools fastq OTX2_as.other.sorted.bam | seqkit stats -a > OTX2_as.other.read.stats.tsv

# file  format  type  num_seqs       sum_len  min_len  avg_len  max_len   Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)
# ROI   FASTQ   DNA         72       271,296      193    3,768   48,190  415  1,585  4,045        0  9,338   64.15   34.94
# other FASTQ   DNA   1,560,686  734,250,179        5    470.5  259,874  385    484    558        0    519   57.96   30.61

#=================================================#
# Assess read length distribution

samtools view OTX2_as.ROI.sorted.bam | awk '{print length($10)}' > OTX2_as.ROI.read.lengths.txt
samtools view OTX2_as.other.sorted.bam | awk '{print length($10)}' > OTX2_as.other.read.lengths.txt

# Load in R

roi <- read_table("./bam_merged/OTX2_as.ROI.read.lengths.txt", col_names = "read_length")%>% mutate(sample = "ROI")
other <- read_table("./bam_merged/OTX2_as.other.read.lengths.txt", col_names = "read_length") %>% mutate(sample = "other")
df <- bind_rows(roi, other) %>% mutate(sample = factor(sample, levels = c("ROI", "other")))

# Compute summary statistics
stats <- df %>%
  group_by(sample) %>%
  summarise(
    N        = n(),
    Mean     = mean(read_length),
    Median   = median(read_length),
    SD       = sd(read_length),
    Min      = min(read_length),
    Max      = max(read_length),
    N50      = sort(read_length)[which(cumsum(sort(read_length)) >= sum(read_length) / 2)[1]],
    .groups  = "drop"
  )
print(stats)

# sample       N  Mean   edian    SD   Min    Max   N50
#    ROI       74 3681.   1555  6877.  193  48190  9338
#  other  1836139  401.    458   479.    1 259874   519

# Colour palette
pal <- c("ROI" = "#4DAF8D", "other" = "#8E6BBF")

# Wilcoxon test
wtest <- wilcox.test(read_length ~ sample, data = df)
p_label <- sprintf("P = %.3e", wtest$p.value)

# Plot
p1 <- ggplot(df, aes(x = sample, y = read_length, fill = sample)) +
  geom_violin(alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black") +
  geom_text(data = stats,
            aes(x = sample,
                y = max(df$read_length) * 2,
                label = paste0("N = ", scales::comma(N),
                               "\nMedian = ", scales::comma(Median),
                               "\nN50 = ", scales::comma(N50))),
            size = 3.5, vjust = 0) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = pal) +
  labs(
    title    = "Read length distribution",
    subtitle = p_label,
    x        = NULL,
    y        = "Read length (bp, log10 scale)",
    fill     = "Region"
  ) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio  = 1.5,
    legend.position = "none",
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

print(p1)

#=================================================#
# Assess read coverage

# Generate coverage tracks using the bamCoverage function fromo deeptools

cd /lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/004_adaptive_sampling

samtools index ./bam_merged/OTX2_as.ROI.sorted.bam
bamCoverage \
    --bam ./bam_merged/OTX2_as.ROI.sorted.bam \
    --outFileName ./bw/OTX2_as.ROI.1kb.bw \
    --outFileFormat bigwig \
    --binSize 1000 \
    --normalizeUsing None

samtools index ./bam_merged/OTX2_as.other.sorted.bam
bamCoverage \
    --bam ./bam_merged/OTX2_as.other.sorted.bam \
    --outFileName ./bw/OTX2_as.other.1kb.bw \
    --outFileFormat bigwig \
    --binSize 1000 \
    --normalizeUsing None

# Load bigWig files
roi.bw   <- import("./bw/OTX2_as.ROI.1kb.bw",   format = "bigWig")
other.bw <- import("./bw/OTX2_as.other.1kb.bw", format = "bigWig")

# Convert to data frames
df_roi   <- as.data.frame(roi.bw)   %>% mutate(sample = "ROI")
df_other <- as.data.frame(other.bw) %>% mutate(sample = "other")

cov_df <- bind_rows(df_roi, df_other) %>%
  mutate(sample = factor(sample, levels = c("ROI", "other"))) %>%
  filter(seqnames == "chr14", start >= 40000000, end <= 70000000)

# Compute summary statistics
cov_stats <- cov_df %>%
  group_by(sample) %>%
  summarise(
    N      = n(),
    Mean   = mean(score),
    Median = median(score),
    SD     = sd(score),
    Max    = max(score),
    .groups = "drop"
  )
print(cov_stats)

# sample          N   Mean.  Median     SD     Max
#    ROI       72    0.172        0  0.808       7
#  other    1525806   1.16        1   7.34    2263

# Wilcoxon test
wtest   <- wilcox.test(score ~ sample, data = cov_df)
p_label <- sprintf("P = %.3e", wtest$p.value)

# Plot
p2 <- cov_df %>%
  ggplot(aes(x = sample, y = score, fill = sample)) +
  geom_violin(alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, colour = "black") +
  geom_text(data = stats,
            aes(x      = sample,
                y      = max(cov_df$score) * 2,
                label  = paste0("N bins = ", scales::comma(N),
                                "\nMean = ",   round(Mean,   2),
                                "\nMedian = ", round(Median, 2))),
            size  = 3.5,
            vjust = 0) +
  scale_fill_manual(values = pal) +
  labs(
    title    = "Coverage distribution",
    subtitle = p_label,
    x        = NULL,
    y        = "Read coverage per 1 kb bin",
    fill     = "Sample"
  ) + ylim(0,8) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio  = 1.5,
    legend.position = "none",
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

print(p2)

#=================================================#
# Combine plots

p1 + p2


