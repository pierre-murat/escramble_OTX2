
#==================================================#
# ATAC-seq analysis in selected clones
#==================================================#

# using the nf-core/atacseq pipeline for ATAC-seq data analysis
# https://github.com/nf-core/atacseq

# on selected clones
# WT, AD, R4, R5, R45

# Sequening data is available at
# /lustre/scratch125/gengen/teams_v2/parts/sequencing/260526_ep10_duck_pm23_otx2_actac_seq

#================================#
# Set up environement

setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2")

packages = c('dplyr', 'tidyr', 'ggplot2', 'scales', 'patchwork', 'rtracklayer', 'GenomicRanges', 'data.table', 'rstatix')
for(p in packages) {
  library(p, character.only = TRUE, verbose = FALSE, quietly = TRUE)
  cat(p, "version:", as.character(packageVersion(p)), "\n")
}

#================================#
# Backup

cd /lustre/scratch126/gengen/projects_v2/escramble
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq_v1_0.R
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts
git add ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/rplots
git commit -m "Last update committed on $(date +'%D %r')"
git push

#==================================================#
# Data conversion [cram to fastq]
#==================================================#

# Prepare lookup table for samples and their corresponding cram files

# Sample  Tag
# WT_1    24
# WT_2    25
# AD_1    3
# AD_2    23
# R4_1    28
# R4_2    2
# R5_1    26
# R5_2    27
# R45_1   4
# R45_2   1

cram.lookup.df <- cbind.data.frame(
  sample = c("WT", "WT", "AD", "AD", "R4", "R4", "R5", "R5", "R45", "R45"),
  rep = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2),
  plex = c(24, 25, 3, 23, 28, 2, 26, 27, 4, 1)
) %>% mutate(
  path.cram = paste0("/lustre/scratch125/gengen/teams_v2/parts/sequencing/260526_ep10_duck_pm23_otx2_actac_seq/52295/plex", plex, "/52295#", plex, ".cram")
)
# Export

write.table(cram.lookup.df,
            "./Revision_R_codes/006_atac_seq/scripts/cram_to_fastq_table.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Run
# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/cram_to_fastq.sh

# Check that all jobs have completed successfully
grep -rl "Successfully completed" *.out | wc -l # 10 out of 10 

# Count number of reads in each fastq file

for FQ in *.fastq.gz; do
BASENAME=$(basename "$FQ" .fastq.gz)
SAMPLE=$(echo "$BASENAME" | sed 's/_rep_.*//')
REP=$(echo "$BASENAME"   | sed 's/.*_rep_//' | cut -d'_' -f1)
TYPE=$(echo "$BASENAME"  | sed "s/${SAMPLE}_rep_${REP}_//")
N_READS=$(zcat "$FQ" | awk 'END {print NR/4}')
echo -e "${SAMPLE}\t${REP}\t${TYPE}\t${N_READS}"
done | sort -k1,1 -k2,2n | cat <(echo -e "sample\treplicate\tread_type\tn_reads") - > read_counts_summary.tsv

read_counts_summary.df <- read.table("./Revision_R_codes/006_atac_seq/fastq/read_counts_summary.tsv", header = TRUE, sep = "\t") %>%
  pivot_wider(names_from = read_type, values_from = n_reads)

#==================================================#
# Data processing with nf-core/atacseq pipeline
#==================================================#

# Prepare lookup table

nfcore.lookup.df <- cram.lookup.df %>% 
  dplyr::select(sample, rep) %>% 
  mutate(fastq_1 = paste0("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/fastq/", sample, "_rep_", rep, "_R1.fastq.gz"),
         fastq_2 = paste0("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/fastq/", sample, "_rep_", rep, "_R2.fastq.gz")
  ) %>% 
  dplyr::select(sample, fastq_1, fastq_2, replicate = rep)

write.table(nfcore.lookup.df,
            "./Revision_R_codes/006_atac_seq/scripts/atac_nfcore/nfcore.lookup.csv",
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# Check the sanger profile for nf-core/atacseq pipeline

nextflow config nf-core/atacseq -profile sanger > sanger.profile.txt

# Create config and sh scripts for running the pipeline

# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/nextflow.config
# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/nextflow_nfcore_atac.sh
# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/run_nextflow_nfcore_atac.sh

# The pipeline stopped because of memory quotas

# Remove work directories for cached/completed processes only, retaining failed ones

module load HGI/common/nextflow/23.10.0
module load singularityce-3.11.5/python-3.11.6

# List all previous run IDs
nextflow log
# want to keep lethal_rubens

# Preview first
nextflow clean -but lethal_rubens -dry-run
# Execute 
nextflow clean -but lethal_rubens -f


# Then clean everything except the most recent run, using its run name
nextflow clean -but <run_name> -dry-run   # preview first
nextflow clean -but <run_name> -f         # execute (-f to force)

# Resume

# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/nextflow_nfcore_atac_resume.sh
# ./001_007_OTX2/Manuscripts/OTX2/Revision_R_codes/006_atac_seq/scripts/atac_nfcore/run_nextflow_nfcore_atac_resume.sh

# Pipeline completed

# COMMENT: R4R5 deletion clones are not correct.
# We ignore these from the analysis

#==================================================#
# Data formating and genomic track overview
#==================================================#

# Normalise ATAC-seq signal by sequencing depth (CPM)

normalise_bw <- function(input_path, output_path) {
  bw <- import(input_path)
  bw.norm <- bw
  bw.norm$score <- bw$score * 1e6 / sum(bw$score)
  export(bw.norm, output_path)
  message("Exported: ", output_path)
}

normalise_bw("./Revision_R_codes/006_atac_seq/atac_nf/bwa/merged_replicate/bigwig/AD.mRp.clN.bigWig",
             "./Revision_R_codes/006_atac_seq/bw/AD.atac.cpm.bigWig")
normalise_bw("./Revision_R_codes/006_atac_seq/atac_nf/bwa/merged_replicate/bigwig/R4.mRp.clN.bigWig",
             "./Revision_R_codes/006_atac_seq/bw/R4.atac.cpm.bigWig")
normalise_bw("./Revision_R_codes/006_atac_seq/atac_nf/bwa/merged_replicate/bigwig/R5.mRp.clN.bigWig",
             "./Revision_R_codes/006_atac_seq/bw/R5.atac.cpm.bigWig")
normalise_bw("./Revision_R_codes/006_atac_seq/atac_nf/bwa/merged_replicate/bigwig/WT.mRp.clN.bigWig",
             "./Revision_R_codes/006_atac_seq/bw/WT.atac.cpm.bigWig")

# Prepare sequencing track overview

# Function transform a bigwig file into a dataframe with scores reported for a given step and a given interval

bigwig_to_table <- function(bigwig_file, chromosome, interval_start, interval_end, step) {
  # Check if the file exists
  if (!file.exists(bigwig_file)) {
    stop("The BigWig file does not exist.")
  }
  # Import the BigWig file for the specified chromosome and interval
  message("Reading BigWig file...")
  bw_data <- import(bigwig_file, which = GRanges(chromosome, IRanges(interval_start, interval_end)))
  # Convert BigWig data to a data.frame for easier handling
  bw_df <- as.data.frame(bw_data)
  # If no data exists in the interval, return an empty data frame
  if (nrow(bw_df) == 0) {
    message("No data found in the specified interval.")
    return(data.frame(Chromosome = character(), Start = integer(), End = integer(), Score = numeric()))
  }
  # Generate the intervals based on the provided step size
  starts <- seq(interval_start, interval_end, by = step)
  ends <- pmin(starts + step - 1, interval_end)
  # Initialize scores
  scores <- numeric(length(starts))
  # Iterate through each interval and calculate the mean score
  for (i in seq_along(starts)) {
    overlapping <- bw_df[bw_df$start <= ends[i] & bw_df$end >= starts[i], ]
    if (nrow(overlapping) > 0) {
      # Calculate the mean score for the overlapping region
      overlap_scores <- overlapping$score
      scores[i] <- mean(overlap_scores, na.rm = TRUE)
    } else {
      # Assign NA if no overlap
      scores[i] <- NA
    }
  }
  # Create the results table
  results <- data.frame(seqnames = chromosome, start = starts, end = ends, score = scores)
  return(results)
}

# Generate table over the entire locus with a step of 50bp
# OTX2 locus
# chr14:56,794,000-56,911,500

GOI.gr <- GRanges(seqnames = "chr14", ranges = IRanges(start = 56794000, end = 56911500))

# Load bigwig tracks and transform into df

AD.atac.bw <- "./Revision_R_codes/006_atac_seq/bw/AD.atac.cpm.bigWig"
  GOI.AD.atac.df <- bigwig_to_table(AD.atac.bw, "chr14", start(GOI.gr), end(GOI.gr), 50) %>% mutate(sample = "AD")

R4.atac.bw <- "./Revision_R_codes/006_atac_seq/bw/R4.atac.cpm.bigWig"
  GOI.R4.atac.df <- bigwig_to_table(R4.atac.bw, "chr14", start(GOI.gr), end(GOI.gr), 50) %>% mutate(sample = "R4")
  
R5.atac.bw <- "./Revision_R_codes/006_atac_seq/bw/R5.atac.cpm.bigWig"
  GOI.R5.atac.df <- bigwig_to_table(R5.atac.bw, "chr14", start(GOI.gr), end(GOI.gr), 50) %>% mutate(sample = "R5")
  
WT.atac.bw <- "./Revision_R_codes/006_atac_seq/bw/WT.atac.cpm.bigWig"
  GOI.WT.atac.df <- bigwig_to_table(WT.atac.bw, "chr14", start(GOI.gr), end(GOI.gr), 50) %>% mutate(sample = "WT")
  
atac.ov.df <- rbind(GOI.AD.atac.df, GOI.R4.atac.df, GOI.R5.atac.df, GOI.WT.atac.df) %>%
  mutate(sample = factor(sample, levels = c("WT", "AD", "R4", "R5")))

# Define markers to delineate the locus and the OTX2 gene

# OTX2
# chr14:56799905-56810479 (-)
# extended by 2kb to include upstream and dowsteam regulatory elements

markers.df <- data.frame(
  name = c("gene_down", "gene_end", "gene_start", "gene_up", "R4_start", "R45_bound", "R5_end"),
  position = c(56799905-2500, 56799905, 56810479, 56810479+2500, 56884635, 56898377, 56907609+1000)/1e6
)

# Define a colour palette
sample_cols <- c(
  "WT"  = "#A7A8A9",
  "AD"  = "#55AC56",
  "R4"  = "#6B8EC9",
  "R5"  = "#A72A2B")

# Plot the ATAC-seq signal across the locus for each sample
   
atac.ov.plot <- atac.ov.df %>%
  ggplot(aes(x = start / 1e6, y = score, fill = sample)) +
  geom_area(col = NA) +
  geom_line(linewidth = 0.25) +
  scale_fill_manual(values = sample_cols) +
  scale_colour_manual(values = "black") +
  xlab("Coordinate on chr14 [hg38]") +
  ylab("ATAC-seq signal [CPM]") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = markers.df$position, linetype = "dashed") +
  facet_wrap(~sample, ncol = 1) +
  theme_bw() +
  theme(aspect.ratio = 0.2,
    panel.grid = element_blank(), legend.position = "none")

pdf("./Revision_R_codes/006_atac_seq/rplots/atac.ov.pdf", width=15, height=8, useDingbats=FALSE)
atac.ov.plot
dev.off()

# Define consensus atac peaks across all samples
# We use the peaks identified by MACS2

macs2.peaks.df <- read.csv("./Revision_R_codes/006_atac_seq/atac_nf/bwa/merged_replicate/macs2/broad_peak/consensus/consensus_peaks.mRp.clN.bed", header = FALSE, sep = "\t") %>%
  setNames(c("chr", "start", "end", "name", "score", "strand")) %>% filter(chr == "chr14" & start >= start(GOI.gr) & end <= end(GOI.gr))

all.peaks.df <- macs2.peaks.df %>% arrange(start) %>% dplyr::select(chr, start, end)
# Export bed file
write.table(all.peaks.df, "./Revision_R_codes/006_atac_seq/bed/otx2.atac.peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#==================================================#
# Signal quantification
#==================================================#

# We quantify signal over three domains: OTX2 gene bodies, R123, R4 and R5 domains

gene.body.gr <- GRanges(seqnames = "chr14", ranges = IRanges(start = 56799905-2500, end = 56810479+2500), strand = "-")
  gene.body.gr$domain <- "OTX2_gene"

# Load atac peak coordinates
atac.peak.gr <- read.csv("./Revision_R_codes/006_atac_seq/bed/otx2.atac.peaks.bed", header = FALSE, sep = "\t") %>%
  setNames(c("seqnames", "start", "end")) %>% GRanges()

R123.peak.gr <- atac.peak.gr[start(atac.peak.gr) >= 56840404 & end(atac.peak.gr) <= 56884635]
  R123.peak.gr$domain <- "R123_domain"

R4.peak.gr <- atac.peak.gr[start(atac.peak.gr) >= 56884635 & end(atac.peak.gr) <= 56898377]
  R4.peak.gr$domain <- "R4_domain"

R5.peak.gr <- atac.peak.gr[start(atac.peak.gr) >= 56898377 & end(atac.peak.gr) <= 56907609+1000]
  R5.peak.gr$domain <- "R5_domain"

domain.gr <- c(gene.body.gr, R123.peak.gr, R4.peak.gr, R5.peak.gr)

# Function to quantify signal over a domain

quantify_signal_func <- function(atac.bw.path, domain_gr, sample) {
  atac.bw <- import(atac.bw.path, which = domain.gr)
  atac.cpm.df <- mergeByOverlaps(atac.bw, domain.gr) %>% as.data.frame() %>%
    dplyr::select(domain, score) %>% mutate(sample = sample)
}

# Compute

WT.atac.df <- quantify_signal_func(WT.atac.bw, domain.gr, "WT")
AD.atac.df <- quantify_signal_func(AD.atac.bw, domain.gr, "AD")
R4.atac.df <- quantify_signal_func(R4.atac.bw, domain.gr, "R4")
R5.atac.df <- quantify_signal_func(R5.atac.bw, domain.gr, "R5")

# Aggregate data
atac.quant.df <- rbind(WT.atac.df, AD.atac.df, R4.atac.df, R5.atac.df) %>%
  mutate(sample = factor(sample, levels = c("WT", "AD", "R4", "R5")))

# Correct values for R4 and R5 domains
atac.quant.df <- atac.quant.df %>%
  mutate(score = case_when(
    sample == "R4" & domain == "R4_domain" ~ 0,
    sample == "R5" & domain == "R5_domain" ~ 0,
    TRUE ~ score
  ))

# Plot

atac.quant.plot <- atac.quant.df %>%
  ggplot(aes(x = sample, y = score, fill = sample)) +
  geom_boxplot(
    outliers = FALSE,
    width         = 0.6,
    linewidth     = 0.4,
    alpha         = 0.85
  ) +
  scale_fill_manual(values   = sample_cols) +
  scale_colour_manual(values = sample_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~ domain, scales = "free_y",
             nrow = 1,
             labeller = as_labeller(c(
               "OTX2_gene" = "OTX2 gene body",
               "R123_domain" = "R123 domain",
               "R4_domain" = "R4 domain",
               "R5_domain" = "R5 domain"
             ))) +
  xlab(NULL) +
  ylab("ATAC-seq signal [CPM]") +
  theme_bw() +
  theme(
    aspect.ratio     = 1.5,
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 9, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 8),
    axis.title.y     = element_text(size = 9),
    legend.position  = "none"
  )

# Compute stats

stat.test <- atac.quant.df %>%
  group_by(domain) %>%
  wilcox_test(score ~ sample, ref.group = "WT") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(
    p.col       = "p.adj",
    cutpoints   = c(0, 0.001, 0.01, 0.05, 1),
    symbols     = c("***", "**", "*", "ns")
  ) %>%
  add_xy_position(x = "sample", dodge = 0.8)

# Compute fold change of mean signal relative to WT
fc <- atac.quant.df %>%
  group_by(domain, sample) %>%
  summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(domain) %>%
  mutate(FC_vs_WT = median_score / median_score[sample == "WT"],
         log2FC_vs_WT = log2(median_score / median_score[sample == "WT"]))
fc

# Export plot

pdf("./Revision_R_codes/006_atac_seq/rplots/atac.quant.pdf", width=6, height=4, useDingbats=FALSE)
atac.quant.plot
dev.off()








