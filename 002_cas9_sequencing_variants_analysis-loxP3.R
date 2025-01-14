
module load HGI/softpack/groups/escramble/eSCRAMBLE/7

cd /lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen

# Define the sample name as a variable
SAMPLE_NAME="iR4_deltaR5"
SAMPLE_NAME="deltaR4_iR5"
SAMPLE_NAME="iR4"
SAMPLE_NAME="iR5"
SAMPLE_NAME="iR45"

# Step 1: Extract read IDs with +-+ and -+- patterns and output to separate files
samtools view -h "${SAMPLE_NAME}.bam" | \
awk '
    $1 !~ /^@/ {
        read_id = $1
        strand_pattern = (and($2, 16) ? "-" : "+")  # Primary alignment strand

        # Look for SA:Z tag for supplementary alignments
        for (i=12; i<=NF; i++) {
            if ($i ~ /^SA:Z:/) {
                n = split(substr($i,6), sa_fields, ";")
                for (j=1; j<=n; j++) {
                    if (sa_fields[j] != "") {
                        split(sa_fields[j], sa_segment, ",")
                        strand_pattern = strand_pattern sa_segment[3]
                    }
                }
            }
        }
        
        # Output read ID based on strand pattern
        if (strand_pattern == "+-+") print read_id > "reads_+-+.txt"
        else if (strand_pattern == "-+-") print read_id > "reads_-+-.txt"
    }
'

# Step 2: Filter original BAM file based on strand pattern read IDs, output to new BAM files, and index them
for pattern in "+-+" "-+-"; do
samtools view -@ 12 -m 2G -bS -N "reads_${pattern}.txt" -o "${SAMPLE_NAME}.${pattern}.bam" /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_LCR_all_reads.bam
samtools index -@ 12 "${SAMPLE_NAME}.${pattern}.bam"
done

# Step 3: Generate coverage files with bamCoverage for each strand direction and pattern
for strand in forward reverse; do
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand ${strand} -b "${SAMPLE_NAME}.+-+.bam" -o "${SAMPLE_NAME}.${strand}.+-+.bw"
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand ${strand} -b "${SAMPLE_NAME}.-+-.bam" -o "${SAMPLE_NAME}.${strand}.-+-.bw"
done

#### R processing

source("./Notebook/Scripts/bigwig_to_table.R")

chromosome <- "chr14"         # Replace with the desired chromosome
interval_start <- 56884000     # Replace with the start of the interval
interval_end <- 56907900       # Replace with the end of the interval
step <- 100

########
# deltaR4_iR5

# Orientation 1
# fwd
deltaR4_iR5.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/deltaR4_iR5.forward.+-+.bw"
deltaR4_iR5.fwd.plus.minus.plus.df <- bigwig_to_table(deltaR4_iR5.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
deltaR4_iR5.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/deltaR4_iR5.reverse.+-+.bw"
deltaR4_iR5.rev.plus.minus.plus.df <- bigwig_to_table(deltaR4_iR5.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
deltaR4_iR5.plus.minus.plus.df <- deltaR4_iR5.fwd.plus.minus.plus.df %>% left_join(deltaR4_iR5.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
deltaR4_iR5.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/deltaR4_iR5.forward.-+-.bw"
deltaR4_iR5.fwd.minus.plus.minus.df <- bigwig_to_table(deltaR4_iR5.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
deltaR4_iR5.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/deltaR4_iR5.reverse.-+-.bw"
deltaR4_iR5.rev.minus.plus.minus.df <- bigwig_to_table(deltaR4_iR5.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
deltaR4_iR5.minus.plus.minus.df <- deltaR4_iR5.fwd.minus.plus.minus.df %>% left_join(deltaR4_iR5.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
deltaR4_iR5.df <- deltaR4_iR5.plus.minus.plus.df %>% left_join(deltaR4_iR5.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "deltaR4_iR5")
plot(deltaR4_iR5.df$start, deltaR4_iR5.df$score)

#########
# iR4_deltaR5

# Orientation 1
# fwd
iR4_deltaR5.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4_deltaR5.forward.+-+.bw"
iR4_deltaR5.fwd.plus.minus.plus.df <- bigwig_to_table(iR4_deltaR5.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR4_deltaR5.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4_deltaR5.reverse.+-+.bw"
iR4_deltaR5.rev.plus.minus.plus.df <- bigwig_to_table(iR4_deltaR5.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR4_deltaR5.plus.minus.plus.df <- iR4_deltaR5.fwd.plus.minus.plus.df %>% left_join(iR4_deltaR5.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR4_deltaR5.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4_deltaR5.forward.-+-.bw"
iR4_deltaR5.fwd.minus.plus.minus.df <- bigwig_to_table(iR4_deltaR5.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR4_deltaR5.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4_deltaR5.reverse.-+-.bw"
iR4_deltaR5.rev.minus.plus.minus.df <- bigwig_to_table(iR4_deltaR5.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR4_deltaR5.minus.plus.minus.df <- iR4_deltaR5.fwd.minus.plus.minus.df %>% left_join(iR4_deltaR5.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR4_deltaR5.df <- iR4_deltaR5.plus.minus.plus.df %>% left_join(iR4_deltaR5.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR4_deltaR5")
plot(iR4_deltaR5.df$start, iR4_deltaR5.df$score)

#########
# iR4

# Orientation 1
# fwd
iR4.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4.forward.+-+.bw"
iR4.fwd.plus.minus.plus.df <- bigwig_to_table(iR4.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR4.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4.reverse.+-+.bw"
iR4.rev.plus.minus.plus.df <- bigwig_to_table(iR4.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR4.plus.minus.plus.df <- iR4.fwd.plus.minus.plus.df %>% left_join(iR4.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR4.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4.forward.-+-.bw"
iR4.fwd.minus.plus.minus.df <- bigwig_to_table(iR4.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR4.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR4.reverse.-+-.bw"
iR4.rev.minus.plus.minus.df <- bigwig_to_table(iR4.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR4.minus.plus.minus.df <- iR4.fwd.minus.plus.minus.df %>% left_join(iR4.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR4.df <- iR4.plus.minus.plus.df %>% left_join(iR4.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR4")
plot(iR4.df$start, iR4.df$score)

########
# iR5

# Orientation 1
# fwd
iR5.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR5.forward.+-+.bw"
iR5.fwd.plus.minus.plus.df <- bigwig_to_table(iR5.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR5.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR5.reverse.+-+.bw"
iR5.rev.plus.minus.plus.df <- bigwig_to_table(iR5.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR5.plus.minus.plus.df <- iR5.fwd.plus.minus.plus.df %>% left_join(iR5.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR5.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR5.forward.-+-.bw"
iR5.fwd.minus.plus.minus.df <- bigwig_to_table(iR5.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR5.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR5.reverse.-+-.bw"
iR5.rev.minus.plus.minus.df <- bigwig_to_table(iR5.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR5.minus.plus.minus.df <- iR5.fwd.minus.plus.minus.df %>% left_join(iR5.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR5.df <- iR5.plus.minus.plus.df %>% left_join(iR5.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR5")
plot(iR5.df$start, iR5.df$score)

########
# iR45

# Orientation 1
# fwd
iR45.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR45.forward.+-+.bw"
iR45.fwd.plus.minus.plus.df <- bigwig_to_table(iR45.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR45.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR45.reverse.+-+.bw"
iR45.rev.plus.minus.plus.df <- bigwig_to_table(iR45.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR45.plus.minus.plus.df <- iR45.fwd.plus.minus.plus.df %>% left_join(iR45.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR45.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR45.forward.-+-.bw"
iR45.fwd.minus.plus.minus.df <- bigwig_to_table(iR45.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR45.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/LCR_screen/iR45.reverse.-+-.bw"
iR45.rev.minus.plus.minus.df <- bigwig_to_table(iR45.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR45.minus.plus.minus.df <- iR45.fwd.minus.plus.minus.df %>% left_join(iR45.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR45.df <- iR45.plus.minus.plus.df %>% left_join(iR45.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR45")
plot(iR45.df$start, iR45.df$score)

# Prepare summary table

loxP3.inversion.cas9.coverage.df <- rbind(deltaR4_iR5.df, iR4_deltaR5.df, iR4.df, iR5.df, iR45.df)
saveRDS(loxP3.inversion.cas9.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/loxP3.inversion.cas9.coverage.df.rds")

# Plot

loxP3.inversion.cas9.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/loxP3.inversion.cas9.coverage.df.rds")

loxP3.inversion.cas9.palette <- wes_palette("Zissou1", length(unique(loxP3.inversion.cas9.coverage.df$sample)), type = "continuous")
loxP3.inversion.cas9.coverage.plot <- loxP3.inversion.cas9.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = sample)) +
  geom_area() +
  geom_line(col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(loxP3.inversion.cas9.palette)) +
  xlim(56884500/1e6, 56907900/1e6) +
  geom_vline(xintercept=56884634/1e6, linetype="dashed") +
  geom_vline(xintercept=56898376/1e6, linetype="dashed") +
  geom_vline(xintercept=56907608/1e6, linetype="dashed") +
  facet_wrap(~ sample, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), strip.background = element_blank(), legend.position="none")
loxP3.inversion.cas9.coverage.plot







