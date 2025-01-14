
# All variants

# all.SV.reads.txt
# deltaR1234.reads.txt
# deltaR123.reads.txt
# deltaR12.reads.txt
# deltaR1.reads.txt
# deltaR34.reads.txt
# deltaR3.reads.txt
# deltaR45.reads.txt
# deltaR5.reads.txt
# 
# deltaR12+iR34.reads.txt
# deltaR12+iR3+deltaR4.reads.txt
# 
# iR12345.reads.txt
# iR1234.reads.txt
# iR123.reads.txt
# iR12.reads.txt
# iR34.reads.txt
# iR3.reads.txt
# iR45.reads.txt
# iR5.reads.txt

# deltaR12+iR34 and deltaR12+iR3+deltaR4.reads are supproting the same genotype which is deltaR12_iR34
cat deltaR12+iR34.reads.txt deltaR12+iR3+deltaR4.reads.txt > deltaR12_iR34.reads.txt

################################
### Define routine to compute coverages for variants with inversion
# Read starting upstream or downstream are processed independently (strand polarity of downstream reads is inverted)

cd /lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen

# Set the sample name as a variable
SAMPLE_NAME="iR12345"
SAMPLE_NAME="iR1234"
SAMPLE_NAME="iR123"
SAMPLE_NAME="iR12"
SAMPLE_NAME="iR34"
SAMPLE_NAME="iR3"
SAMPLE_NAME="iR45"
SAMPLE_NAME="iR5"
SAMPLE_NAME="deltaR12_iR34"

# Filter reads
samtools view -@ 12 -m 2G -bS -N ${SAMPLE_NAME}.reads.txt -o ${SAMPLE_NAME}.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_loxP6_all_reads.bam
samtools index -@ 12 ${SAMPLE_NAME}.bam

# Process reads starting upstream
samtools view -@ 12 ${SAMPLE_NAME}.bam chr14:1-56840000 | cut -f1 > ${SAMPLE_NAME}.up.read_names.txt
# Step 1: Extract primary reads from the BAM file using read names
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep '^@' > ${SAMPLE_NAME}.header
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep -Ff ${SAMPLE_NAME}.up.read_names.txt | grep -v -P '\t(2048|2049)\t' > ${SAMPLE_NAME}.up.primary_reads.sam
cat ${SAMPLE_NAME}.header ${SAMPLE_NAME}.up.primary_reads.sam > ${SAMPLE_NAME}.up.primary_reads.header.sam
# Step 2: Extract supplementary alignments for those reads
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep -Ff ${SAMPLE_NAME}.up.read_names.txt | grep -P '\t(2048|2049)\t' > ${SAMPLE_NAME}.up.supplementary_reads.sam
cat ${SAMPLE_NAME}.header ${SAMPLE_NAME}.up.supplementary_reads.sam > ${SAMPLE_NAME}.up.supplementary_reads.header.sam
# Step 3: Combine primary and supplementary reads into one BAM file
samtools view -@ 12 -Sb ${SAMPLE_NAME}.up.primary_reads.header.sam > ${SAMPLE_NAME}.up.primary_reads.bam
samtools view -@ 12 -Sb ${SAMPLE_NAME}.up.supplementary_reads.header.sam > ${SAMPLE_NAME}.up.supplementary_reads.bam
# Step 4: merge and index
samtools merge -@ 12 -f -o ${SAMPLE_NAME}.up.bam ${SAMPLE_NAME}.up.primary_reads.bam ${SAMPLE_NAME}.up.supplementary_reads.bam
samtools index -@ 12 ${SAMPLE_NAME}.up.bam
# Step 5: compute coverage
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ${SAMPLE_NAME}.up.bam -o ${SAMPLE_NAME}.up.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ${SAMPLE_NAME}.up.bam -o ${SAMPLE_NAME}.up.rev.bw
# Clean
rm -f ${SAMPLE_NAME}.up.read_names.txt
rm -f ${SAMPLE_NAME}.header
rm -f ${SAMPLE_NAME}.up.primary_reads.sam
rm -f ${SAMPLE_NAME}.up.primary_reads.header.sam
rm -f ${SAMPLE_NAME}.up.supplementary_reads.sam
rm -f ${SAMPLE_NAME}.up.supplementary_reads.header.sam
rm -f ${SAMPLE_NAME}.up.primary_reads.bam
rm -f ${SAMPLE_NAME}.up.supplementary_reads.bam
rm -f ${SAMPLE_NAME}.up.supplementary_reads.bam
rm -f ${SAMPLE_NAME}.down.bam

# Process reads starting downstream
samtools view -@ 12 ${SAMPLE_NAME}.bam chr14:56908000-99999999 | cut -f1 > ${SAMPLE_NAME}.down.read_names.txt
# Step 1: Extract primary reads from the BAM file using read names
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep '^@' > ${SAMPLE_NAME}.header
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep -Ff ${SAMPLE_NAME}.down.read_names.txt | grep -v -P '\t(2048|2049)\t' > ${SAMPLE_NAME}.down.primary_reads.sam
cat ${SAMPLE_NAME}.header ${SAMPLE_NAME}.down.primary_reads.sam > ${SAMPLE_NAME}.down.primary_reads.header.sam
# Step 2: Extract supplementary alignments for those reads
samtools view -@ 12 -h ${SAMPLE_NAME}.bam | grep -Ff ${SAMPLE_NAME}.down.read_names.txt | grep -P '\t(2048|2049)\t' > ${SAMPLE_NAME}.down.supplementary_reads.sam
cat ${SAMPLE_NAME}.header ${SAMPLE_NAME}.down.supplementary_reads.sam > ${SAMPLE_NAME}.down.supplementary_reads.header.sam
# Step 3: Combine primary and supplementary reads into one BAM file
samtools view -@ 12 -Sb ${SAMPLE_NAME}.down.primary_reads.header.sam > ${SAMPLE_NAME}.down.primary_reads.bam
samtools view -@ 12 -Sb ${SAMPLE_NAME}.down.supplementary_reads.header.sam > ${SAMPLE_NAME}.down.supplementary_reads.bam
# Step 4: merge and index
samtools merge -@ 12 -f -o ${SAMPLE_NAME}.down.bam ${SAMPLE_NAME}.down.primary_reads.bam ${SAMPLE_NAME}.down.supplementary_reads.bam
samtools index -@ 12 ${SAMPLE_NAME}.down.bam
# Step 5: compute coverage
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ${SAMPLE_NAME}.down.bam -o ${SAMPLE_NAME}.down.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ${SAMPLE_NAME}.down.bam -o ${SAMPLE_NAME}.down.rev.bw
# Clean
rm -f ${SAMPLE_NAME}.down.read_names.txt
rm -f ${SAMPLE_NAME}.header
rm -f ${SAMPLE_NAME}.down.primary_reads.sam
rm -f ${SAMPLE_NAME}.down.primary_reads.header.sam
rm -f ${SAMPLE_NAME}.down.supplementary_reads.sam
rm -f ${SAMPLE_NAME}.down.supplementary_reads.header.sam
rm -f ${SAMPLE_NAME}.down.primary_reads.bam
rm -f ${SAMPLE_NAME}.down.supplementary_reads.bam
rm -f ${SAMPLE_NAME}.down.supplementary_reads.bam
rm -f ${SAMPLE_NAME}.down.bam

################################
# Manual validation

# deltaR12_iR34 is correct, but the previous approach does not allow representating it as all reads cover the entire domain

# We use the approach develop for the LCR variant based on strand polarity +-+ vs -+-

# Define the sample name as a variable
SAMPLE_NAME="deltaR12_iR34"

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
samtools view -@ 12 -m 2G -bS -N "reads_${pattern}.txt" -o "${SAMPLE_NAME}.${pattern}.bam" /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_loxP6_all_reads.bam
samtools index -@ 12 "${SAMPLE_NAME}.${pattern}.bam"
done

# Step 3: Generate coverage files with bamCoverage for each strand direction and pattern
for strand in forward reverse; do
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand ${strand} -b "${SAMPLE_NAME}.+-+.bam" -o "${SAMPLE_NAME}.${strand}.+-+.bw"
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand ${strand} -b "${SAMPLE_NAME}.-+-.bam" -o "${SAMPLE_NAME}.${strand}.-+-.bw"
done

################################
#### R processing

source("./Notebook/Scripts/bigwig_to_table.R")

chromosome <- "chr14"         # Replace with the desired chromosome
interval_start <- 56837000     # Replace with the start of the interval
interval_end <- 56912000       # Replace with the end of the interval
step <- 100

# iR123

# Orientation 1
# fwd
iR123.fwd.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR123.down.fwd.bw"
iR123.fwd.down.df <- bigwig_to_table(iR123.fwd.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR123.rev.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR123.down.rev.bw"
iR123.rev.down.df <- bigwig_to_table(iR123.rev.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR123.down.df <- iR123.fwd.down.df %>% left_join(iR123.rev.down.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR123.fwd.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR123.up.fwd.bw"
iR123.fwd.up.df <- bigwig_to_table(iR123.fwd.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)# rev
iR123.rev.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR123.up.rev.bw"
iR123.rev.up.df <- bigwig_to_table(iR123.rev.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR123.up.df <- iR123.fwd.up.df %>% left_join(iR123.rev.up.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR123.df <- iR123.down.df %>% left_join(iR123.up.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR123")
plot(iR123.df$start, iR123.df$score)

# iR1234

# Orientation 1
# fwd
iR1234.fwd.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR1234.down.fwd.bw"
iR1234.fwd.down.df <- bigwig_to_table(iR1234.fwd.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR1234.rev.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR1234.down.rev.bw"
iR1234.rev.down.df <- bigwig_to_table(iR1234.rev.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR1234.down.df <- iR1234.fwd.down.df %>% left_join(iR1234.rev.down.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR1234.fwd.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR1234.up.fwd.bw"
iR1234.fwd.up.df <- bigwig_to_table(iR1234.fwd.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)# rev
iR1234.rev.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR1234.up.rev.bw"
iR1234.rev.up.df <- bigwig_to_table(iR1234.rev.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR1234.up.df <- iR1234.fwd.up.df %>% left_join(iR1234.rev.up.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR1234.df <- iR1234.down.df %>% left_join(iR1234.up.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR1234")
plot(iR1234.df$start, iR1234.df$score)

# iR34

# Orientation 1
# fwd
iR34.fwd.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR34.down.fwd.bw"
iR34.fwd.down.df <- bigwig_to_table(iR34.fwd.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR34.rev.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/iR34.down.rev.bw"
iR34.rev.down.df <- bigwig_to_table(iR34.rev.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR34.down.df <- iR34.fwd.down.df %>% left_join(iR34.rev.down.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.1 = score)
# No information for orientation 2
# combine
iR34.df <- iR34.down.df %>% dplyr::select(start, score = score.1) %>% mutate(sample = "iR34")
plot(iR34.df$start, iR34.df$score)

# deltaR12_iR34

# Orientation 1
# fwd
deltaR12_iR34.fwd.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/deltaR12_iR34.forward.+-+.bw"
deltaR12_iR34.fwd.plus.minus.plus.df <- bigwig_to_table(deltaR12_iR34.fwd.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
deltaR12_iR34.rev.plus.minus.plus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/deltaR12_iR34.reverse.+-+.bw"
deltaR12_iR34.rev.plus.minus.plus.df <- bigwig_to_table(deltaR12_iR34.rev.plus.minus.plus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
deltaR12_iR34.plus.minus.plus.df <- deltaR12_iR34.fwd.plus.minus.plus.df %>% left_join(deltaR12_iR34.rev.plus.minus.plus.df) %>% mutate(score = score.rev-score.fwd) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
deltaR12_iR34.fwd.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/deltaR12_iR34.forward.-+-.bw"
deltaR12_iR34.fwd.minus.plus.minus.df <- bigwig_to_table(deltaR12_iR34.fwd.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
deltaR12_iR34.rev.minus.plus.minus.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_clonal_screen/deltaR12_iR34.reverse.-+-.bw"
deltaR12_iR34.rev.minus.plus.minus.df <- bigwig_to_table(deltaR12_iR34.rev.minus.plus.minus.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
deltaR12_iR34.minus.plus.minus.df <- deltaR12_iR34.fwd.minus.plus.minus.df %>% left_join(deltaR12_iR34.rev.minus.plus.minus.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
deltaR12_iR34.df <- deltaR12_iR34.plus.minus.plus.df %>% left_join(deltaR12_iR34.minus.plus.minus.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "deltaR12_iR34")
plot(deltaR12_iR34.df$start, deltaR12_iR34.df$score)

# Prepare summary table

loxP6.clonal.inversion.cas9.coverage.df <- rbind(iR123.df, iR1234.df, iR34.df, deltaR12_iR34.df)
saveRDS(loxP6.clonal.inversion.cas9.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/loxP6.clonal.inversion.cas9.coverage.df.rds")

# Plot

loxP6.clonal.inversion.cas9.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/loxP6.clonal.inversion.cas9.coverage.df.rds")

loxP6.clonal.inversion.cas9.palette <- wes_palette("Zissou1", length(unique(loxP6.clonal.inversion.cas9.coverage.df$sample)), type = "continuous")
loxP6.clonal.inversion.cas9.coverage.plot <- loxP6.clonal.inversion.cas9.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = sample)) +
  geom_area() +
  geom_line(col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(loxP6.clonal.inversion.cas9.palette)) +
  xlim(56837000/1e6, 56912000/1e6) +
  geom_vline(xintercept=56840404/1e6, linetype="dashed") +
  geom_vline(xintercept=56851954/1e6, linetype="dashed") +
  geom_vline(xintercept=56864719/1e6, linetype="dashed") +
  geom_vline(xintercept=56884635/1e6, linetype="dashed") +
  geom_vline(xintercept=56898377/1e6, linetype="dashed") +
  geom_vline(xintercept=56907609/1e6, linetype="dashed") +
  facet_wrap(~ sample, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), strip.background = element_blank(), legend.position="none")
loxP6.clonal.inversion.cas9.coverage.plot



