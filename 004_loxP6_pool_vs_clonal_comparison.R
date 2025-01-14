
#########################################################################
# Read coverage

# for the loxPsym 6 pool transfection data

# Select read names

loxP6.pool.read.names.df <- read.table("/lustre/scratch126/gengen/projects/escramble/Data/JK_data/read_names/read_names_trans.tsv", header = T, sep = "\t")
unique(loxP6.pool.read.names.df$name) # 32 architectures
# Mutate names for file manipulation
loxP6.pool.read.names.df <- loxP6.pool.read.names.df %>% mutate(file.name = gsub("∆", "delta", name))

# Compute the number of reads associated with each architecture
loxP6.pool.read.count.df <- loxP6.pool.read.names.df %>% group_by(name) %>% summarise(read.count = n())

# Prepare read name txt files

for (i in 1:length(unique(loxP6.pool.read.count.df$name))) {
  genotype.i <- unique(loxP6.pool.read.count.df$name)[i]
  print(genotype.i)
  read.names.i <- (loxP6.pool.read.names.df %>% filter(name == genotype.i))$fiber
  name.i <- unique((loxP6.pool.read.names.df %>% filter(name == genotype.i))$file.name)
  path.i <- paste0("/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/", name.i, ".reads.txt")
  write.table(read.names.i, path.i, sep="\t", col.names = F, row.names = F, quote = F)
}

# WT architectures reads are defined by the group of reads that do not support a structural variant
# Save all reads supporting a structural variant
SV.read.names <- unique(loxP6.pool.read.names.df$fiber)
path <- paste0("/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/", "all.SV", ".reads.txt")
write.table(SV.read.names, path, sep="\t", col.names = F, row.names = F, quote = F)

# List of detected variants

# all.SV
# deltaR1
# deltaR12
# deltaR123
# deltaR1234
# deltaR12345
# deltaR2
# deltaR23
# deltaR234
# deltaR2345
# deltaR3
# deltaR34
# deltaR345
# deltaR4
# deltaR45
# deltaR5
# iR1
# iR12
# iR123
# iR1234
# iR12345
# iR2
# iR23
# iR234
# iR2345
# iR3
# iR34
# iR345
# iR4
# iR5
# deltaR12+iR34
# deltaR12+iR3+deltaR4
# deltaR12+iR3+deltaR45

# Combine bam file from all sorting gates

cd /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2

# loxP6

samtools merge -@ 12 \
-o OTX2_trans_all_reads.bam \
OTX2_trans_bright_R1_OTX2.bam OTX2_trans_bright_R2_OTX2.bam OTX2_trans_dark_R1_OTX2.bam OTX2_trans_dark_R2_OTX2.bam \
OTX2_trans_dim_R1_OTX2.bam OTX2_trans_dim_R2_OTX2.bam OTX2_trans_medium_R1_OTX2.bam OTX2_trans_medium_R2_OTX2.bam
samtools index -@ 12 OTX2_trans_all_reads.bam

# Compute coverages

cd /lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/
  module load HGI/softpack/groups/escramble/eSCRAMBLE/8

# Simple deletion

samtools view /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam | grep -v -F -f all.SV.reads.txt | cut -f 1 > WT.reads.txt

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -h /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam | grep -v -F -f all.SV.reads.txt | samtools view -b -o WT.bam
samtools index -@ 12 WT.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b WT.bam -o WT.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR1.reads.txt -o deltaR1.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR1.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR1.bam -o deltaR1.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR12.reads.txt -o deltaR12.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR12.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR12.bam -o deltaR12.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR123.reads.txt -o deltaR123.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR123.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR123.bam -o deltaR123.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR1234.reads.txt -o deltaR1234.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR1234.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR1234.bam -o deltaR1234.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR12345.reads.txt -o deltaR12345.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR12345.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR12345.bam -o deltaR12345.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR2.reads.txt -o deltaR2.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR2.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR2.bam -o deltaR2.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR23.reads.txt -o deltaR23.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR23.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR23.bam -o deltaR23.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR234.reads.txt -o deltaR234.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR234.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR234.bam -o deltaR234.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR2345.reads.txt -o deltaR2345.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR2345.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR2345.bam -o deltaR2345.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR3.reads.txt -o deltaR3.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR3.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR3.bam -o deltaR3.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR34.reads.txt -o deltaR34.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR34.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR34.bam -o deltaR34.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR345.reads.txt -o deltaR345.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR345.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR345.bam -o deltaR345.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR4.reads.txt -o deltaR4.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR4.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR4.bam -o deltaR4.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR45.reads.txt -o deltaR45.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR45.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR45.bam -o deltaR45.bw'

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o coverage.log -J coverage \
'samtools view -@ 12 -m 2G -bS -N deltaR5.reads.txt -o deltaR5.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
samtools index -@ 12 deltaR5.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b deltaR5.bam -o deltaR5.bw'

# All deletion variants are correct

# Plot deletion variants coverage plots

source("./Notebook/Scripts/bigwig_to_table.R")

chromosome <- "chr14"         
interval_start <- 56837000
interval_end <- 56912000   
step <- 100

WT.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/WT.bw"
WT.loxp6.pool.df <- bigwig_to_table(WT.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "WT")

deltaR1.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR1.bw"
deltaR1.loxp6.pool.df <- bigwig_to_table(deltaR1.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR1")

deltaR12.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR12.bw"
deltaR12.loxp6.pool.df <- bigwig_to_table(deltaR12.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR12")

deltaR123.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR123.bw"
deltaR123.loxp6.pool.df <- bigwig_to_table(deltaR123.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR123")

deltaR1234.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR1234.bw"
deltaR1234.loxp6.pool.df <- bigwig_to_table(deltaR1234.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR1234")

deltaR12345.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR12345.bw"
deltaR12345.loxp6.pool.df <- bigwig_to_table(deltaR12345.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR12345")

deltaR2.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR2.bw"
deltaR2.loxp6.pool.df <- bigwig_to_table(deltaR2.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR2")

deltaR23.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR23.bw"
deltaR23.loxp6.pool.df <- bigwig_to_table(deltaR23.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR23")

deltaR234.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR234.bw"
deltaR234.loxp6.pool.df <- bigwig_to_table(deltaR234.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR234")

deltaR2345.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR2345.bw"
deltaR2345.loxp6.pool.df <- bigwig_to_table(deltaR2345.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR2345")

deltaR3.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR3.bw"
deltaR3.loxp6.pool.df <- bigwig_to_table(deltaR3.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR3")

deltaR34.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR34.bw"
deltaR34.loxp6.pool.df <- bigwig_to_table(deltaR34.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR34")

deltaR345.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR345.bw"
deltaR345.loxp6.pool.df <- bigwig_to_table(deltaR345.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR345")

deltaR4.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR4.bw"
deltaR4.loxp6.pool.df <- bigwig_to_table(deltaR4.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR4")

deltaR45.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR45.bw"
deltaR45.loxp6.pool.df <- bigwig_to_table(deltaR45.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR45")

deltaR5.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/deltaR5.bw"
deltaR5.loxp6.pool.df <- bigwig_to_table(deltaR5.path, chromosome, interval_start, interval_end, step) %>% mutate(variant = "deltaR5")

# Prepare summary table

loxP6.pool.deletion.cas9.coverage.df <- rbind(WT.loxp6.pool.df, deltaR1.loxp6.pool.df, deltaR12.loxp6.pool.df, deltaR123.loxp6.pool.df, deltaR1234.loxp6.pool.df, deltaR12345.loxp6.pool.df, deltaR2.loxp6.pool.df, deltaR23.loxp6.pool.df, deltaR234.loxp6.pool.df, deltaR2345.loxp6.pool.df, deltaR3.loxp6.pool.df, deltaR34.loxp6.pool.df, deltaR345.loxp6.pool.df, deltaR4.loxp6.pool.df, deltaR45.loxp6.pool.df, deltaR5.loxp6.pool.df)
saveRDS(loxP6.pool.deletion.cas9.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/005_OTX2_loxP6_pool_vs_clonal_comparison/rds/loxP6.pool.deletion.cas9.coverage.df.rds")

# Plot

loxP6.pool.deletion.cas9.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/005_OTX2_loxP6_pool_vs_clonal_comparison/rds/loxP6.pool.deletion.cas9.coverage.df.rds")

loxP6.pool.deletion.cas9.palette <- c(wes_palette("Zissou1", length(unique(loxP6.pool.deletion.cas9.coverage.df$variant))-1, type = "continuous"), "#BCBEC0")
loxP6.pool.deletion.cas9.coverage.plot <- loxP6.pool.deletion.cas9.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = variant)) +
  geom_area() +
  geom_line(col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(loxP6.pool.deletion.cas9.palette)) +
  xlim(56837000/1e6, 56912000/1e6) +
  geom_vline(xintercept=56840404/1e6, linetype="dashed") +
  geom_vline(xintercept=56851954/1e6, linetype="dashed") +
  geom_vline(xintercept=56864719/1e6, linetype="dashed") +
  geom_vline(xintercept=56884635/1e6, linetype="dashed") +
  geom_vline(xintercept=56898377/1e6, linetype="dashed") +
  geom_vline(xintercept=56907609/1e6, linetype="dashed") +
  facet_wrap(~ variant, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.minor = element_blank(), strip.background = element_blank(), legend.position="none")
loxP6.pool.deletion.cas9.coverage.plot

pdf("/lustre/scratch126/gengen/projects/escramble/Notebook/Rplots/OTX2/loxP6.pool.deletion.coverage.pdf", width=10, height=10, useDingbats=FALSE)
loxP6.pool.deletion.cas9.coverage.plot
dev.off()

# Simple inversion

##########################################################################################################################
### Define routine to compute coverages for variants with inversion
# Read starting upstream or downstream are processed independently (strand polarity of downstream reads is inverted)

cd /lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/
  module load HGI/softpack/groups/escramble/eSCRAMBLE/8

# Set the sample name as a variable
SAMPLE_NAME="iR1234"
SAMPLE_NAME="iR34"
SAMPLE_NAME="iR5"
SAMPLE_NAME="iR2"
SAMPLE_NAME="deltaR12+iR3+deltaR4"
SAMPLE_NAME="deltaR12+iR3+deltaR45"

# Filter reads
samtools view -@ 12 -m 2G -bS -N ${SAMPLE_NAME}.reads.txt -o ${SAMPLE_NAME}.bam /lustre/scratch126/gengen/projects/escramble/Data/JK_data/bam_2/OTX2_trans_all_reads.bam
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

# Only iR34 and iR1234 are correct or have enough reads to support them

# Plot deletion variants coverage plots

# iR34

# Orientation 1
# fwd
iR34.fwd.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR34.down.fwd.bw"
iR34.fwd.down.df <- bigwig_to_table(iR34.fwd.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR34.rev.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR34.down.rev.bw"
iR34.rev.down.df <- bigwig_to_table(iR34.rev.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR34.down.df <- iR34.fwd.down.df %>% left_join(iR34.rev.down.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR34.fwd.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR34.up.fwd.bw"
iR34.fwd.up.df <- bigwig_to_table(iR34.fwd.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)# rev
iR34.rev.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR34.up.rev.bw"
iR34.rev.up.df <- bigwig_to_table(iR34.rev.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR34.up.df <- iR34.fwd.up.df %>% left_join(iR34.rev.up.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR34.df <- iR34.down.df %>% left_join(iR34.up.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR34")
plot(iR34.df$start, iR34.df$score)

# iR1234

# Orientation 1
# fwd
iR1234.fwd.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR1234.down.fwd.bw"
iR1234.fwd.down.df <- bigwig_to_table(iR1234.fwd.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)
# rev
iR1234.rev.down.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR1234.down.rev.bw"
iR1234.rev.down.df <- bigwig_to_table(iR1234.rev.down.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR1234.down.df <- iR1234.fwd.down.df %>% left_join(iR1234.rev.down.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.1 = score)
# Orientation 2
# fwd
iR1234.fwd.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR1234.up.fwd.bw"
iR1234.fwd.up.df <- bigwig_to_table(iR1234.fwd.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.fwd = score)# rev
iR1234.rev.up.path <- "/lustre/scratch126/gengen/projects/escramble/Data/JK_data/coverage/loxP6_pool_screen/iR1234.up.rev.bw"
iR1234.rev.up.df <- bigwig_to_table(iR1234.rev.up.path, chromosome, interval_start, interval_end, step) %>% dplyr::rename(score.rev = score)
# combine
iR1234.up.df <- iR1234.fwd.up.df %>% left_join(iR1234.rev.up.df) %>% mutate(score = score.fwd-score.rev) %>% dplyr::select(start, score.2 = score)
# combine both strand orientation
iR1234.df <- iR1234.down.df %>% left_join(iR1234.up.df) %>% mutate(score=score.1+score.2) %>% dplyr::select(start, score) %>% mutate(sample = "iR1234")
plot(iR1234.df$start, iR1234.df$score)

# Prepare summary table

loxP6.pool.inversion.cas9.coverage.df <- rbind(iR34.df, iR1234.df)
saveRDS(loxP6.pool.inversion.cas9.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/005_OTX2_loxP6_pool_vs_clonal_comparison/rds/loxP6.pool.inversion.cas9.coverage.df.rds")

# Plot

loxP6.pool.inversion.cas9.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/005_OTX2_loxP6_pool_vs_clonal_comparison/rds/loxP6.pool.inversion.cas9.coverage.df.rds")

loxP6.pool.inversion.cas9.palette <- wes_palette("Zissou1", length(unique(loxP6.pool.inversion.cas9.coverage.df$sample)), type = "continuous")
loxP6.pool.inversion.cas9.coverage.plot <- loxP6.pool.inversion.cas9.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = sample)) +
  geom_area() +
  geom_line(col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(loxP6.pool.inversion.cas9.palette)) +
  xlim(56837000/1e6, 56912000/1e6) +
  geom_vline(xintercept=56840404/1e6, linetype="dashed") +
  geom_vline(xintercept=56851954/1e6, linetype="dashed") +
  geom_vline(xintercept=56864719/1e6, linetype="dashed") +
  geom_vline(xintercept=56884635/1e6, linetype="dashed") +
  geom_vline(xintercept=56898377/1e6, linetype="dashed") +
  geom_vline(xintercept=56907609/1e6, linetype="dashed") +
  facet_wrap(~ sample, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.minor = element_blank(), strip.background = element_blank(), legend.position="none")
loxP6.pool.inversion.cas9.coverage.plot

pdf("/lustre/scratch126/gengen/projects/escramble/Notebook/Rplots/OTX2/loxP6.pool.inversion.coverage.pdf", width=10, height=10, useDingbats=FALSE)
loxP6.pool.inversion.cas9.coverage.plot
dev.off()

#########################################################################

# Compare SV diversity for SV originating form the loxP6 cells - clonal vs pool transfection

# Load results

# clonal

SV.6.loxPsym.clone.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/SV.6.loxPsym.clone.wrap.df.rds")

# pool

SV.summary.df <- read.table("/lustre/scratch126/gengen/projects/escramble/Data/JK_data/summary_table/OTX2_scramble_outcomes.tsv", sep = "\t", header = T)
unique(SV.summary.df$experiment)
# 3_loxPsym_clone, 6_loxPsym_transfection, 6_loxPsym_clone 

SV.6.loxPsym.pool.df <- SV.summary.df %>% filter(experiment == "6_loxPsym_transfection")
# wrap data - Filling missing observations for plotting and score computation
genotype <- unique(SV.6.loxPsym.pool.df$variant_name)
wrap.df <-   tibble(variant_name = rep(genotype, 8)) %>% arrange(variant_name) %>% 
  mutate(gate = rep(c("dark", "dim", "medium", "bright"), 2*length(unique(SV.6.loxPsym.pool.df$variant_name)))) %>% 
  mutate(replicate = rep(c(rep(1,4), rep(2,4)), length(unique(SV.6.loxPsym.pool.df$variant_name))))
SV.6.loxPsym.pool.wrap.df <- SV.6.loxPsym.pool.df %>% full_join(wrap.df, by = c("variant_name", "gate", "replicate")) %>%
  mutate(perc_reads = if_else(is.na(perc_reads), 0, perc_reads))
saveRDS(SV.6.loxPsym.pool.wrap.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/SV.6.loxPsym.pool.wrap.df.rds")

SV.6.loxPsym.pool.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/SV.6.loxPsym.pool.wrap.df.rds")

SV.6.loxPsym.pool.clean.df <- SV.6.loxPsym.pool.df %>%
  dplyr::select(name = variant_name, genotype = variant_definition, replicate, gate, read_count = supp_reads, read_percentage = perc_reads) %>%
  mutate(read_count = replace_na(read_count, 0)) %>% arrange(name)
write.table(SV.6.loxPsym.pool.clean.df, "./Notebook/Summary/OTX2/loxPsym.6.pool.cas9.summary.tsv", sep="\t", col.names = T, row.names = F, quote = F)

#####
# Filter variants

# We exclude variants with less than observations in 2 gates and with a minimum of 25 total reads from the analysis
# As done previously

SV.6.loxPsym.clone.summary.df <- SV.6.loxPsym.clone.df %>% group_by(variant_name, replicate) %>% summarise(na_count = sum(is.na(supp_reads))) %>% filter(na_count <= 2)
SV.6.loxPsym.clone.summary.df.2 <- SV.6.loxPsym.clone.df %>% group_by(variant_name, replicate) %>% summarise(count = sum(supp_reads, na.rm = T)) %>% filter(count >= 25)
SV.6.loxPsym.clone.variant.filter <- intersect(unique(SV.6.loxPsym.clone.summary.df$variant_name), unique(SV.6.loxPsym.clone.summary.df.2$variant_name))

SV.6.loxPsym.pool.summary.df <- SV.6.loxPsym.pool.df %>% group_by(variant_name, replicate) %>% summarise(na_count = sum(is.na(supp_reads))) %>% filter(na_count <= 2)
SV.6.loxPsym.pool.summary.df.2 <- SV.6.loxPsym.pool.df %>% group_by(variant_name, replicate) %>% summarise(count = sum(supp_reads, na.rm = T)) %>% filter(count >= 25)
SV.6.loxPsym.pool.variant.filter <- intersect(unique(SV.6.loxPsym.pool.summary.df$variant_name), unique(SV.6.loxPsym.pool.summary.df.2$variant_name))

#####
# Combine and format data

SV.6.loxPsym.clone.format.df <- SV.6.loxPsym.clone.df %>% dplyr::select(variant_name, gate, replicate, supp_reads, perc_reads) %>% mutate(across(everything(), ~ replace_na(., 0))) %>% 
  group_by(variant_name, gate) %>% summarise(perc_reads = mean(perc_reads, na.rm = T)) %>% ungroup() %>% 
  group_by(gate) %>% mutate(perc_reads_scaled = perc_reads / sum(perc_reads) * 100) %>% ungroup() %>% 
  mutate(exp = "clone")

SV.6.loxPsym.pool.format.df <- SV.6.loxPsym.pool.df %>% dplyr::select(variant_name, gate, replicate, supp_reads, perc_reads) %>% mutate(across(everything(), ~ replace_na(., 0))) %>% 
  group_by(variant_name, gate) %>% summarise(perc_reads = mean(perc_reads, na.rm = T)) %>% ungroup() %>% 
  group_by(gate) %>% mutate(perc_reads_scaled = perc_reads / sum(perc_reads) * 100) %>% ungroup() %>% 
  mutate(exp = "pool")

SV.6.loxPsym.SV.diversity.df <- rbind(SV.6.loxPsym.clone.format.df, SV.6.loxPsym.pool.format.df)
saveRDS(SV.6.loxPsym.SV.diversity.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/SV.6.loxPsym.SV.diversity.df.rds")

#####
# Plot read distributions for the pool transfection experiment

SV.6.loxPsym.pool.filter.df <- SV.6.loxPsym.pool.df %>% filter(variant_name %in% SV.6.loxPsym.pool.variant.filter)

SV.6.loxPsym.pool.palette <- wes_palette("Zissou1", length(unique(SV.6.loxPsym.pool.filter.df$variant_name)), type = "continuous")
SV.6.loxPsym.pool.gate.plot <- SV.6.loxPsym.pool.filter.df %>% 
  mutate(SAMPLE = fct_relevel(gate, "dark", "dim", "medium", "bright")) %>%
  ggplot(aes(x=SAMPLE, y=perc_reads, group=1, col = variant_name)) +
  geom_point(aes(x = SAMPLE, y = perc_reads), size = 1) +
  geom_point(stat='summary', fun.y=sum, size = 2) +
  geom_line(stat='summary', fun.y=sum) +
  ylab("Percent of events in sorting gate") +
  scale_colour_manual(values=c(SV.6.loxPsym.pool.palette)) +
  facet_wrap(~ variant_name, scales = "free", ncol = 5) +
  theme_bw() + theme(aspect.ratio=1, strip.background = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none")
SV.6.loxPsym.pool.gate.plot

pdf("./Notebook/Rplots/OTX2/SV.6.loxPsym.pool.gate.plot.pdf", width=10, height=10, useDingbats=FALSE)
SV.6.loxPsym.pool.gate.plot
dev.off()

#####
# Compute gene expression score for the identified variants

SV.6.loxPsym.pool.score.df <- SV.6.loxPsym.pool.filter.df %>%
  group_by(variant_name, replicate) %>%
  mutate(scale = perc_reads/sum(perc_reads)) %>% 
  dplyr::select(variant_name, gate, replicate, scale) %>% ungroup() %>% 
  pivot_wider(names_from = gate, values_from = scale) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(score = (-1*dark)+(-0.5*dim)+(0.5*medium)+(1*bright))
# Scale scores to WT and full deletion levels
SV.6.loxPsym.pool.means.df <- SV.6.loxPsym.pool.score.df %>% group_by(variant_name) %>% summarise(mean=mean(score, na.rm = T))
# WT = 0.0160, ∆R12345 = -0.639
SV.6.loxPsym.pool.score.df <- SV.6.loxPsym.pool.score.df %>% mutate(score.norm = ((score+0.639)/(0.0160+0.639))*100)

# Plot

# Compute means and differences

SV.6.loxPsym.pool.score.stats.df <- SV.6.loxPsym.pool.score.df %>% group_by(variant_name, replicate) %>% summarise(mean=mean(score.norm, na.rm = T), diff=max(score.norm, na.rm = T)-mean(score.norm, na.rm = T))

SV.6.loxPsym.pool.score.plot <- SV.6.loxPsym.pool.score.df %>% 
  ggplot() + 
  geom_bar(aes(x=reorder(variant_name, score.norm), y=score.norm), position = "dodge", stat = "summary", fun.y = "mean", col = "black", size = 0.25, width = 0) +
  geom_errorbar(data=SV.6.loxPsym.pool.score.stats.df, aes(x=variant_name, ymin=mean-diff, ymax=mean+diff), linewidth = 0.5) +
  geom_point(aes(x=reorder(variant_name, -score.norm), y=score.norm), size = 1.5) +
  xlab("Gene expression score") +
  geom_hline(yintercept=0, linetype="dashed", linewidth = 0.25) +
  geom_hline(yintercept=100, linetype="dashed", linewidth = 0.25) +
  scale_y_continuous("Relative gene expression", breaks = seq(0, 300, by = 50), limits = c(-20, 250)) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.title.y=element_blank()) + coord_flip()
SV.6.loxPsym.pool.score.plot

pdf("./Notebook/Rplots/OTX2/SV.6.loxPsym.pool.score.pdf", width=6, height=4, useDingbats=FALSE)
SV.6.loxPsym.pool.score.plot
dev.off()

######
# Plot SV distribution

SV.6.loxPsym.SV.diversity.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/003_OTX2_superenhancer_scramble_reanalysis/rds/SV.6.loxPsym.SV.diversity.df.rds")

SV.div.clone.plot <- SV.6.loxPsym.SV.diversity.df %>% filter(exp == "clone") %>%
  filter(variant_name %in% SV.6.loxPsym.clone.variant.filter) %>% 
  filter(variant_name != "no_sv_evidence") %>% 
  mutate(GATE = fct_relevel(gate, "dark", "dim", "medium", "bright")) %>%
  ggplot(aes(x=GATE, y=perc_reads_scaled, fill=variant_name)) +
  geom_bar(stat="identity") + ylab("Percentage of observed reads") + ggtitle("clonal") + ylim(0,40) +
  scale_fill_manual(values = c("#5BBDD7", "#00A08A", "#ED2024", "#F1AD1D", "#F5841F", "#D1D3D4")) +
  theme_bw() + theme(aspect.ratio=2, panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
SV.div.clone.plot

SV.div.pool.plot <- SV.6.loxPsym.SV.diversity.df %>% filter(exp == "pool") %>%
  filter(variant_name %in% SV.6.loxPsym.pool.variant.filter) %>% 
  mutate(variant_name = case_when(variant_name %in% SV.6.loxPsym.clone.variant.filter ~ variant_name, T ~ "others")) %>%
  filter(variant_name != "no sv evidence") %>% 
  mutate(GATE = fct_relevel(gate, "dark", "dim", "medium", "bright")) %>%
  ggplot(aes(x=GATE, y=perc_reads_scaled, fill=variant_name)) +
  geom_bar(stat="identity") + ylab("Percentage of observed reads") + ggtitle("pool transfection") + ylim(0,30) +
  scale_fill_manual(values = c("#5BBDD7", "#00A08A", "#ED2024", "#F1AD1D", "#D1D3D4")) +
  theme_bw() + theme(aspect.ratio=2, panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
SV.div.pool.plot

pdf("./Notebook/Rplots/OTX2/Comp_clonal_pool.pdf", width=6, height=4, useDingbats=FALSE)
ggarrange(SV.div.clone.plot, SV.div.pool.plot, nrow = 1)
dev.off()

# Assess contribution of full deletion to total signal

full.deletion.clone.df <- SV.6.loxPsym.SV.diversity.df %>% filter(exp == "clone") %>%
  filter(variant_name == "∆R12345") %>% pull(perc_reads_scaled)
full.deletion.clone.contribution <- sum(full.deletion.clone.df)/4

full.deletion.pool.df <- SV.6.loxPsym.SV.diversity.df %>% filter(exp == "pool") %>%
  filter(variant_name == "∆R12345") %>% pull(perc_reads_scaled)
full.deletion.pool.contribution <- sum(full.deletion.pool.df)/4
