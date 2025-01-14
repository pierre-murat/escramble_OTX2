
setwd("/lustre/scratch126/gengen/projects/escramble")

cd /lustre/scratch126/gengen/projects/escramble/Data/ONT/240408

scp -r farm:/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW /Users/pm23/Desktop/Projects/EnhancerScramble/Data/ONT/240408
scp -r farm:/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BAMsupp /Users/pm23/Desktop/Projects/EnhancerScramble/Data/ONT/240408

# module load HGI/softpack/groups/escramble/eSCRAMBLE/7

library(dplyr)
library(wesanderson)

############################################################################################
# Load savana SV information

all.savana.LCR.SV.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.scramble.savana.genotype.read.summary.df.rds")

# Load read length and sorting gates information
all.savana.LCR.read.length.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/OTX2.LCR.read.length.df.rds")
nrow(all.savana.LCR.read.length.df)
# 154,360 reads
all.savana.LCR.read.length.format.df <- all.savana.LCR.read.length.df %>% select(read_names = read.id, read.length)

############################################################################################
# Extract bam files related to SVs for assiging genotypes and computing coverages

# WT
# Identify reads with full length with no mapped SVs and remove mapped INV SVs

full.length.reads.df <- all.savana.LCR.read.length.format.df %>% filter(as.numeric(read.length) >= 23137 & as.numeric(read.length) <= 23537) %>% distinct()
# 3372 reads
full.length.savana.LCR.df <- all.savana.LCR.SV.df %>% separate_rows(read_names, sep = "_") %>% filter(read_names %in% full.length.reads.df$read_names)
# 173 reads

WT.savana.reads.df <- full.length.reads.df %>% filter(!read_names %in% full.length.savana.LCR.df$read_names) %>% filter(!read_names %in% INV_AC.savana.reads) %>% filter(!read_names %in% INV_AB.savana.reads) %>% filter(!read_names %in% INV_BC.savana.reads)
WT.savana.reads <- unique(unlist(strsplit(WT.savana.reads.df$read_names, "_"))) # 3193 reads
write.table(WT.savana.reads, "./Data/ONT/240408/BAMsupp/WT.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/WT.reads.txt -o ./BAMsupp/WT.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/WT.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b ./BAMsupp/WT.bam -o ./BW/WT.bw'

# DEL_A-C

DEL_AC.savana.reads.df <- all.savana.LCR.SV.df %>% filter(genotype == "DEL_A-C") %>% separate_rows(read_names, sep = "_") %>% left_join(all.savana.LCR.read.length.format.df)
DEL_AC.savana.reads <- unique(unlist(strsplit(DEL_AC.savana.reads.df$read_names, "_"))) # 109830 reads
write.table(DEL_AC.savana.reads, "./Data/ONT/240408/BAMsupp/DEL_AC.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/DEL_AC.reads.txt -o ./BAMsupp/DEL_AC.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/DEL_AC.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b ./BAMsupp/DEL_AC.bam -o ./BW/DEL_AC.bw'

# DEL_A-B

DEL_AB.savana.reads.df <- all.savana.LCR.SV.df %>% filter(genotype == "DEL_A-B") %>% separate_rows(read_names, sep = "_") %>% left_join(all.savana.LCR.read.length.format.df)
DEL_AB.savana.reads <- unique(DEL_AB.savana.reads.df$read_names, "_") # 16197 reads
write.table(DEL_AB.savana.reads, "./Data/ONT/240408/BAMsupp/DEL_AB.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/DEL_AB.reads.txt -o ./BAMsupp/DEL_AB.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/DEL_AB.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b ./BAMsupp/DEL_AB.bam -o ./BW/DEL_AB.bw'

# INV_A-C-INV_B-C
# Corrected genotype: DEL_A-B_INV_B-C

INV_AC_INV_BC.savana.reads.df <- all.savana.LCR.SV.df %>% filter(genotype == "INV_A-C-INV_B-C") %>% separate_rows(read_names, sep = "_") %>% left_join(all.savana.LCR.read.length.format.df)
INV_AC_INV_BC.savana.reads <- unique(INV_AC_INV_BC.savana.reads.df$read_names, "_") # 3185 reads
write.table(INV_AC_INV_BC.savana.reads, "./Data/ONT/240408/BAMsupp/INV_AC_INV_BC.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/INV_AC_INV_BC.reads.txt -o ./BAMsupp/INV_AC_INV_BC.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/INV_AC_INV_BC.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ./BAMsupp/INV_AC_INV_BC.bam -o ./BW/INV_AC_INV_BC.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ./BAMsupp/INV_AC_INV_BC.bam -o ./BW/INV_AC_INV_BC.rev.bw'

# DEL_B-C

DEL_BC.savana.reads.df <- all.savana.LCR.SV.df %>% filter(genotype == "DEL_B-C") %>% separate_rows(read_names, sep = "_") %>% left_join(all.savana.LCR.read.length.format.df)
DEL_BC.savana.reads <- unique(DEL_BC.savana.reads.df$read_names, "_") # 1470 reads
write.table(DEL_BC.savana.reads, "./Data/ONT/240408/BAMsupp/DEL_BC.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/DEL_BC.reads.txt -o ./BAMsupp/DEL_BC.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/DEL_BC.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None -b ./BAMsupp/DEL_BC.bam -o ./BW/DEL_BC.bw'

# INV_A-C-INV_A-B
# Corrected genotype: INV_A-B_DEL_B-C

INV_AC_INV_AB.savana.reads.df <- all.savana.LCR.SV.df %>% filter(genotype == "INV_A-C-INV_A-B") %>% separate_rows(read_names, sep = "_") %>% left_join(all.savana.LCR.read.length.format.df)
INV_AC_INV_AB.savana.reads <- unique(INV_AC_INV_AB.savana.reads.df$read_names, "_") # 814 reads
write.table(INV_AC_INV_AB.savana.reads, "./Data/ONT/240408/BAMsupp/INV_AC_INV_AB.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/INV_AC_INV_AB.reads.txt -o ./BAMsupp/INV_AC_INV_AB.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/INV_AC_INV_AB.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ./BAMsupp/INV_AC_INV_AB.bam -o ./BW/INV_AC_INV_AB.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ./BAMsupp/INV_AC_INV_AB.bam -o ./BW/INV_AC_INV_AB.rev.bw'

# INV_A-C
# NOPE

INV_AC.savana.reads.df <- LCR.scramble.savana.genotype.df %>% filter(SVTYPE == "INV") %>% filter(SVLEN >= 22000 & genotype == "INV_A-C") %>% left_join(all.savana.LCR.read.length.format.df) %>% filter(read.length > 23000 & read.length < 23700)
INV_AC.savana.reads <- unique(INV_AC.savana.reads.df$read_names, "_") # 30 reads
write.table(INV_AC.savana.reads, "./Data/ONT/240408/BAMsupp/INV_AC.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)
hist(INV_AC.savana.reads.df$read.length, breaks = 100)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/INV_AC.reads.txt -o ./BAMsupp/INV_AC.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/INV_AC.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ./BAMsupp/INV_AC.bam -o ./BW/INV_AC.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ./BAMsupp/INV_AC.bam -o ./BW/INV_AC.rev.bw'

# INV_A-B

INV_AB.savana.reads.df <- LCR.scramble.savana.genotype.df %>% filter(SVTYPE == "INV") %>% filter(genotype == "INV_A-B") %>% left_join(all.savana.LCR.read.length.format.df) %>% filter(read.length > 23000 & read.length < 23700)
INV_AB.savana.reads <- unique(INV_AB.savana.reads.df$read_names, "_") # 44 reads
write.table(INV_AB.savana.reads, "./Data/ONT/240408/BAMsupp/INV_AB.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)
hist(INV_AB.savana.reads.df$read.length, breaks = 100)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/INV_AB.reads.txt -o ./BAMsupp/INV_AB.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/INV_AB.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ./BAMsupp/INV_AB.bam -o ./BW/INV_AB.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ./BAMsupp/INV_AB.bam -o ./BW/INV_AB.rev.bw'

# INV_B-C

INV_BC.savana.reads.df <- LCR.scramble.savana.genotype.df %>% filter(SVTYPE == "INV") %>% filter(genotype == "INV_B-C") %>% left_join(all.savana.LCR.read.length.format.df) %>% filter(read.length > 23000 & read.length < 23700)
INV_BC.savana.reads <- unique(INV_BC.savana.reads.df$read_names, "_") # 65 reads
write.table(INV_BC.savana.reads, "./Data/ONT/240408/BAMsupp/INV_BC.reads.txt", sep="\t", col.names = F, row.names = F, quote = F)
hist(INV_BC.savana.reads.df$read.length, breaks = 100)

bsub -q normal -n 12 -R'select[mem>32000] rusage[mem=32000] span[hosts=1]' -M32000 -o ./BAMsupp/samtools.log -J samtools \
'samtools view -@ 12 -m 2G -bS -N ./BAMsupp/INV_BC.reads.txt -o ./BAMsupp/INV_BC.bam ./BAM/OTX2_LCR_scramble_all_reads.bam
samtools index -@ 12 ./BAMsupp/INV_BC.bam
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand reverse -b ./BAMsupp/INV_BC.bam -o ./BW/INV_BC.fwd.bw
bamCoverage -p 12 -bs 100 --normalizeUsing None --filterRNAstrand forward -b ./BAMsupp/INV_BC.bam -o ./BW/INV_BC.rev.bw'

############################################################################################
# Plot distribution of reads within sorting gates

# Wrapping function for data plotting

LCR.wrap.df <- function(genotype, df) {
  tibble(genotype = rep(genotype, 8)) %>%
    mutate(sample = rep(c("very_dim", "dim", "bright", "very_bright"), 2)) %>% 
    mutate(rep = rep(c(rep(1,4), rep(2,4)), 1)) %>% 
    left_join(df, by = c("genotype", "sample", "rep")) %>% mutate(across(where(is.numeric), ~replace_na(., 0)))
}

# Define number of reference reads

# very.dim.rep.1.read.count <- 4935
# dim.rep.1.read.count <- 8675
# bright.rep.1.read.count <- 5421
# very.bright.rep.1.read.count <- 4455
# 
# very.dim.rep.2.read.count <- 5766
# dim.rep.2.read.count <- 3527
# bright.rep.2.read.count <- 3096
# very.bright.rep.2.read.count <- 5977
# 
# ref.read.count.df <- tibble(sample.id = rep(c("very_dim", "dim", "bright", "very_bright"), 2), replicate = c(rep(1,4), rep(2,4)),
#                             count.ref = c(very.dim.rep.1.read.count, dim.rep.1.read.count, bright.rep.1.read.count, very.bright.rep.1.read.count, very.dim.rep.2.read.count, dim.rep.2.read.count, bright.rep.2.read.count, very.bright.rep.2.read.count))
# 

ref.read.count.df <- all.savana.LCR.read.length.df %>% group_by(sample.id, replicate) %>% summarise(count.ref = dplyr::n())

# WT

WT.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% WT.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "WT") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
WT.LCR.wrap.df <- LCR.wrap.df("WT", WT.LCR.df)

# DEL_A-C

DEL_AC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% DEL_AC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "DEL_A-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
DEL_AC.LCR.wrap.df <- LCR.wrap.df("DEL_A-C", DEL_AC.LCR.df)

# DEL_A-B

DEL_AB.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% DEL_AB.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "DEL_A-B") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
DEL_AB.LCR.wrap.df <- LCR.wrap.df("DEL_A-B", DEL_AB.LCR.df)

# DEL_B-C

DEL_BC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% DEL_BC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "DEL_B-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
DEL_BC.LCR.wrap.df <- LCR.wrap.df("DEL_B-C", DEL_BC.LCR.df)

# INV_A-C

INV_AC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_AC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_A-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
INV_AC.LCR.wrap.df <- LCR.wrap.df("INV_A-C", INV_AC.LCR.df)

# INV_A-B

INV_AB.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_AB.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_A-B") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
INV_AB.LCR.wrap.df <- LCR.wrap.df("INV_A-B", INV_AB.LCR.df)

# INV_B-C

INV_BC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_BC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_B-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
INV_BC.LCR.wrap.df <- LCR.wrap.df("INV_B-C", INV_BC.LCR.df)

# INV_A-C-INV_B-C
# Corrected genotype: DEL_A-B_INV_B-C

DELAB_INVBC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_AC_INV_BC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_A-C-INV_B-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
DELAB_INVBC.LCR.wrap.df <- LCR.wrap.df("INV_A-C-INV_B-C", DELAB_INVBC.LCR.df) %>% mutate(genotype = "DEL_A-B_INV_B-C")

# INV_A-C-INV_B-C
# Corrected genotype: DEL_A-B_INV_B-C

DELAB_INVBC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_AC_INV_BC.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_A-C-INV_B-C") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
DELAB_INVBC.LCR.wrap.df <- LCR.wrap.df("INV_A-C-INV_B-C", DELAB_INVBC.LCR.df) %>% mutate(genotype = "DEL_A-B_INV_B-C")

# INV_A-C-INV_A-B
# Corrected genotype: INV_A-B_DEL_B-C

INVAB_DELBC.LCR.df <- all.savana.LCR.read.length.df %>% filter(read.id %in% INV_AC_INV_AB.savana.reads) %>% group_by(sample.id, replicate) %>% summarise(count=dplyr::n()) %>% left_join(ref.read.count.df, by = c("sample.id", "replicate")) %>% mutate(perc = (count/count.ref)*100) %>% mutate(genotype = "INV_A-C-INV_A-B") %>% dplyr::select(genotype, sample = sample.id, rep = replicate, count, perc)
INVAB_DELBC.LCR.wrap.df <- LCR.wrap.df("INV_A-C-INV_A-B", INVAB_DELBC.LCR.df) %>% mutate(genotype = "INV_A-B_DEL_B-C")

############################################################################################
# Combine all data

LCR.savana.SV.summary <- rbind(WT.LCR.wrap.df, DEL_AC.LCR.wrap.df, DEL_AB.LCR.wrap.df, DEL_BC.LCR.wrap.df, INV_AC.LCR.wrap.df, INV_AB.LCR.wrap.df, INV_BC.LCR.wrap.df, DELAB_INVBC.LCR.wrap.df, DELAB_INVBC.LCR.wrap.df, INVAB_DELBC.LCR.wrap.df)
saveRDS(LCR.savana.SV.summary, "/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.SV.summary.rds")

############################################################################################
# Plot data

LCR.savana.SV.summary <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.SV.summary.rds")

LCR.palette <- wes_palette("Zissou1", length(unique(LCR.savana.SV.summary$genotype)), type = "continuous")
LCR.savana.SV.plot <- LCR.savana.SV.summary %>% 
  mutate(sample = case_when(sample == "very_dim" ~ "dark", sample == "dim" ~ "dim", sample == "bright" ~ "medium", sample == "very_bright" ~ "bright")) %>% 
  mutate(SAMPLE = fct_relevel(sample, "dark", "dim", "medium", "bright")) %>%
  ggplot(aes(x=SAMPLE, y=perc, group=1, col = genotype)) +
  geom_point(aes(x = SAMPLE, y = perc), size = 1) +
  geom_point(stat='summary', fun.y=sum, size = 2) +
  geom_line(stat='summary', fun.y=sum) +
  ylab("Percent of events in sorting gate") +
  scale_colour_manual(values=c(LCR.palette)) +
  facet_wrap(~ genotype, scales = "free", ncol = 5) +
  theme_bw() + theme(aspect.ratio=1, strip.background = element_blank(), panel.grid.major = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none")
LCR.savana.SV.plot

############################################################################################
# Format coverage data

# Define region to plot
chr <- "chr14"
from <- 56884500
to <- 56907900

# Prepare empty table for filling missing values
LCR.coord.df <- tibble(start = seq(from+1, to+1, 100))

# Deletion

WT.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/WT.bw")
WT.df <- as.data.frame(WT.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "WT")

DEL_AC.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/DEL_AC.bw"  )
DEL_AC.df <- as.data.frame(DEL_AC.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "DEL_A-C")

DEL_AB.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/DEL_AB.bw"  )
DEL_AB.df <- as.data.frame(DEL_AB.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "DEL_A-B")

DEL_BC.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/DEL_BC.bw"  )
DEL_BC.df <- as.data.frame(DEL_BC.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "DEL_B-C")

# Deletion and/or inversion

INVAC.fwd.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC.fwd.bw")
INVAC.rev.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC.rev.bw")
INVAC.fw.df <- as.data.frame(INVAC.fwd.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "INV_A-C", strand = "+")
INVAC.rev.df <- as.data.frame(INVAC.rev.bw) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(score = -score, variant = "INV_A-C", strand = "-")
INVAC.df <- rbind(INVAC.fw.df, INVAC.rev.df)

INVAB.fwd.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AB.fwd.bw")
INVAB.rev.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AB.rev.bw")
INVAB.fw.df <- as.data.frame(INVAB.fwd.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "INV_A-B", strand = "+")
INVAB.rev.df <- as.data.frame(INVAB.rev.bw) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(score = -score, variant = "INV_A-B", strand = "-")
INVAB.df <- rbind(INVAB.fw.df, INVAB.rev.df)

INVBC.fwd.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_BC.fwd.bw")
INVBC.rev.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_BC.rev.bw")
INVBC.fw.df <- as.data.frame(INVBC.fwd.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "INV_B-C", strand = "+")
INVBC.rev.df <- as.data.frame(INVBC.rev.bw) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(score = -score, variant = "INV_B-C", strand = "-")
INVBC.df <- rbind(INVBC.fw.df, INVBC.rev.df)

DELAB_INVBC.fwd.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC_INV_BC.fwd.bw")
DELAB_INVBC.rev.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC_INV_BC.rev.bw")
DELAB_INVBC.fw.df <- as.data.frame(DELAB_INVBC.fwd.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "DEL_A-B_INV_B-C", strand = "+")
DELAB_INVBC.rev.df <- as.data.frame(DELAB_INVBC.rev.bw) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(score = -score, variant = "DEL_A-B_INV_B-C", strand = "-")
DELAB_INVBC.df <- rbind(DELAB_INVBC.fw.df, DELAB_INVBC.rev.df)

INVAB_DELBC.fwd.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC_INV_AB.fwd.bw")
INVAB_DELBC.rev.bw <- import("/lustre/scratch126/gengen/projects/escramble/Data/ONT/240408/BW/INV_AC_INV_AB.rev.bw")
INVAB_DELBC.fw.df <- as.data.frame(INVAB_DELBC.fwd.bw) %>% filter(seqnames == chr & start >= from & end <= to) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(variant = "INV_A-B_DEL_B-C", strand = "+")
INVAB_DELBC.rev.df <- as.data.frame(INVAB_DELBC.rev.bw) %>% mutate_at(ncol(.), funs(ifelse(row_number() == n(), 0, .))) %>% right_join(LCR.coord.df) %>% arrange(start) %>% fill(score) %>% mutate(score = -score, variant = "INV_A-B_DEL_B-C", strand = "-")
INVAB_DELBC.df <- rbind(INVAB_DELBC.fw.df, INVAB_DELBC.rev.df)

############################################################################################
# Combine data

# Deletion

LCR.savana.deletion.coverage.df <- rbind(WT.df, DEL_AC.df, DEL_AB.df, DEL_BC.df)
saveRDS(LCR.savana.deletion.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.deletion.coverage.df.rds")

LCR.savana.inversion.coverage.df <- rbind(INVAC.df, INVAB.df, INVBC.df, DELAB_INVBC.df, INVAB_DELBC.df)
saveRDS(LCR.savana.inversion.coverage.df, "/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.inversion.coverage.df.rds")

############################################################################################
# Plot coverage data

LCR.savana.deletion.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.deletion.coverage.df.rds")
LCR.savana.inversion.coverage.df <- readRDS("/lustre/scratch126/gengen/projects/escramble/Notebook/004_OTX2_LCR_scramble/rds/LCR.savana.inversion.coverage.df.rds")

LCR.palette <- wes_palette("Zissou1", length(unique(LCR.savana.SV.summary$genotype)), type = "continuous")

# Deletion

LCR.deletion.palette <- wes_palette("Zissou1", length(unique(LCR.savana.deletion.coverage.df$variant)), type = "continuous")
LCR.savana.deletion.coverage.plot <- LCR.savana.deletion.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = variant)) +
  geom_area() +
  geom_line(col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(LCR.deletion.palette)) +
  xlim(from/1e6, to/1e6) +
  geom_vline(xintercept=56884634/1e6, linetype="dashed") +
  geom_vline(xintercept=56898376/1e6, linetype="dashed") +
  geom_vline(xintercept=56907608/1e6, linetype="dashed") +
  facet_wrap(~ variant, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), strip.background = element_blank(), legend.position="none")
LCR.savana.deletion.coverage.plot

LCR.inversion.palette <- wes_palette("Zissou1", length(unique(LCR.savana.inversion.coverage.df$variant)), type = "continuous")
LCR.savana.inversion.coverage.plot <- LCR.savana.inversion.coverage.df %>% ggplot(aes(x = start/1e6, y = score, group=1, fill = variant)) +
  geom_area(aes(group = strand)) +
  geom_line(aes(group = strand), col = "black", linewidth = 0.25) +
  ylab("Read count") + xlab("Coordinate [Mb]") +
  scale_fill_manual(values=c(LCR.inversion.palette)) +
  xlim(from/1e6, to/1e6) +
  geom_hline(yintercept=0, linetype="solid") +
  geom_vline(xintercept=56884634/1e6, linetype="dashed") +
  geom_vline(xintercept=56898376/1e6, linetype="dashed") +
  geom_vline(xintercept=56907608/1e6, linetype="dashed") +
  facet_wrap(~ variant, scales = "free", ncol = 3) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), strip.background = element_blank(), legend.position="none")
LCR.savana.inversion.coverage.plot










