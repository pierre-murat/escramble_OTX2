
setwd("/lustre/scratch126/gengen/projects_v2/escramble/001_007_OTX2")

#==========================================#
# Load libraries

packages <- c('GenomicRanges', 'rtracklayer', 'dplyr', 'purrr', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome')
for(p in packages){ library(p, character.only = TRUE) }

#==========================================#
# Define tiles

OTX2_tiles <- GRanges(
  seqnames = rep("chr14", 12),
  ranges = IRanges(
    start = c(56886000,56888000,56890000,56892000,56894000,56896000,56898000,56900000,56902000,56904000, 56846374, 56855162),
    end = c(56887999,56889999,56891999,56893999,56895999,56897999,56899999,56901999,56903999,56905999, 56848373, 56857161)),
  strand = "*",
  tile_id = c(paste0("tile_", 1:10), paste0("up_", 1:2)))
OTX2_tiles
OTX2_tiles.seq <- DNAStringSet(getSeq(BSgenome.Hsapiens.UCSC.hg38, OTX2_tiles))
names(OTX2_tiles.seq) <- paste0("OTX2_", OTX2_tiles$tile_id)
writeXStringSet(OTX2_tiles.seq, filepath = "./009_Pilot_design/fasta/OTX2.superenhance.tile.seq.fasta", format = "fasta")

# Generate fasta files reporting the super-enhancer domains analysed in the OTX2 manuscript
# pegRNA positions
# chr14	56884634	56884635
# chr14	56898376	56898377
# chr14	56907608	56907609
# chr14	56893129	56893130
# chr14	56901292	56901293
# chr14	56903183	56903184
# chr14	56904792	56904793
pos <- sort(c(56898376, 56901292, 56903183, 56904792, 56907608))
OTX2_domains <- GRanges(
  seqnames = "chr14",
  ranges = IRanges(
    start = pos[-length(pos)],
    end   = pos[-1]
  ),
  strand = "*",
  names = c("R5_1", "R5_2", "R5_3", "R5_4")
)
OTX2_domains.seq <- DNAStringSet(getSeq(BSgenome.Hsapiens.UCSC.hg38, OTX2_domains))
names(OTX2_domains.seq) <- OTX2_domains$names
writeXStringSet(OTX2_domains.seq, filepath = "./009_Pilot_design/fasta/OTX2.superenhance.domain.seq.fasta", format = "fasta")

#==========================================#
# Recover OTX2 enhancer candidates sequences

# Load HAP1 DHS peaks

HAP1.DHS.peak.df <- read.table("/lustre/scratch125/gengen/teams_v2/parts/resources/HAP1/ENCODE/HAP1.DHS.seq.peaks.bed") %>% 
  dplyr::select(V1, V2, V3)
colnames(HAP1.DHS.peak.df) <- c("seqnames", "start", "end")
HAP1.DHS.peak.gr <- HAP1.DHS.peak.df %>% GRanges()

# Subset OTX2 peaks

OTX2.GOI.gr <- GRanges(seqnames = "chr14",
                   ranges = IRanges(start = 56884600,
                                    end = 56907640),
                   strand = "*",
                   gene = "OTX2")
OTX2.DHS.peak.gr <- subsetByOverlaps(HAP1.DHS.peak.gr, OTX2.GOI.gr)

# Export bed file
write.table(as.data.frame(OTX2.DHS.peak.gr), "./Ressources/bed/OTX2.DHS.peak.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Complete with manually curated peaks

add.peak.gr <- GRanges(seqnames = "chr14",
                       ranges = IRanges(start = c(56895217, 56899861, 56901448),
                                        end = c(56895688, 56900210, 56901760)),
                       strand = "*")
OTX2.DHS.peak.2.gr <- sort(c(OTX2.DHS.peak.gr, add.peak.gr))

# Export bed file
write.table(as.data.frame(OTX2.DHS.peak.2.gr), "./Ressources/bed/OTX2.DHS.peak.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Resize all peaks to 450 nt
OTX2.DHS.peak.resize.gr <- resize(OTX2.DHS.peak.2.gr, 450, fix="center")

# Export bed file
write.table(as.data.frame(OTX2.DHS.peak.resize.gr), "./Ressources/bed/OTX2.DHS.peak.bed", sep="\t", col.names = F, row.names = F, quote = F)

# 12 enhancers

# Recover sequences
OTX2.DHS.peak.seq <- getSeq(Hsapiens, OTX2.DHS.peak.resize.gr)
names(OTX2.DHS.peak.seq) <- c("R4-1_A", "R4-1_B", "R4-1_C", "R4-2_A", "R4-2_B", "R4-2_C", "R5-1", "R5-2_A", "R5-2_B", "R5-2_C", "R5-3", "R5-4")
OTX2.DHS.peak.resize.gr$tile_id <- c("R4-1_A", "R4-1_B", "R4-1_C", "R4-2_A", "R4-2_B", "R4-2_C", "R5-1", "R5-2_A", "R5-2_B", "R5-2_C", "R5-3", "R5-4")

# The peak of interest is R5-2_C and overlap chr14:56,902,255-56,902,721
OTX2.DHS.peak.seq[["R5-2_C"]]
R52C.seq <- as.character(OTX2.DHS.peak.seq[["R5-2_C"]])

# Combine with tiles and export bed file for visualisation
all <- c(OTX2_tiles, OTX2.DHS.peak.resize.gr)
names(all) <- all$tile_id
export(all, "./000_OTX2_enhancer_sequences/bed/OTX2_tiles.bed", format = "bed")

#==========================================#
# Mutate R5-2_C TF binding sites

# Define TFBS (motif +- 1 bp)

FOXA1_1 <- "TCTTTATT"
FOXA1_1_mut <- "TT"

FOXA1_2 <- "GTAAATAC"
FOXA1_2_mut <- "GC"

FOXA1_3 <- "ATTTTTAC"
FOXA1_3_mut <- "AC"

FOXA1_4 <- "AAAGATGT"
FOXA1_4_mut <- "AT"

FOXA1_5 <- "CTGTTTAC"
FOXA1_5_mut <- "CC"

FOXA1_6 <- "TTATTAAC"
FOXA1_6_mut <- "TC"

YY1 <- "GCCATTTTC"
YY1_mut <- "GC"

RFX5 <- "TTTCCATAGC"
RFX5_mut <- "TC"

MYB_1 <- "TCAGTCGC"
MYB_1_mut <- "TC"

MYB_2 <- "CCAACTGG"
MYB_2_mut <- "CG"

MYB_3 <- "GCAACTGC"
MYB_3_mut <- "GC"

MYB_4 <- "TCAGTTGA"
MYB_4_mut <- "TA"

# Prepare sequences
# to individually mutate FOXA1-1, FOXA1-2, MYB-1, FOXA1-3, FOXA1-4, FOXA1-5, MYB-2 (7 sequences)
# mutate at the same time FOXA1-1, FOXA1-2, MYB-1, FOXA1-3, FOXA1-4, FOXA1-5, MYB-2 (1 sequence)
# mutate all FOXA1 sites (1 sequence)
# mutate all MYB sites (1 sequence)
# mutate RFX5 site (1 sequence)
# mutate YY1 site (1 sequence)

# individual motifs

R52C.FOXA1_1.seq <- gsub(FOXA1_1, FOXA1_1_mut, R52C.seq)
R52C.FOXA1_2.seq <- gsub(FOXA1_2, FOXA1_2_mut, R52C.seq)
R52C.MYB_1.seq <- gsub(MYB_1, MYB_1_mut, R52C.seq)
R52C.FOXA1_3.seq <- gsub(FOXA1_3, FOXA1_3_mut, R52C.seq)
R52C.FOXA1_4.seq <- gsub(FOXA1_4, FOXA1_4_mut, R52C.seq)
R52C.FOXA1_5.seq <- gsub(FOXA1_5, FOXA1_5_mut, R52C.seq)
R52C.MYB_2.seq <- gsub(MYB_2, MYB_2_mut, R52C.seq)

# core motifs

R52C.core.seq <-
  R52C.seq |>
  gsub(FOXA1_1, FOXA1_1_mut, x = _) |>
  gsub(FOXA1_2, FOXA1_2_mut, x = _) |>
  gsub(MYB_1, MYB_1_mut, x = _) |>
  gsub(FOXA1_3, FOXA1_3_mut, x = _) |>
  gsub(FOXA1_4, FOXA1_4_mut, x = _) |>
  gsub(FOXA1_5, FOXA1_5_mut, x = _) |>
  gsub(MYB_2, MYB_2_mut, x = _)
  
# all FOXA1 motifs

R52C.all.FOXA1.seq <-
  R52C.seq |>
  gsub(FOXA1_1, FOXA1_1_mut, x = _) |>
  gsub(FOXA1_2, FOXA1_2_mut, x = _) |>
  gsub(FOXA1_3, FOXA1_3_mut, x = _) |>
  gsub(FOXA1_4, FOXA1_4_mut, x = _) |>
  gsub(FOXA1_5, FOXA1_5_mut, x = _) |>
  gsub(FOXA1_6, FOXA1_6_mut, x = _)

# all MYB motifs

R52C.all.MYB.seq <-
  R52C.seq |>
  gsub(MYB_1, MYB_1_mut, x = _) |>
  gsub(MYB_2, MYB_2_mut, x = _) |>
  gsub(MYB_3, MYB_3_mut, x = _) |>
  gsub(MYB_4, MYB_4_mut, x = _)

# RFX5 site 

R52C.RFX5.seq <- gsub(RFX5, RFX5_mut, R52C.seq)

# YY1 site 

R52C.YY1.seq <- gsub(YY1, YY1_mut, R52C.seq)

# Combine all sequences

TFBS.seqs <- c(
  `R5-2_C_FOXA1_1`             = R52C.FOXA1_1.seq,
  `R5-2_C_FOXA1_2`             = R52C.FOXA1_2.seq,
  `R5-2_C_MYB_1`               = R52C.MYB_1.seq,
  `R5-2_C_FOXA1_3`             = R52C.FOXA1_3.seq,
  `R5-2_C_FOXA1_4`             = R52C.FOXA1_4.seq,
  `R5-2_C_FOXA1_5`             = R52C.FOXA1_5.seq,
  `R5-2_C_MYB_2`               = R52C.MYB_2.seq,
  `R5-2_C_core_TFBS`           = R52C.core.seq,
  `R5-2_C_all_FOXA1`           = R52C.all.FOXA1.seq,
  `R5-2_C_all_MYB`             = R52C.all.MYB.seq,
  `R5-2_C_RFX5`                = R52C.RFX5.seq,
  `R5-2_C_YY1`                 = R52C.YY1.seq
)
R52C.TFBS.seqs <- DNAStringSet(TFBS.seqs)
#
OTX2.enhancer.seqs <- c(OTX2.DHS.peak.seq, TFBS.seqs)
#
length(unique(as.character(OTX2.enhancer.seqs))) == length(OTX2.enhancer.seqs)
# TRUE
length(OTX2.enhancer.seqs)
# 24
sort(unique(nchar(as.character(OTX2.enhancer.seqs))))
# 408 414 426 442 443 444 450

#==========================================#
# Export

writeXStringSet(OTX2.enhancer.seqs, "./Ressources/fasta/OTX2.enhancer.seqs.fasta")

names(OTX2.enhancer.seqs)






