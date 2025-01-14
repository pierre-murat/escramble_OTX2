# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(rstatix)
library(GenomicRanges)
library(ggrastr)

theme_sv <-   theme_bw(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black"))

loxPsym_sites = tibble(precise = c(56840403, 56851953, 56864718, 56884634, 56898376, 56907608)) %>% mutate(chr = "chr14", start = precise - 100, end = precise + 100) %>% GRanges()

variant_annotation <- tibble(
  name = c("∆R4", "∆R5", "∆R45", "iR4", "iR5", "iR45", "∆R4+iR5", "iR4+∆R5", "wildtype"),
  variant = c("DEL13742", "DEL9232", "DEL22974", "INV13742", "INV9232", "INV22974", "INV22974_INV9232", "INV13742_INV22974", "WT"),
  order = c("4_DEL13742", "2_DEL9232", "DEL22974", "INV13742", "INV9232", "INV22974", "5_INV22974_INV9232", "3_INV13742_INV22974", "1_WT"))

read_names_loxP6 <- read_tsv("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/read_names/read_names_loxP6.tsv")
read_names_LCR <- read_tsv("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/read_names/read_names_LCR.tsv")


# Read the BED12 file
process_bed12 <- function(path, modality, read_names) {
  input <- read_tsv(path, col_names = FALSE)
  colnames(input) <- c("chr", "read_start", "read_end", "fiber", "score", "strand", "thick_start", "thick_end", "itemRgb", "block_count", "block_size", "block_starts")
  
  output <- input %>%
    mutate(
      block_size = str_split(block_size, ","),
      block_starts = str_split(block_starts, ",")
    ) %>%
    unnest(c(block_size, block_starts)) %>%
    mutate(
      block_size = as.numeric(block_size),
      block_starts = as.numeric(block_starts)
    ) %>% filter(block_size > 0) %>%
    dplyr::select(chr, read_start, read_end, fiber, strand, block_size, block_starts) %>%
    mutate(start = read_start + block_starts, end = start + block_size, read_length = read_end - read_start, modality = modality) %>% ungroup() %>%
    group_by(fiber) %>% mutate(perc_m6A = dplyr::n()/read_length*100*4) %>% ungroup() %>%
    left_join(read_names, by = "fiber") %>%
    dplyr::select(-block_size, -block_starts)
}

process_bed12_mCpG <- function(path, modality, read_names) {
  input <- read_tsv(path, col_names = FALSE)
  colnames(input) <- c("chr", "read_start", "read_end", "fiber", "score", "strand", "thick_start", "thick_end", "itemRgb", "block_count", "block_size", "block_starts")
  
  output <- input %>%
    mutate(
      block_size = str_split(block_size, ","),
      block_starts = str_split(block_starts, ",")
    ) %>%
    unnest(c(block_size, block_starts)) %>%
    mutate(
      block_size = as.numeric(block_size),
      block_starts = as.numeric(block_starts)
    ) %>% filter(block_size > 0) %>%
    dplyr::select(chr, read_start, read_end, fiber, strand, block_size, block_starts) %>%
    mutate(start = read_start + block_starts, end = start + block_size, read_length = read_end - read_start, modality = modality) %>% ungroup() %>%
    group_by(fiber) %>% mutate(perc_mCpG = dplyr::n()/read_length*100) %>% ungroup() %>%
    left_join(read_names, by = "fiber") %>%
    dplyr::select(-block_size, -block_starts)
}


LCR_dark_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dark_R2_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_dim_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dim_R2_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_medium_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_medium_R1_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_medium_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_medium_R2_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_bright_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_bright_R1_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_bright_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_bright_R2_m6A.bed", "m6A", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)

#trans_brilliant_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_brilliant_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "brilliant") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_bright_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_bright_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_bright_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_bright_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dim_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dim_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dim_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dim_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_medium_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_medium_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_medium_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_medium_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dark_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dark_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dark_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dark_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)

loxP6_dark_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dark_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dark_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dark_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dim_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dim_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dim_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dim_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_medium_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_medium_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_medium_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_medium_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_bright_R1_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_bright_R1_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_bright_R2_m6a <- process_bed12("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_bright_R2_m6A.bed", "m6A", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)

LCR_dark_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dark_R1_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_dark_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dark_R2_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_dim_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dim_R1_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_dim_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_dim_R2_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_medium_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_medium_R1_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_medium_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_medium_R2_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_bright_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_bright_R1_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
LCR_bright_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_LCR_bright_R2_mCpG.bed", "mCpG", read_names_LCR) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)

trans_dark_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dark_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dark_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dark_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dim_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dim_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_dim_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_dim_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_medium_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_medium_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_medium_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_medium_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_bright_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_bright_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
trans_bright_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_trans_bright_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)

loxP6_dark_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dark_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dark_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dark_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dark") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dim_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dim_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_dim_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_dim_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "dim") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_medium_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_medium_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_medium_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_medium_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "medium") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_bright_R1_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_bright_R1_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)
loxP6_bright_R2_mCpG <- process_bed12_mCpG("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/bed/OTX2_loxP6_bright_R2_mCpG.bed", "mCpG", read_names_loxP6) %>% mutate(name = ifelse(is.na(name), "No sv evidence", name), SVLEN = ifelse(is.na(SVLEN), 0, SVLEN), gate = "bright") %>% mutate(length_for_fl = 68897-1000-SVLEN)

loxP6_m6A <- filter(bind_rows(trans_dark_R1_m6a, trans_dark_R2_m6a, trans_dim_R1_m6a, trans_dim_R2_m6a, trans_medium_R1_m6a, trans_medium_R2_m6a, trans_bright_R1_m6a, trans_bright_R2_m6a, loxP6_dark_R1_m6a, loxP6_dark_R2_m6a, loxP6_dim_R1_m6a, loxP6_dim_R2_m6a, loxP6_medium_R1_m6a, loxP6_medium_R2_m6a, loxP6_bright_R1_m6a, loxP6_bright_R2_m6a), name %in% c("∆R234", "∆R1234", "∆R2345", "∆R34", "∆R23"), perc_m6A < 8)
LCR_m6A <- bind_rows(LCR_dark_R2_m6a, LCR_dim_R2_m6a, LCR_medium_R1_m6a, LCR_medium_R2_m6a, LCR_bright_R1_m6a, LCR_bright_R2_m6a) %>% filter(perc_m6A < 8)

LCR_m6A_WT <- LCR_m6A %>% filter(name == "wildtype" | is.na(variant) & read_length > 24000, read_length < 26000)
loxP6_m6A_WT <- filter(bind_rows(trans_dark_R1_m6a, trans_dark_R2_m6a, trans_dim_R1_m6a, trans_dim_R2_m6a, trans_medium_R1_m6a, trans_medium_R2_m6a, trans_bright_R1_m6a, trans_bright_R2_m6a, loxP6_dark_R1_m6a, loxP6_dark_R2_m6a, loxP6_dim_R1_m6a, loxP6_dim_R2_m6a, loxP6_medium_R1_m6a, loxP6_medium_R2_m6a, loxP6_bright_R1_m6a, loxP6_bright_R2_m6a), name == c("No sv evidence"), read_length > length_for_fl, perc_m6A < 8)

loxP6_mCpG <- bind_rows(filter(trans_medium_R1_mCpG, name %in% c("∆R234", "∆R1234", "∆R2345", "∆R34", "∆R23")), filter(trans_medium_R1_mCpG, name == c("No sv evidence"), read_length > length_for_fl))
LCR_mCpG <- bind_rows(LCR_dark_R1_mCpG, LCR_dark_R2_mCpG, LCR_dim_R1_mCpG, LCR_dim_R2_mCpG, LCR_medium_R1_mCpG, LCR_medium_R2_mCpG, LCR_bright_R1_mCpG, LCR_bright_R2_mCpG)

LCR_mCpG_WT <- LCR_mCpG %>% filter(name == "wildtype" | is.na(variant) & read_length > 24000, read_length < 26000)
loxP6_mCpG_WT <- filter(bind_rows(trans_dim_R1_mCpG, trans_dim_R2_mCpG, trans_medium_R1_mCpG, trans_medium_R2_mCpG, loxP6_dark_R1_mCpG, loxP6_dark_R2_mCpG, loxP6_dim_R1_mCpG, loxP6_dim_R2_mCpG, loxP6_medium_R1_mCpG, loxP6_medium_R2_mCpG, loxP6_bright_R1_mCpG, loxP6_bright_R2_mCpG), name == c("No sv evidence"), read_length > length_for_fl)

# (1) Plot the m6A methylation differences across sorting gates
# mean levels
p <- LCR_m6A_WT %>% mutate(sample = paste(gate, replicate, sep = "_")) %>% group_by(gate, fiber) %>% summarise(perc_m6A = mean(perc_m6A)) %>% 
  ggplot(aes(x = factor(gate, levels = c("bright", "medium", "dim", "dark")), y = perc_m6A, fill = factor(gate, levels = c("bright", "medium", "dim", "dark")))) +
  geom_boxplot(show.legend = F, linewidth = 0.15, outlier.size = 0.1, outlier.stroke = 0.1) +
  scale_fill_manual("Sorting gate", values = c( "#008080", "#800000", "#808000", "#6A5ACD")) +
  labs(x = "Sorting gate", y = "% of adenines methylated") +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/m6A_gates_LCR_boxplot.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2.5, "cm")))


p <- loxP6_m6A_WT %>% filter(perc_m6A < 8) %>% group_by(gate, fiber) %>% summarise(perc_m6A = mean(perc_m6A)) %>% 
  ggplot(aes(x = factor(gate, levels = c("bright", "medium", "dim", "dark")), y = perc_m6A, fill = factor(gate, levels = c("bright", "medium", "dim", "dark")))) +
  geom_boxplot(show.legend = F, linewidth = 0.15, outlier.size = 0.1, outlier.stroke = 0.1) +
  scale_fill_manual("Sorting gate", values = c( "#008080", "#800000", "#808000", "#6A5ACD")) +
  labs(x = "Sorting gate", y = "% of adenines methylated") +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/m6A_gates_loxP6_boxplot.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2.5, "cm")))


levels = LCR_m6A_WT %>% mutate(gate = factor(gate, levels = c("bright", "medium", "dim", "dark"))) %>% arrange(desc(gate), read_length) %>% `$`(fiber) %>% unique()
p <- LCR_m6A_WT %>% mutate(start = start/1000000) %>%
  ggplot(aes(x = start, y = factor(fiber, levels = levels), col = factor(gate, levels = c("bright", "medium", "dim", "dark")))) +
  rasterize(geom_point(size = 0.15, shape = 16, stroke = 0, alpha = 0.2, show.legend = F), dpi = 400) +  # Rasterize the points
  coord_cartesian(xlim = c(56.885800, 56.907500)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual("Sorting gate", values = c( "#008080", "#800000", "#808000", "#6A5ACD")) +
  labs(x = "Position on chromosome 14", y = "Sequencing reads") +
  theme_sv +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/m6A_gates_LCR.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(3, "cm"), height=unit(3, "cm")))


levels = loxP6_m6A_WT %>% mutate(gate = factor(gate, levels = c("bright", "medium", "dim", "dark"))) %>% arrange(desc(gate), read_length) %>% `$`(fiber) %>% unique()
p <- loxP6_m6A_WT %>% mutate(start = start/1000000) %>%
  ggplot(aes(x = start, y = factor(fiber, levels = levels), col = factor(gate, levels = c("bright", "medium", "dim", "dark")))) +
  rasterize(geom_point(size = 0.15, shape = 16, stroke = 0, alpha = 0.2, show.legend = F), dpi = 400) +  # Rasterize the points
  coord_cartesian(xlim = c(56.840593, 56.907500)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual("Sorting gate", values = c( "#008080", "#800000", "#808000", "#6A5ACD")) +
  labs(x = "Position on chromosome 14", y = "Sequencing reads") +
  theme_sv +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/m6A_gates_loxP6.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(9.249811, "cm"), height=unit(3, "cm")))



# (2) Plot fiber seq data for different structural variants =====
LCR_variants <- bind_rows(LCR_dark_R2_m6a, LCR_dim_R2_m6a, LCR_medium_R1_m6a, LCR_medium_R2_m6a, LCR_bright_R1_m6a, LCR_bright_R2_m6a,
                          trans_bright_R1_m6a, trans_bright_R2_m6a, trans_dim_R1_m6a, trans_dim_R2_m6a, trans_medium_R1_m6a, trans_medium_R2_m6a, trans_dark_R1_m6a, trans_dark_R2_m6a,
                          loxP6_dark_R1_m6a, loxP6_dark_R2_m6a, loxP6_dim_R1_m6a, loxP6_dim_R2_m6a, loxP6_medium_R1_m6a, loxP6_medium_R2_m6a, loxP6_bright_R1_m6a, loxP6_bright_R2_m6a) %>%
  filter(name %in% c("∆R5", "∆R4", "∆R4+iR5", "iR4+∆R5", "∆R1234", "∆R123"), perc_m6A < 8, read_length > 9000, start != read_start, end != read_end) %>% # the software always put a modified at read starts and end
  bind_rows(mutate(LCR_m6A_WT, name = "wildtype"), mutate(loxP6_m6A_WT, name = "wildtype"))


levels = LCR_variants %>% mutate(name = factor(name, levels = c("wildtype", "∆R5", "iR4+∆R5", "∆R4", "∆R4+iR5", "∆R123", "∆R1234"))) %>% arrange(desc(name), read_length) %>% `$`(fiber) %>% unique()
p <- LCR_variants %>% filter(start > 56885800, start < 56907500) %>% mutate(start = start/1000000) %>%
  ggplot(aes(x = start, y = factor(fiber, levels = levels), col = factor(name, levels = c("wildtype", "∆R5", "iR4+∆R5", "∆R4", "∆R4+iR5", "∆R123", "∆R1234")))) +
  rasterize(geom_point(size = 0.1, shape = 16, stroke = 0, alpha = 0.2, show.legend = F), dpi = 400) +  # Rasterize the points
  coord_cartesian(xlim = c(56.885800, 56.907500)) +
  scale_color_manual("Variants", values = c("#808080", "#6495ed", "#4169e1", "#a52a2a", "#f08080", "#55AA55", "darkgreen")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Position on chromosome 14", y = "Sequencing reads") +
  theme_sv +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/fiber_LCR.pdf", egg::set_panel_size(p, width=unit(4, "cm"), height=unit(5, "cm")))


# Quantify on enhancers
DNase_peaks <- read_tsv("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/OTX2_locus/LCR/DNase_ENCFF289MVI.bed", col_names = c("chr", "start", "end", "emssion", "X5", "X6", "intensity")) %>%
  filter(chr == "chr14", start > 56840593, start < 56909300, intensity > 1) %>% dplyr::select(chr, start, end, intensity) %>% mutate(width = end - start)
DNase_peaks <- read_tsv("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/OTX2_locus/LCR/DNase_ENCFF289MVI.bed", col_names = c("chr", "start", "end", "emssion", "X5", "X6", "intensity")) %>%
  filter(chr == "chr14", start > 56840593, start < 56909300, intensity > 2.5) %>% dplyr::select(chr, start, end, intensity) %>% mutate(width = end - start)
DNase_peaks$AT = c(309, 152, 178, 143, 130)
DNase_peaks$region_name = c("DHS1", "DHS2", "DHS3", "DHS4", "DHS5")
write_tsv(DNase_peaks, "/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/OTX2_locus/LCR/DNase_peaks_OTX2.tsv")

p <- DNase_peaks %>% filter(intensity > 1) %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 50 , ymax = 100)) +
  geom_rect(fill = "gray") +
  coord_cartesian(xlim = c(56885800, 56907500)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_sv +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/DNase_peaks_LCR.pdf", egg::set_panel_size(p, width=unit(4, "cm"), height=unit(0.25, "cm")))

p <- DNase_peaks %>% filter(intensity > 1) %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 50 , ymax = 100)) +
  geom_rect(fill = "gray") +
  coord_cartesian(xlim = c(56840593, 56906000)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_sv +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/DNase_peaks.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(8, "cm"), height=unit(0.5, "cm")))

# stats
condense_peaks <- function(x, variant) {
  variant_name = as.character(variant)
  print(variant_name)
  peak_signal <- tibble()
  search <- c("wildtype", "No sv evidence", variant)
  for(i in 1:nrow(DNase_peaks)) {
    # Note all fibers that were there in the beginning to figure out 0s
    unique_fibers <- filter(x, name %in% search, read_start < DNase_peaks$start[i], read_end > DNase_peaks$end[i]) %>% dplyr::select(fiber, variant, name) %>% 
      mutate(region = paste(DNase_peaks$start[i], DNase_peaks$end[i], sep = "-"), region_name = DNase_peaks$region_name[i]) %>% unique()
    tmp <- filter(x, start > DNase_peaks$start[i], start < DNase_peaks$end[i], name %in% search,
                  read_start < DNase_peaks$start[i], read_end > DNase_peaks$end[i]) %>% group_by(fiber, variant, name) %>% 
      summarise(perc_mod = mean(dplyr::n()/DNase_peaks$AT[i]*100)) %>% 
      mutate(region = paste(DNase_peaks$start[i], DNase_peaks$end[i], sep = "-"), region_name = DNase_peaks$region_name[i]) %>%
      full_join(unique_fibers, by = c("fiber", "variant", "name", "region", "region_name")) %>% 
      mutate(perc_mod = ifelse(is.na(perc_mod), 0, perc_mod))
    peak_signal <- bind_rows(peak_signal, tmp)
  }
  peak_signal <- mutate(peak_signal, name = ifelse(name == "wildtype", "wildtype", variant_name))
  return(peak_signal)
}

# calculate methylation status across peaks
DEL_R4_peaks <- bind_rows(condense_peaks(LCR_variants, "∆R4"), condense_peaks(LCR_variants, "∆R4+iR5")) %>% filter(region %in% c("56902440-56902700", "56903400-56903620"))
DEL_R5_peaks <- bind_rows(condense_peaks(LCR_variants, "∆R5"), condense_peaks(LCR_variants, "iR4+∆R5"), condense_peaks(LCR_variants, "∆R123")) %>% filter(!region %in% c("56902440-56902700", "56903400-56903620"))
DEL_R1234_peaks <- condense_peaks(LCR_variants, "∆R1234") %>% filter(region %in% c("56902440-56902700", "56903400-56903620"))
DEL_R123_peaks <- condense_peaks(LCR_variants, "∆R123")

p <- DEL_R5_peaks %>% filter(!region %in% c("56902440-56902700", "56903400-56903620"), name != "iR4+∆R5") %>%
  ggplot(aes(x = factor(name, levels = c("∆R5", "wildtype", "∆R123")), y = perc_mod, fill = factor(name, levels = c("∆R5", "wildtype", "∆R123")))) +
  geom_boxplot(outlier.size = 0.01, outlier.stroke = 0.1, outlier.color = "gray", size = 0.2, show.legend = F) +
  labs(x = "Locus architecture", y = "% adenines methylated") +
  facet_wrap(~region_name) +
  coord_cartesian(ylim = c(0, 20)) +
  stat_compare_means(comparisons = list(c("∆R5", "wildtype"), c("∆R5", "∆R123")), method = "t.test", label.y = c(14, 16), size = 1.7, bracket.size = 0.2, vjust = -0.2, tip.length = 0.01) +
  scale_fill_manual("Variants", values = c("#a52a2a", "darkgray", "#55AA55")) +
  theme_sv +
  theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/DHS_dR5.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(1.2, "cm"), height=unit(2.5, "cm")))

p <- bind_rows(DEL_R4_peaks, DEL_R1234_peaks, DEL_R123_peaks) %>% filter(region %in% c("56902440-56902700", "56903400-56903620"), name != "∆R4+iR5") %>%
  ggplot(aes(x = factor(name, levels = c("∆R1234", "∆R123", "∆R4", "wildtype")), y = perc_mod, fill = factor(name, levels = c("∆R1234", "∆R123", "∆R4", "wildtype")))) +
  geom_boxplot(outlier.size = 0.01, outlier.stroke = 0.1, size = 0.2, show.legend = F) +
  labs(x = "Locus architecture", y = "% adenines methylated") +
  facet_wrap(~region_name) +
  coord_cartesian(ylim = c(0, 20)) +
  stat_compare_means(comparisons = list(c("∆R4", "wildtype"), c("∆R123", "wildtype"), c("∆R1234", "wildtype")), method = "t.test", label.y = c(12, 14, 16), size = 1.7, bracket.size = 0.2, vjust = -0.2, tip.length = 0.02) +
  scale_fill_manual("Variants", values = c("darkgreen", "#55AA55", "#6495ed", "#808080")) +
  theme_sv +
  theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/DHS_LCR_moving.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(1.7, "cm"), height=unit(2.5, "cm")))



## calculate statistics for CpG
total_reads_LCR <- filter(LCR_variants_mCpG, name %in% c("wildtype", "∆R5", "∆R4")) %>%
  group_by(name, strand) %>%
  summarize(total_reads = n_distinct(fiber))

# Calculate the count of CpG positions for each name at each position
cpg_counts_LCR <- filter(LCR_variants_mCpG, name %in% c("wildtype", "∆R5", "∆R4")) %>%
  group_by(start, name, strand) %>%
  summarize(cpg_count = n()) %>%
  left_join(total_reads_LCR, by = c("name", "strand")) %>%
  group_by(start, name) %>%
  mutate(score = cpg_count / total_reads, keep = ifelse(score == max(score), T, F)) %>%
  group_by(start) %>%
  mutate(total_score = sum(cpg_count, na.rm = T)/sum(total_reads, na.rm = T)) %>%
  ungroup() %>%
  filter(keep == T, cpg_count < total_reads, total_score > 0.05)

total_reads_moving <- filter(LCR_moving_mCpG, name %in% c("wildtype", "∆R123", "∆R1234")) %>%
  group_by(name, strand) %>%
  summarize(total_reads = n_distinct(fiber))

cpg_counts_moving <- filter(LCR_moving_mCpG, name %in% c("wildtype", "∆R123", "∆R1234")) %>%
  group_by(start, name, strand) %>%
  summarize(cpg_count = n()) %>%
  left_join(total_reads_moving, by = c("name", "strand")) %>%
  group_by(start, name) %>%
  mutate(score = cpg_count / total_reads, keep = ifelse(score == max(score), T, F)) %>%
  group_by(start) %>%
  mutate(total_score = sum(cpg_count, na.rm = T)/sum(total_reads, na.rm = T)) %>%
  ungroup() %>%
  filter(keep == T, cpg_count < total_reads, total_score > 0.05)

# dR5
annotate_variant <- function(counts, variant) {
  cpg_count_col <- sym(paste0("cpg_count_", variant))
  total_reads_col <- sym(paste0("total_reads_", variant))
  
  cpg_variant <- counts %>%
    dplyr::select(start, name, cpg_count, total_reads) %>%
    pivot_wider(names_from = name, values_from = c(cpg_count, total_reads), names_sep = "_") %>%
    dplyr::select(start, cpg_count_wildtype, !!cpg_count_col, total_reads_wildtype, !!total_reads_col) %>%
    mutate(total_score = (cpg_count_wildtype + !!cpg_count_col) / 
             (total_reads_wildtype + !!total_reads_col)) %>%
    filter(!is.na(cpg_count_wildtype), !is.na(!!cpg_count_col), total_score > 0.05)
  
  # Function to perform Fisher's exact test
  p_values <- numeric(nrow(cpg_variant))
  
  # Perform Fisher's exact test for each row
  for (i in 1:nrow(cpg_variant)) {
    contingency_table <- matrix(
      c(
        cpg_variant$cpg_count_wildtype[i], cpg_variant$total_reads_wildtype[i] - cpg_variant$cpg_count_wildtype[i],
        cpg_variant[[cpg_count_col]][i], cpg_variant[[total_reads_col]][i] - cpg_variant[[cpg_count_col]][i]
      ),
      nrow = 2
    )
    p_values[i] <- fisher.test(contingency_table)$p.value
  }
  
  cpg_variant$p_value = p_values
  cpg_variant <- mutate(cpg_variant, padj = p.adjust(p_value, method = "BH"), significance = ifelse(padj < 0.05, "significant", "non-significant"),
                        score_wildtype = cpg_count_wildtype/total_reads_wildtype, score_variant = !!cpg_count_col/!!total_reads_col, enrichment = ifelse(score_variant > score_wildtype, variant, "wildtype"), variant = variant)
  
}


cpg_dR4 <- annotate_variant(cpg_counts_LCR, variant = "∆R4")
cpg_dR5 <- annotate_variant(cpg_counts_LCR, variant = "∆R5")
cpg_dR123 <- annotate_variant(cpg_counts_moving, variant = "∆R123")
cpg_dR1234 <- annotate_variant(cpg_counts_moving, variant = "∆R1234")

p <- cpg_counts_LCR %>%
  ggplot(aes(xmin = start, xmax = start+1, ymin = 0, ymax = score, col = name, fill = name)) +
  geom_rect(linewidth = 0.1, show.legend = F) +
  facet_wrap(~name, nrow = 3, scales = "free_y") +
  scale_color_manual("Variant", values = c("#6495ed", "#a52a2a", "#808080")) +
  scale_fill_manual("Variant", values = c("#6495ed", "#a52a2a", "#808080")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = c(0,1)) +
  coord_cartesian(xlim = c(c(56885800, 56907500))) +
  theme_sv
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/OTX2_LCR_mCpG_aggregated.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(4, "cm"), height=unit(0.5, "cm")))


p <- cpg_counts_moving %>%
  ggplot(aes(xmin = start, xmax = start+1, ymin = 0, ymax = score, col = name, fill = name)) +
  geom_rect(linewidth = 0.1, show.legend = F) +
  facet_wrap(~name, nrow = 3, scales = "free_y") +
  scale_color_manual("Variant", values = c("#808080", "#55AA55", "darkgreen")) +
  scale_fill_manual("Variant", values = c("#808080", "#55AA55", "darkgreen")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = c(0,1)) +
  coord_cartesian(xlim = c(56885800, 56907500)) +
  theme_sv
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/OTX2_moving_mCpG_aggregated.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(4, "cm"), height=unit(0.5, "cm")))

p <- ggplot(bind_rows(filter(cpg_dR5, start < 56898376, padj < 0.01), filter(cpg_dR4, start > 56898376, padj < 0.01)), aes(x = start, y = 0.5, col = enrichment)) +
  geom_point(shape = 8, size = 0.5, stroke = 0.3) +
  coord_cartesian(xlim = c(56885800, 56907500)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual("Variant", values = c("#6495ed", "#a52a2a", "#808080")) +
  theme_sv +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/mCpG_LCR_pval.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(4, "cm"), height=unit(0.5, "cm")))


p <- ggplot(bind_rows(filter(cpg_dR1234, start > 56898376, padj < 0.01), filter(cpg_dR123, start > 56885800, padj < 0.01)), aes(x = start, y = 0.5, col = enrichment)) +
  geom_point(shape = 8, size = 0.5, stroke = 0.3) +
  facet_wrap(~variant, nrow = 2) +
  coord_cartesian(xlim = c(56885800, 56907500)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual("Variant", values = c("#55AA55", "darkgreen", "#808080")) +
  theme_sv +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/mCpG_loxP6_pval.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(4, "cm"), height=unit(0.5, "cm")))
