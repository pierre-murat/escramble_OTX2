# Read vcf files
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(Repitools)
library(tidyverse)
library(plyranges)
rc <- function(x) {toupper(spgs::reverseComplement(x))}

theme_sv <-   theme_bw(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.25, "cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black"))


# ==== Generating annotations  ====
chr_list <- sprintf("chr%s",c(seq(1,22,1), "X", "Y"))
loxPsym_sites = tibble(precise = c(56840403, 56851953, 56864718, 56884634, 56898376, 56907608)) %>% mutate(chr = "chr14", start = precise - 100, end = precise + 100) %>% GRanges()

sv_annotations_enhancers <- tibble(
  orientation = c("+-c", "++c", "--c", "-+c", "+-t", "++t", "--t", "-+t", "+-l", "++l", "--l"),
  SVTYPE = c("DEL", "INV", "INV", "DUP", "TRA", "TRA", "TRA", "TRA", "INS", "FB", "FB"),
  category = c("Deletion", "Inversion", "Inversion", "Duplication", "Translocation", "Translocation", "Translocation", "Translocation", "Insertion", "Fold back", "Fold back"))

loxP6_variants <- read_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/annotation_variants.txt") %>% mutate(name = str_replace_all(name, "deletion", "∆"))
loxP6_variants_trans <- read_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/annotation_variants_trans.txt") %>% mutate(name = str_replace_all(name, "deletion", "∆"))

variant_types <- tibble(
  name = rep(c("∆R4", "∆R5", "∆R45", "iR4", "iR5", "iR45", "∆R4+iR5", "iR4+∆R5", "wildtype"), 4),
  variant = rep(c("DEL13742", "DEL9232", "DEL22974", "INV13742", "INV9232", "INV22974", "INV22974_INV9232", "INV13742_INV22974", "WT"), 4),
  gate = rep(c("dark", "dim", "medium", "bright"), each = 9)
)
variant_types_rep <- bind_rows(mutate(variant_types, replicate = 1), mutate(variant_types, replicate = 2))
levels <- c("wildtype", "∆R45", "∆R4", "∆R5", "iR4+∆R5", "iR4", "iR5", "iR45", "∆R4+iR5")

# ==== 1 Define functions ====

compute_outcomes <- function(x) {
  variants <- filter(x, SVTYPE != "INS") %>% separate_rows(read_names, sep = ",") %>% group_by(read_names) %>% mutate(distinct_occurences = n())
  simple <- filter(variants, distinct_occurences == 1, SVTYPE == "DEL") %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(start, end, variant) %>% summarise(supp_reads = n(), read_names = paste0(read_names, collapse = "_"))
  complex <- filter(variants, distinct_occurences > 1, SVTYPE != "INS") %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(read_names) %>% summarise(variant = paste0(variant, collapse = "_")) %>% 
    mutate(variant = ifelse(variant %in% c("INV22974_INV22974", "INV13742_INV13742", "INV9232_INV9232", "INV22974_INV9232", "INV13742_INV22974"), variant, "other")) %>% group_by(variant) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  WT =  x %>% separate_rows(read_names, sep = ",") %>% group_by(read_names) %>% mutate(distinct_occurences = n()) %>%
    filter(distinct_occurences == 3, SVTYPE == "INS") %>% mutate(variant = paste0(SVTYPE, start)) %>% group_by(read_names) %>% summarise(variant = paste0(variant, collapse = "_")) %>% group_by(variant) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_")) %>%
    filter(variant == "INS56884584_INS56898326_INS56907558") %>% mutate(variant = "WT")
  variants <- bind_rows(simple, complex, WT) %>% ungroup() %>% mutate(perc_reads = supp_reads/sum(supp_reads)*100)
}

load_vcf_savana = function(path, path_reads) {
  print("Reading vcf file")
  vcf = VariantAnnotation::readVcf(path, "GRCh38")
  print("Processing vcf file")
  gr = GRanges(vcf@rowRanges, MATEID = vcf@info$MATEID, normal_reads = vcf@info$NORMAL_SUPPORT, supp_reads = vcf@info$TUMOUR_SUPPORT, NORMAL_DP = vcf@info$NORMAL_DP, TUMOR_DP = vcf@info$TUMOUR_DP, SVLEN = vcf@info$SVLEN, STRANDS = vcf@info$BP_NOTATION,
               REF = vcf@fixed$REF, ALT = vcf@fixed$ALT, VARIANT_ID = names(vcf)) %>% gintools::unique_granges()
  read_supp <- read_tsv(path_reads, col_types = c("ccc")) %>% dplyr::select(VARIANT_ID, "read_names" = "TUMOUR_SUPPORTING_READS") %>% mutate(read_names = str_replace_all(read_names, ",", "_"))
  
  df <- tibble(annoGR2DF(gr)) %>% mutate(ALT = as.character(flatten(ALT))) %>%
    mutate(ALT = str_remove_all(ALT, "^[ACGT]|\\[|\\]|[ACGT]$"), MATEID = as.character(MATEID), TUMOR_DP = as.character(TUMOR_DP), VARIANT_ID = str_remove(VARIANT_ID, "_1$")) %>%
    separate(ALT, into = c("chr_2", "end"), sep = ":", remove = F) %>%
    separate(MATEID, into = c("ID", "temp", "Mate", sep = "_"), convert = T) %>%
    separate(TUMOR_DP, into = c("depth_start", "depth_end"), sep = ",") %>%
    mutate(depth_start = as.numeric(str_remove_all(depth_start, "c\\(")), depth_end = as.numeric(str_remove_all(depth_start, "\\)")),
           end = as.numeric(str_remove_all(end, "N")), type = ifelse(chr != chr_2, "t", ifelse(abs(end - start) < 50, "l", "c")), orientation = paste0(STRANDS, type)) %>% 
    left_join(sv_annotations_enhancers, by = "orientation") %>% left_join(read_supp, by = "VARIANT_ID") %>%
    filter(normal_reads == 0, supp_reads > 0, SVLEN > 10, !category %in% c("Translocation", "Insertion"), chr == "chr14", start < end, start > 56839500, end < 56909000) %>%
    dplyr::select(chr, start, chr_2, end, SVLEN, SVTYPE, category, orientation, STRANDS, VARIANT_ID, depth_start, depth_end, supp_reads, read_names) %>% distinct()
  
  df_ins <- tibble(annoGR2DF(gr)) %>% filter(chr == "chr14", STRANDS == "<INS>", supp_reads > 5, start > 56839500, end < 56909000, SVLEN > 10, SVLEN < 100) %>% 
    dplyr::select(chr, start, end, supp_reads, SVLEN, VARIANT_ID) %>% mutate(chr_2 = "chr14", SVTYPE = "INS", category = "Insertion", VARIANT_ID = str_remove(VARIANT_ID, "_1$")) %>%
    left_join(read_supp, by = "VARIANT_ID")
  df_all <- bind_rows(df, df_ins)
  return(df_all)
}


call_enhancer_scramble_snf = function(x, insertion_sites, liftover_sites) {
  # find start and end sites that are overlapping with loxPsym site insertions
  print("find start and end sites that are overlapping with loxPsym site insertions")
  start_sites <- x %>% dplyr::select(chr, start) %>% mutate(end = start)
  end_sites <- x %>% dplyr::select("chr" = "chr_2", "start" = "end") %>% mutate(end = start)
  gr_start <- makeGRangesFromDataFrame(start_sites, keep.extra.columns=T)
  gr_end <- makeGRangesFromDataFrame(end_sites, keep.extra.columns=T)
  overlaps_start <- subsetByOverlaps(gr_start, insertion_sites) %>% annoGR2DF() %>% mutate(region = paste0(chr, ":", start)) %>% `$`(region)
  overlaps_end <- subsetByOverlaps(gr_end, insertion_sites) %>% annoGR2DF() %>% mutate(region = paste0(chr, ":", start)) %>% `$`(region)
  
  print("filtering variants")
  passing_rearrangements <- x %>% mutate(region_start = paste0(chr, ":", start), region_end = paste0(chr_2, ":", end)) %>% filter(region_start %in% overlaps_start & region_end %in% overlaps_end, chr %in% chr_list, chr_2 %in% chr_list)
  passing_rearrangements_lifted <- liftover_coordinates_snf(passing_rearrangements, insertion_sites = liftover_sites)
  print(paste(nrow(passing_rearrangements_lifted), "rearrangements harmonized"))
  passing_rearrangements_lifted <- mutate(passing_rearrangements_lifted, SVLEN = end - start)
  return(passing_rearrangements_lifted)
  #passing_rearrangements <- passing_rearrangements %>% left_join(dplyr::select(passing_rearrangements_lifted, id, "start_lift" = "start", "end_lift" = "end"), by = "id") 
  #  mutate(start = ifelse(!is.na(start_lift), start_lift, start), end = ifelse(!is.na(end_lift), end_lift, end)) %>% dplyr::select(-start_lift, -end_lift)
}

liftover_coordinates_snf <- function(x, insertion_sites) {
  print("Harmonizing coordinates") 
  x <- mutate(x, id = 1:nrow(x))
  # liftover start
  sv_liftover_start <- dplyr::select(x, chr, start) %>% mutate(end = start, sv_start = start) %>% GRanges()
  sv_liftover_start <- join_overlap_left(insertion_sites, sv_liftover_start) %>% filter(!is.na(sv_start)) %>% annoGR2DF() %>% distinct()
  sv_liftover_start <- left_join(dplyr::select(sv_liftover_start, chr, start, sv_start), dplyr::select(x, chr,"sv_start" = "start", chr_2, end, SVLEN, supp_reads, SVTYPE, read_names, id), by = c("chr", "sv_start"), relationship = "many-to-many") %>% distinct()
  # liftover end
  sv_liftover_end <- dplyr::select(x, "chr" = "chr_2", "start" = "end") %>% mutate(end = start, sv_end = start) %>% GRanges()
  sv_liftover_end <- join_overlap_left(insertion_sites, sv_liftover_end) %>% filter(!is.na(sv_end)) %>% annoGR2DF() %>% distinct()
  sv_liftover <- left_join(dplyr::select(sv_liftover_end, "chr_2" = "chr", "end" = "start", sv_end), dplyr::select(sv_liftover_start, chr, start, "sv_end" = "end", chr_2, end, supp_reads, SVTYPE, read_names, id), by = c("chr_2", "sv_end"), relationship = "many-to-many") %>% filter(!is.na(supp_reads)) %>% distinct()
  # final cleaning to combine reads that are now the same
  passing_rearrangements <- sv_liftover %>% group_by(chr, start, chr_2, end, SVTYPE) %>% summarise(supp_reads = sum(supp_reads), id = id[1], read_names = paste0(read_names, collapse = "_")) %>% mutate(start = start + 50, end = end + 50) %>% ungroup()
}

compute_outcomes_snf <- function(x) {
  variants <- filter(x, SVTYPE != "INS") %>% separate_rows(read_names, sep = "_") %>% group_by(read_names) %>% mutate(distinct_occurences = dplyr::n())
  simple <- filter(variants, distinct_occurences == 1, SVTYPE %in% c("DEL", "INV")) %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(start, end, variant) %>% summarise(supp_reads = n(), read_names = paste0(read_names, collapse = "_"))
  complex <- filter(variants, distinct_occurences > 1, SVTYPE != "INS") %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(read_names) %>% summarise(variant = paste0(variant, collapse = "_")) %>% 
    mutate(variant = ifelse(variant %in% c("INV22974_INV22974", "INV13742_INV13742", "INV9232_INV9232", "INV22974_INV9232", "INV13742_INV22974"), variant, "other")) %>% group_by(variant) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  WT =  x %>% separate_rows(read_names, sep = "_") %>% group_by(read_names) %>% mutate(distinct_occurences = n()) %>%
    filter(distinct_occurences == 3, SVTYPE == "INS") %>% mutate(variant = paste0(SVTYPE, start)) %>% group_by(read_names) %>% summarise(variant = paste0(variant, collapse = "_")) %>% group_by(variant) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_")) %>%
    filter(variant == "INS56884584_INS56898326_INS56907558") %>% 
    mutate(variant = "WT")
  variants <- bind_rows(simple, complex, WT) %>% ungroup() %>% mutate(perc_reads = supp_reads/sum(supp_reads)*100)
}

compute_outcomes_loxP6 <- function(mini, bwa, max_cov) {
  x <- full_join(mini, bwa, by = c("chr", "start", "chr_2", "end", "SVTYPE", "SVLEN")) %>%
    mutate(supp_reads.x = ifelse(is.na(supp_reads.x), 0, supp_reads.x), supp_reads.y = ifelse(is.na(supp_reads.y), 0, supp_reads.y)) %>%
    mutate(supp_reads = ifelse(supp_reads.x > supp_reads.y, supp_reads.x, supp_reads.y), read_names = ifelse(supp_reads.x > supp_reads.y, read_names.x, read_names.y))
  
  all_sv_reads <- filter(x, SVTYPE != "INS", SVLEN > 0) %>% separate_rows(read_names, sep = "_") %>% `$`(read_names) %>% unique() %>% length()
  variants <- filter(x, SVTYPE != "INS", SVLEN > 0) %>% separate_rows(read_names, sep = "_") %>% group_by(read_names) %>% mutate(distinct_occurences = dplyr::n())
  simple <- filter(variants, distinct_occurences == 1, SVTYPE %in% c("DEL", "INV")) %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(start, end, variant) %>% summarise(supp_reads = n(), read_names = paste0(read_names, collapse = "_"))
  complex <- filter(variants, distinct_occurences > 1) %>% mutate(variant = paste0(SVTYPE, SVLEN)) %>% group_by(read_names) %>% summarise(variant = paste0(variant, collapse = "_")) %>% group_by(variant) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  no_evidence <- tibble(variant = "no_sv_evidence", supp_reads = max_cov - all_sv_reads)
  variants <- complex %>% ungroup() %>% mutate(perc_reads = supp_reads/max_cov*100)
  variants <- bind_rows(simple, complex, no_evidence) %>% mutate(perc_reads = supp_reads/max_cov*100)
}


# 2 ==== LCR screen ====

# LCR screen
LCR_bright_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_bright_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_bright_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_medium_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_medium_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_medium_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_dim_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dim_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dim_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_dark_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dark_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dark_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

LCR_bright_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_bright_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_bright_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_medium_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_medium_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_medium_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_dim_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dim_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dim_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
LCR_dark_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dark_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_LCR_dark_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

# Compute outcomes
LCR_bright_R1_o <- compute_outcomes_snf(LCR_bright_R1) %>% mutate(gate = "bright")
LCR_medium_R1_o <- compute_outcomes_snf(LCR_medium_R1) %>% mutate(gate = "medium")
LCR_dim_R1_o <- compute_outcomes_snf(LCR_dim_R1) %>% mutate(gate = "dim")
LCR_dark_R1_o <- compute_outcomes_snf(LCR_dark_R1) %>% mutate(gate = "dark")

LCR_bright_R2_o <- compute_outcomes_snf(LCR_bright_R2) %>% mutate(gate = "bright")
LCR_medium_R2_o <- compute_outcomes_snf(LCR_medium_R2) %>% mutate(gate = "medium")
LCR_dim_R2_o <- compute_outcomes_snf(LCR_dim_R2) %>% mutate(gate = "dim")
LCR_dark_R2_o <- compute_outcomes_snf(LCR_dark_R2) %>% mutate(gate = "dark")


# Combine
LCR_combined_R1 <- bind_rows(LCR_bright_R1_o, LCR_medium_R1_o, LCR_dim_R1_o, LCR_dark_R1_o) %>% mutate(replicate = 1)
LCR_combined_R2 <- bind_rows(LCR_bright_R2_o, LCR_medium_R2_o, LCR_dim_R2_o, LCR_dark_R2_o) %>% mutate(replicate = 2)
LCR_combined <- bind_rows(LCR_combined_R1, LCR_combined_R2) %>% mutate(SVLEN = end - start) %>% dplyr::select(variant, supp_reads, perc_reads, gate, replicate, read_names, SVLEN) %>% full_join(variant_types_rep, by = c("variant", "gate", "replicate")) %>% filter(!is.na(name)) %>% mutate(perc_reads = ifelse(is.na(perc_reads), min(perc_reads, na.rm = T), perc_reads), supp_reads = ifelse(is.na(supp_reads), 0, supp_reads))


# Combine outcomes
LCR_combined %>% group_by(gate, replicate) %>% summarise(total_reads = sum(supp_reads))

p <- LCR_combined %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = perc_reads, col = factor(name, levels = levels), group = factor(name, levels = levels))) +
  stat_summary(geom = "line", fun = mean) +
  stat_summary(geom = "point", fun = mean) +
  geom_point(size = 0.5) +
  scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ factor(name, levels = levels), scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/lines_LCR.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2, "cm")))

p <- LCR_combined %>% group_by(replicate, name, gate) %>% summarise(mean_perc_reads = mean(perc_reads)) %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = mean_perc_reads, col = factor(name, levels = levels), fill = factor(name, levels = levels))) +
  geom_col(position = "fill") +
  scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  scale_fill_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/plots/event_frequencies_stack.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2, "cm")))

# plot enrichments
Cas9_bright_dim <- LCR_combined %>% mutate(classification = ifelse(gate %in% c("dark", "dim"), "dim", "bright")) %>% group_by(classification, replicate, name) %>% 
  summarise(perc_reads = sum(perc_reads, na.rm = T)) %>%
  pivot_wider(id_cols = c("replicate", "name"), names_from = classification, values_from = perc_reads) %>%
  mutate(log2fc = log2(bright/dim), method = "Cas9") %>% group_by(name) %>% mutate(mean_log2fc = mean(log2fc))
levels <- arrange(Cas9_bright_dim, mean_log2fc) %>% `$`(name) %>% unique()

p <- Cas9_bright_dim %>%
  ggplot(aes(x = factor(name, levels = levels), y = log2fc, col = name)) +
  geom_point(aes(y = mean_log2fc), size = 1.5) +
  geom_point(size = 0.4) +
  geom_segment(aes(x = name, xend = name, y = 0, yend = mean_log2fc), linewidth = 0.3) +
  scale_color_manual("Variant", values = c("#f08080", "#9b72aa", "#a52a2a", "#cd5c5c", "#b0c4de", "#6a4574", "#6495ed", "#4169e1", "#808080")) +
  labs(y = "Log2 fold change\nbright vs dim gates") +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/LCR_enrichments.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(3, "cm"), height=unit(3, "cm")))



# 6 ====  loxP6 screen ====
loxP6_bright_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_bright_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_bright_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_medium_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_medium_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_medium_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_dim_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dim_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dim_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_dark_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dark_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dark_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

loxP6_bright_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_bright_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_bright_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_medium_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_medium_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_medium_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_dim_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dim_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dim_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
loxP6_dark_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dark_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_loxP6_dark_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

# Compute outcomes
loxP6_bright_R1_o <- compute_outcomes_loxP6(loxP6_bright_R1, loxP6_bright_R1, 1147) %>% mutate(gate = "bright")
loxP6_medium_R1_o <- compute_outcomes_loxP6(loxP6_medium_R1, loxP6_medium_R1, 205) %>% mutate(gate = "medium")
loxP6_dim_R1_o <- compute_outcomes_loxP6(loxP6_dim_R1, loxP6_dim_R1, 1309) %>% mutate(gate = "dim")
loxP6_dark_R1_o <- compute_outcomes_loxP6(loxP6_dark_R1, loxP6_dark_R1, 518) %>% mutate(gate = "dark")

loxP6_bright_R2_o <- compute_outcomes_loxP6(loxP6_bright_R2, loxP6_bright_R2, 492) %>% mutate(gate = "bright")
loxP6_medium_R2_o <- compute_outcomes_loxP6(loxP6_medium_R2, loxP6_medium_R2, 714) %>% mutate(gate = "medium")
loxP6_dim_R2_o <- compute_outcomes_loxP6(loxP6_dim_R2, loxP6_dim_R2, 1522) %>% mutate(gate = "dim")
loxP6_dark_R2_o <- compute_outcomes_loxP6(loxP6_dark_R2, loxP6_dark_R2, 1141) %>% mutate(gate = "dark")

# Combine
loxP6_combined_R1 <- bind_rows(loxP6_bright_R1_o, loxP6_medium_R1_o, loxP6_dim_R1_o, loxP6_dark_R1_o) %>% mutate(replicate = 1)
loxP6_combined_R2 <- bind_rows(loxP6_bright_R2_o, loxP6_medium_R2_o, loxP6_dim_R2_o, loxP6_dark_R2_o) %>% mutate(replicate = 2)
loxP6_combined <- bind_rows(loxP6_combined_R1, loxP6_combined_R2) %>% ungroup() %>% mutate(SVLEN = end - start) %>% dplyr::select(variant, supp_reads, perc_reads, gate, replicate, read_names, SVLEN) %>% 
  full_join(loxP6_variants, by = c("variant", "gate", "replicate")) %>% filter(!is.na(name)) %>% mutate(perc_reads = ifelse(is.na(perc_reads), 0.1, perc_reads), supp_reads = ifelse(is.na(supp_reads), 0, supp_reads))

loxP6_bright_dim <- loxP6_combined %>% mutate(classification = ifelse(gate %in% c("dark", "dim"), "dim", "bright")) %>% group_by(classification, replicate, name) %>% 
  summarise(perc_reads = sum(perc_reads, na.rm = T), supp_reads = sum(supp_reads, na.rm = T)) %>%
  pivot_wider(id_cols = c("replicate", "name"), names_from = classification, values_from = c(perc_reads, supp_reads)) %>%
  group_by(name) %>% mutate(supp_reads = supp_reads_bright + supp_reads_dim, supp_reads_rep = sum(supp_reads)) %>% filter(supp_reads_rep > 10) %>%
  mutate(perc_reads_bright = ifelse(is.na(perc_reads_bright), 0.1, perc_reads_bright), perc_reads_dim = ifelse(is.na(perc_reads_dim), 0.1, perc_reads_dim), log2fc = log2(perc_reads_bright/perc_reads_dim), method = "Cas9") %>% group_by(name) %>% mutate(mean_log2fc = mean(log2fc))
levels <- arrange(loxP6_bright_dim, mean_log2fc) %>% `$`(name) %>% unique()


p <- loxP6_combined %>% group_by(name) %>% mutate(total_reads = sum(supp_reads)) %>% filter(total_reads > 15) %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = perc_reads, group = factor(name))) +
  stat_summary(geom = "line", fun = mean, col = "steelblue") +
  stat_summary(geom = "point", fun = mean, col = "steelblue") +
  geom_point(size = 0.5, col = "steelblue") +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ factor(name), scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/lines_loxP6.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2, "cm")))

loxP6_bright_dim %>% group_by(name) %>% mutate(total_reads = sum(supp_reads)) %>% filter(total_reads > 15) %>%
  ggplot(aes(x = factor(name, levels = levels), y = log2fc, col = name)) +
  geom_point(aes(y = mean_log2fc), size = 1.5) +
  geom_point(size = 0.4) +
  geom_segment(aes(x = name, xend = name, y = 0, yend = mean_log2fc), linewidth = 0.3) +
  # scale_color_manual("Variant", values = c("#f08080", "#9b72aa", "#a52a2a", "#cd5c5c", "#b0c4de", "#6a4574", "#6495ed", "#4169e1", "#808080")) +
  labs(y = "Log2 fold change\nbright vs dim gates") +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))


# TRANSFECTIONS
# load files
trans_brilliant_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_brilliant_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_brilliant_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_medium_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_medium_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_medium_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_medium_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_medium_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_medium_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_dim_R1 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dim_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dim_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_dim_R2 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dim_R2_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dim_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_dark_R1 <- load_vcf_savana("/Users/jk24/Desktop/vcf/OTX2_trans_dark_R1_OTX2.sv_breakpoints.vcf", "/Users/jk24/Desktop/vcf/OTX2_trans_dark_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
trans_dark_R2 <- load_vcf_savana("/Users/jk24/Desktop/vcf/OTX2_trans_dark_R2_OTX2.sv_breakpoints.vcf", "/Users/jk24/Desktop/vcf/OTX2_trans_dark_R2_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

trans_dark_R10 <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dark_R1_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_trans_dark_R1_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

# using full reads
trans_brilliant_R1_o <- compute_outcomes_loxP6(trans_brilliant_R1, trans_brilliant_R1, 379) %>% mutate(gate = "brilliant", replicate = 1)
trans_medium_R1_o <- compute_outcomes_loxP6(trans_medium_R1, trans_medium_R1, 2744) %>% mutate(gate = "medium", replicate = 1)
trans_medium_R2_o <- compute_outcomes_loxP6(trans_medium_R2, trans_medium_R2, 1173) %>% mutate(gate = "medium", replicate = 2)
trans_dim_R1_o <- compute_outcomes_loxP6(trans_dim_R1, trans_dim_R1, 1072) %>% mutate(gate = "dim", replicate = 1)
trans_dim_R2_o <- compute_outcomes_loxP6(trans_dim_R2, trans_dim_R2, 591) %>% mutate(gate = "dim", replicate = 2)
trans_dark_R1_o <- compute_outcomes_loxP6(trans_dark_R1, trans_dark_R1, 1377) %>% mutate(gate = "dark", replicate = 1)
trans_dark_R2_o <- compute_outcomes_loxP6(trans_dark_R2, trans_dark_R2, 182) %>% mutate(gate = "dark", replicate = 1)

trans_dark_R10_o <- compute_outcomes_loxP6(trans_dark_R10, trans_dark_R10, 400) %>% mutate(gate = "dark", replicate = 2)


t4_bright <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_bright_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_bright_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
t4_medium <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_medium_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_medium_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
t4_dim <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_dim_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_dim_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)
t4_dark <- load_vcf_savana("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_dark_OTX2.sv_breakpoints.vcf", "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/vcf/OTX2_t4_dark_OTX2.sv_breakpoints_read_support.tsv") %>% call_enhancer_scramble_snf(loxPsym_sites, loxPsym_sites)

t4_bright_o <- compute_outcomes_loxP6(t4_bright, t4_bright, 381) %>% mutate(gate = "bright")
t4_medium_o <- compute_outcomes_loxP6(t4_medium, t4_medium, 446) %>% mutate(gate = "medium")
t4_dim_o <- compute_outcomes_loxP6(t4_dim, t4_dim, 42) %>% mutate(gate = "dim")
t4_dark_o <- compute_outcomes_loxP6(t4_dark, t4_dark, 530) %>% mutate(gate = "dark")


trans_combined <- bind_rows(trans_brilliant_R1_o, trans_medium_R1_o, trans_medium_R2_o, trans_dim_R1_o, trans_dim_R2_o, trans_dark_R1_o, trans_dark_R2_o) %>% ungroup() %>% mutate(SVLEN = end - start) %>% dplyr::select(variant, supp_reads, perc_reads, gate, replicate, read_names, SVLEN) %>% 
  full_join(loxP6_variants_trans, by = c("variant", "gate", "replicate")) %>% filter(!is.na(name)) %>% 
  mutate(perc_reads = ifelse(is.na(perc_reads), 0.1, perc_reads), supp_reads = ifelse(is.na(supp_reads), 0, supp_reads)) %>% dplyr::select(-replicate) %>% distinct() %>%
  group_by(variant) %>% mutate(n = sum(supp_reads), n_zero = sum(perc_reads == 0.1)) %>% ungroup()

t4_combined <- bind_rows(t4_bright_o, t4_medium_o, t4_dark_o) %>% dplyr::select(variant, supp_reads, perc_reads, gate) %>% full_join(loxP6_variants, by = c("variant", "gate")) %>% filter(!is.na(name)) %>% mutate(perc_reads = ifelse(is.na(perc_reads), 0.1, perc_reads), supp_reads = ifelse(is.na(supp_reads), 0, supp_reads)) %>% select(-replicate) %>% distinct() %>%
  group_by(variant) %>% mutate(n = sum(supp_reads))

trans_combined %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "brilliant")), y = perc_reads, col = variant, group = variant)) +
  stat_summary(geom = "line", fun = mean) +
  stat_summary(geom = "point", fun = mean) +
  geom_point(size = 0.5) +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ variant, scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))

p <- trans_combined %>% group_by(name) %>% mutate(total_reads = sum(supp_reads)) %>% ungroup() %>% group_by(name) %>% mutate(n = dplyr::n()) %>% filter(total_reads > 20, n_zero < 3) %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "brilliant")), y = perc_reads, group = factor(name))) +
  stat_summary(geom = "line", fun = mean, col = "steelblue") +
  stat_summary(geom = "point", fun = mean, col = "steelblue") +
  geom_point(size = 0.5, col = "steelblue") +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ factor(name), scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/plots/lines_trans.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2, "cm"), height=unit(2, "cm")))

t4_combined %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = perc_reads, col = variant, group = variant)) +
  stat_summary(geom = "line", fun = mean) +
  stat_summary(geom = "point", fun = mean) +
  geom_point(size = 0.5) +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ variant, scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))


loxP6_combined %>% group_by(name) %>% mutate(total_reads = sum(supp_reads)) %>% filter(total_reads > 15) %>%
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = perc_reads, col = variant, group = variant)) +
  stat_summary(geom = "line", fun = mean) +
  stat_summary(geom = "point", fun = mean) +
  geom_point(size = 0.5) +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ variant, scales = "free", ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))


# Save variant reads
read_names_loxP6 <- bind_rows(filter(trans_combined, gate != "dark"), loxP6_combined) %>% filter(!is.na(read_names)) %>% separate_longer_delim(read_names, delim = "_") %>% dplyr::select("fiber" = "read_names", variant, name, SVLEN, gate)
read_names_LCR <- LCR_combined %>% filter(!is.na(read_names)) %>% separate_longer_delim(read_names, delim = "_") %>% dplyr::select("fiber" = "read_names", variant, name, SVLEN, gate, replicate)
write_tsv(read_names_loxP6, "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/read_names/read_names_loxP6.tsv")
write_tsv(read_names_LCR, "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/read_names/read_names_LCR.tsv")


outcomes_combined <- bind_rows(LCR_combined_R1, LCR_combined_R2)
outcomes_combined %>% filter(variant == "DEL13742") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/DEL13742_all.tsv", col_names = F)
outcomes_combined %>% filter(variant == "DEL9232") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/DEL9232_all.tsv", col_names = F)
outcomes_combined %>% filter(variant == "INV22974_INV9232") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/INV22974_INV9232_all.tsv", col_names = F)
outcomes_combined %>% filter(variant == "INV13742_INV22974") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/INV13742_INV22974_all.tsv", col_names = F)

trans_combined %>% filter(name == "∆R234") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/trans_delR234.tsv", col_names = F)
trans_combined %>% filter(name == "∆R1234") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/trans_delR1234.tsv", col_names = F)
trans_combined %>% filter(name == "∆R2345") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/trans_delR2345.tsv", col_names = F)
trans_combined %>% filter(name == "∆R1234") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/trans_delR1234.tsv", col_names = F)
trans_combined %>% filter(name == "no sv evidence") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/trans_no_sv_ev.tsv", col_names = F)

outcomes_combined %>% filter(variant == "WT") %>% separate_longer_delim(read_names, delim = "_") %>% separate_longer_delim(read_names, delim = ",") %>% dplyr::select(read_names) %>% write_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/WT_all.tsv", col_names = F)

write_tsv(read_names, "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/read_names/read_names.tsv", col_names = F)

# ==== Depth based analysis ====
bg_loxP6_files <- read_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/loxP6_screen/coverage/samples_bedgraphs.txt", col_names = "file")

bg_annotation <- read_tsv("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/LCR_screen/coverage/loxP6_regions.bed.txt", col_names = c("chr", "start", "end", "region"))
read_bedgraphs <- function(files, path) {
  files <- mutate(files, path = paste0(path, file), sample = str_remove(file, "_region_coverage.bedgraph"))
  data = list()
  
  for(i in 1:nrow(files)) {
    print(files$path[i])
    data[[i]] = read_tsv(files$path[i], col_names = c("chr", "start", "end", "coverage")) %>% mutate(sample = files$sample[[i]])
  }
  return(bind_rows(data))
}
bg_loxP6 <- read_bedgraphs(bg_loxP6_files, "/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/loxP6_screen/coverage/") %>% mutate(width = end - start) %>% mutate(cov = coverage/width) %>% left_join(bg_annotation, by = c("chr", "start", "end"))
max_covs <- filter(bg_loxP6, region %in% c("R0", "R6")) %>% group_by(sample) %>% summarise(flanking_cov = max(cov))

bg_loxP6 <- bg_loxP6 %>% filter(!region %in% c("R0", "R6")) %>% left_join(max_covs, by = "sample") %>% 
  separate(sample, into = c("gene", "experiment", "gate", "replicate"), sep = "_") %>% mutate(perc_flanking = cov/flanking_cov*100, replicate = ifelse(is.na(replicate), "R1", replicate))

regions_bright_dim <- bg_loxP6 %>% filter(gate != "dark") %>% mutate(classification = ifelse(gate %in% c("dark", "dim"), "dim", "bright")) %>% 
  group_by(classification, replicate, region, experiment) %>% 
  summarise(perc_flanking = mean(perc_flanking, na.rm = T)) %>%
  pivot_wider(id_cols = c("replicate", "region", "experiment"), names_from = classification, values_from = perc_flanking) %>%
  mutate(log2fc = log2(bright/dim), fc = bright/dim, method = "Cas9") %>% group_by(region, experiment) %>% 
  mutate(mean_log2fc = mean(log2fc), mean_fc = mean(fc))


bg_loxP6 %>% 
  ggplot(aes(x = factor(gate, levels = c("dark", "dim", "medium", "bright")), y = perc_flanking, group = factor(region))) +
  stat_summary(geom = "line", fun = mean, col = "steelblue") +
  stat_summary(geom = "point", fun = mean, col = "steelblue") +
  geom_point(size = 0.5, col = "steelblue") +
  #scale_color_manual("Variant", values = c("#808080", "#a52a2a","#f08080", "#cd5c5c", "#6a4574", "#b0c4de", "#6495ed", "#4169e1","#9b72aa")) +
  labs(y = "Percent of events in sorting gate", x = "Sorting gate") +
  facet_wrap(~ factor(region), ncol = 5) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))

p <- regions_bright_dim %>%
  ggplot(aes(x = factor(region), y = fc)) +
  geom_point(aes(y = mean_fc), size = 1.5, col = "gray30") +
  geom_point(size = 0.4, col = "gray30") +
  geom_segment(aes(x = region, xend = region, y = 0, yend = mean_fc), col = "gray30", linewidth = 0.3) +
  labs(y = "Log2 fold change\nbright vs dim gates") +
  facet_wrap(~experiment) +
  theme_sv +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1))
ggsave("/Users/jonas/Library/CloudStorage/OneDrive-Personal/PhD/OTX2/data/loxP6_screen/plots/region_enrichments.pdf", device = cairo_pdf, egg::set_panel_size(p, width=unit(2.5, "cm"), height=unit(3, "cm")))



# 1 ST
# 2 4
# 3 H
# 4 3
# 5 2
# 6 LT
# 7 HtSide
# 8 LtR
# 9 LtRtL
# 10 LtRtLtH
# 11 All toes

set.seed(9)
x <- tibble(loc = round(runif(11, 1, 11), 0), times = round(runif(11, 1, 4), 0)-1, mult = round(sample(c(rep(1,9), 9, 10), 11), 0))
runif(11, 1, 11)
sample(c(rep(1,9), 9, 10), 11)



