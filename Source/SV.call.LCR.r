
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(VariantAnnotation)
library(Repitools)
library(stringr)
library(plyranges)

# Function for formatting and analysing structural variants at recombination sites

###############
# OTX2 SV caller

# Define loxPsym sites coordinates (+- 500 nt)

loxPsym.LCR.coordinates <- c(56884634, 56898376, 56907608)
loxPsym.LCR.sorted.coordinates <- loxPsym.LCR.coordinates[order(loxPsym.LCR.coordinates)]
loxPsym.LCR.sites.gr <- tibble(precise = loxPsym.LCR.sorted.coordinates) %>% mutate(chr = "chr14", start = precise - 500, end = precise + 500, sites = c("A", "B", "C")) %>% GRanges()

# Compute distances from all loxPsym sites
loxPsym.LCR.sites.dist.df <- expand_grid(sites = c("A", "B", "C"), sites.2 = c("A", "B", "C")) %>% left_join(as.data.frame(loxPsym.LCR.sites.gr) %>%
                                                                                       dplyr::select(sites, pos.1 = precise)) %>% left_join(as.data.frame(loxPsym.LCR.sites.gr) %>%
                                                                                                                                              dplyr::select(sites.2 = sites, pos.2 = precise)) %>% mutate(dist = pos.1 - pos.2) %>% filter(dist > 0) %>% arrange(-dist)

################
# savana

# Format savana vcf output

OTX2.SV.call.savana <- function(path.vcf, path.tsv) {
  # Load vcf
  vcf <- VariantAnnotation::readVcf(path.vcf, "GRCh38")
  # Format gr
  gr <- GRanges(vcf@rowRanges, SVLEN = vcf@info$SVLEN, SVTYPE = vcf@info$SVTYPE, BP_NOTATION = vcf@info$BP_NOTATION, REF = vcf@fixed$REF, ALT = unlist(vcf@fixed$ALT), FILTER = vcf@fixed$FILTER, normal_support = vcf@info$NORMAL_SUPPORT, tumor_support = vcf@info$TUMOUR_SUPPORT)
  VARIANT_ID <- names(gr)
  VARIANT_ID <- gsub("_1", "", VARIANT_ID)
  VARIANT_ID <- gsub("_2", "", VARIANT_ID)
  gr$VARIANT_ID <- VARIANT_ID
  # Filter gr
  gr <- gr[seqnames(gr) == "chr14"]
  # Filter df for tumor SVs
  df <- as.data.frame(gr) %>% filter(tumor_support >= 1)
  # Annotate SVs
  df <- df %>% mutate(SVTYPE = case_when(BP_NOTATION == "--" ~ "INV",
                                         BP_NOTATION == "-+" ~ "DUP",
                                         BP_NOTATION == "+-" ~ "DEL",
                                         BP_NOTATION == "++" ~ "INV",
                                         T ~ "INS"))
  # Add read information
  tsv <- read.table(path.tsv, header = T) %>% dplyr::select(VARIANT_ID, TUMOUR_SUPPORTING_READS)
  df <- df %>% left_join(tsv, by = "VARIANT_ID") %>% drop_na() %>% dplyr::rename(read_names = TUMOUR_SUPPORTING_READS)
  # Extract SV ends
  df <- df %>% dplyr::select(-end) %>% 
    mutate(ALT = str_remove_all(ALT, "\\["), ALT = str_remove_all(ALT, "\\]"), ALT = str_remove_all(ALT, "N")) %>%
    separate(ALT, into = c("chr_2", "end"), sep = ":", remove = F) %>%
    mutate(end = as.numeric(gsub("[^0-9.-]", "", end))) %>% 
    mutate(chr_2 = as.character(gsub("[^0-9.-]", "", chr_2))) %>% 
    mutate(chr_2 = paste0("chr", chr_2)) %>% filter(chr_2 == "chr14") %>% dplyr::rename(chr = seqnames)
  # reformat df
  out.df <- tibble(df) %>% mutate(read_names = gsub(",", "_", read_names)) %>% 
    mutate(region_start = paste(chr, start, sep = ':'), region_end = paste(chr_2, end, sep = ':')) %>%
    dplyr::select(chr, start, chr_2, end, SVLEN, SVTYPE, supp_reads = tumor_support, read_names, region_start, region_end) %>% distinct()
  # Identify SV sites
  SV.start.df <- out.df %>% dplyr::select(seqnames = chr, start = start, end = start, region_start)
  SV.start.gr <- makeGRangesFromDataFrame(SV.start.df, keep.extra.columns = T)
  SV.start.loxPsym.LCR.dist <- distanceToNearest(SV.start.gr, loxPsym.LCR.sites.gr)
  SV.start.df <- cbind.data.frame(region_start = SV.start.gr$region_start, loxPsym_start = loxPsym.LCR.sites.gr[SV.start.loxPsym.LCR.dist@to]$sites, loxPsym_start_dist = SV.start.loxPsym.LCR.dist@elementMetadata$distance) %>% 
    mutate(loxPsym_start = case_when(loxPsym_start_dist != 0 ~ NA, T ~ loxPsym_start)) %>% dplyr::select(-loxPsym_start_dist) %>% distinct()
  SV.end.df <- out.df %>% dplyr::select(seqnames = chr_2, start = end, end = end, region_end)
  SV.end.gr <- makeGRangesFromDataFrame(SV.end.df, keep.extra.columns = T)
  SV.end.loxPsym.LCR.dist <- distanceToNearest(SV.end.gr, loxPsym.LCR.sites.gr)
  SV.end.df <- cbind.data.frame(region_end = SV.end.gr$region_end, loxPsym_end = loxPsym.LCR.sites.gr[SV.end.loxPsym.LCR.dist@to]$sites, loxPsym_end_dist = SV.end.loxPsym.LCR.dist@elementMetadata$distance) %>% 
    mutate(loxPsym_end = case_when(loxPsym_end_dist != 0 ~ NA, T ~ loxPsym_end)) %>% dplyr::select(-loxPsym_end_dist) %>% distinct()
  # Add information to SV dataframe
  out.loxPsym.LCR.df <- out.df %>% left_join(SV.start.df, by = "region_start") %>% left_join(SV.end.df, by = "region_end") %>% drop_na()  
  # Extract information for individual reads
  read.SV.df <- out.loxPsym.LCR.df %>% separate_rows(read_names, sep = "_")
  # Compute number of SV per read
  SV.occurence.df <- read.SV.df %>% group_by(read_names) %>% summarise(SV.count = n())
  # Add information
  read.SV.df <- read.SV.df %>% left_join(SV.occurence.df, by = "read_names")
  # Generate genotype info
  genotype.df <- read.SV.df %>% mutate(genotype = paste0(SVTYPE, "_", loxPsym_start, "-", loxPsym_end))
  # Filter short event based on loxPsym sites
  genotype.df <- genotype.df %>% filter(loxPsym_start < loxPsym_end)
  return(genotype.df)
}

################
# Call SV from savana output

# 1 break point

OTX2.complex.1.BP.SV.call.savana <- function(df) { # df is the output of the OTX2.SV.call.savana function
  complex.1.BP.df <- df %>% filter(SV.count == 2) %>% group_by(read_names) %>% summarise(complex.genotype = paste0(genotype, collapse = "-")) %>%
    group_by(complex.genotype) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  return(complex.1.BP.df)
}

# 2 break points

OTX2.complex.2.BP.SV.call.savana <- function(df) { # df is the output of the OTX2.SV.call.savana function
  complex.2.BP.df <- df %>% filter(SV.count == 4) %>% group_by(read_names) %>% summarise(complex.genotype = paste0(genotype, collapse = "-")) %>%
    group_by(complex.genotype) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  return(complex.2.BP.df)
}

# 3 break points

OTX2.complex.3.BP.SV.call.savana <- function(df) { # df is the output of the OTX2.SV.call.savana function
  complex.3.BP.df <- df %>% filter(SV.count == 6) %>% group_by(read_names) %>% summarise(complex.genotype = paste0(genotype, collapse = "-")) %>%
    group_by(complex.genotype) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  return(complex.3.BP.df)
}

# 4 break points

OTX2.complex.4.BP.SV.call.savana <- function(df) { # df is the output of the OTX2.SV.call.savana function
  complex.4.BP.df <- df %>% filter(SV.count == 8) %>% group_by(read_names) %>% summarise(complex.genotype = paste0(genotype, collapse = "-")) %>%
    group_by(complex.genotype) %>% summarise(supp_reads = dplyr::n(), read_names = paste0(read_names, collapse = "_"))
  return(complex.4.BP.df)
}

# Other functions

unique_granges <- function(sites, sum.cols = FALSE, 
                           rm.cols = NULL, rm.dup.cols = NULL){
  # Checks and balance
  if( !class(sites) == "GRanges" ){
    stop("\n  Sites object is not a GRanges class.\n")
  }
  if( any(sum.cols != FALSE) ){
    if( all(sum.cols == TRUE) ){
      counts_col <- "counts"
    }else{
      counts_col <- sum.cols
    }
    if( any(!counts_col %in% names(GenomicRanges::mcols(sites))) ){
      missing_cols <- counts_col[
        which(!counts_col %in% names(GenomicRanges::mcols(sites)))
      ]
      stop(
        "\n  Could not find column name(s) in sites object: ",
        paste(missing_cols, collapse = ", "), 
        "\n"
      )
    }
  }
  if( !is.null(rm.cols) ){
    if( any(!rm.cols %in% names(GenomicRanges::mcols(sites))) ){
      missing_cols <- rm.cols[
        which(!rm.cols %in% names(GenomicRanges::mcols(sites)))
      ]
      message(
        "Could not find column name(s) to remove in sites object: ",
        paste(missing_cols, collapse = ", ")
      )
      rm.cols <- rm.cols[rm.cols %in% missing_cols]
    }
    if( length(rm.cols) > 0 ){
      GenomicRanges::mcols(sites)[, rm.cols] <- NULL
    }
  }
  if( !is.null(rm.dup.cols) ){
    col_names <- names(GenomicRanges::mcols(sites))
    if( is.character(rm.dup.cols) ){
      col_names <- stringr::str_extract(col_names, rm.dup.cols)
    }
    dup_cols <- names(table(col_names))[table(col_names) > 1]
    dup_idxs <- lapply(dup_cols, function(x) which(col_names == x))
    rm_col_idx <- unlist(lapply(dup_idxs, function(x){
      col_digests <- unlist(lapply(x, function(y){
        digest::digest(GenomicRanges::mcols(sites)[,y, drop = TRUE])
      }))
      x[which(duplicated(col_digests))]
    }))
    GenomicRanges::mcols(sites)[,rm_col_idx] <- NULL
    if( is.character(rm.dup.cols) ){
      names(GenomicRanges::mcols(sites)) <- stringr::str_extract(
        names(GenomicRanges::mcols(sites)), rm.dup.cols
      )
    }
  }
  # Convert sites to a data.frame and remove duplicates
  if( !length(names(sites)) == length(unique(names(sites))) ){
    message("Dropping rownames for data.frame conversion.")
    df <- GenomicRanges::as.data.frame(sites, row.names = NULL)
  }else{
    df <- GenomicRanges::as.data.frame(sites)
  }
  cols <- names(df)
  if( any(sum.cols != FALSE) ){
    counts_pos <- match(counts_col, cols)
  }
  # Sum counts if needed
  if( !any(sum.cols != FALSE) ){
    df <- dplyr::distinct(df)
  }else{
    groups <- lapply(cols[-counts_pos], as.symbol)
    df <- dplyr::group_by(df, .dots = groups) %>%
      dplyr::summarise_all(.funs = "sum") %>%
      dplyr::ungroup()
    df <- df[,cols]
  }
  # Rebuild GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = df$seqnames,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand,
    seqinfo = GenomicRanges::seqinfo(sites)
  )
  GenomicRanges::mcols(gr) <- dplyr::select(df, 6:length(df))
  gr
}



