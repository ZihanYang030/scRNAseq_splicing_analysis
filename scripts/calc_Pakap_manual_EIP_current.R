options(stringsAsFactors = FALSE)

personal_dir <- Sys.getenv("PROJECT_DIRECTORY")
if (!nzchar(personal_dir)) stop("Set PROJECT_DIRECTORY")
regions_file <- Sys.getenv("TARGET_REGIONS_FILE", file.path(personal_dir, "config", "target_regions.tsv"))
regions <- read.delim(regions_file, stringsAsFactors = FALSE)
get_region <- function(name) {
  row <- regions[regions$region == name, ]
  if (nrow(row) != 1) stop(paste("Missing or duplicate region:", name))
  row
}

runs <- list(
  # May 1, 2026 - focus on neuron / non-neuron (L23 IT CTX vs oligo)
  # may switch to larger classes later
  list(run_name = "L23ITCTX_matched", class_label = "Glutamatergic", subclass_label = "L2/3 IT CTX"),
  list(run_name = "Oligo_matched", class_label = "Non-neuronal", subclass_label = "Oligo")
)

#out <- num / den
# to avoid den = 0:
safe_divide <- function(num, den) {
  out <- num / den
  out[den == 0] <- NA
  out
}

sum_reads_by_srr <- function(dat, sample_ids) {
  if (nrow(dat) == 0) return(rep(0, length(sample_ids)))
  out <- tapply(dat$unique_reads, dat$SRR, sum)
  vals <- out[sample_ids]
  vals[is.na(vals)] <- 0
  as.numeric(vals)
}

for (cfg in runs) {
  run_name <- cfg$run_name
  result_dir <- file.path(personal_dir, "results", run_name)
  manual_dir <- file.path(result_dir, "manual_Pakap_current_exons")

  junction_file <- file.path(result_dir, "all_Pakap_manual_region_junctions.csv")
  manual_exon_file <- file.path(manual_dir, "Pakap_manual_exons_current_combined.csv")

  out_event_summary <- file.path(result_dir, "Pakap_manual_current_EIP_event_summary.csv")
  out_by_cell <- file.path(result_dir, "Pakap_manual_current_EIP_by_cell.csv")
  out_by_subclass <- file.path(result_dir, "Pakap_manual_current_EIP_by_subclass.csv")
  out_exon_status <- file.path(result_dir, "Pakap_manual_current_exon_status.csv")
  out_exon6_by_cell <- file.path(result_dir, "Pakap_exon6_current_targeted_by_cell.csv")
  out_exon6_summary <- file.path(result_dir, "Pakap_exon6_current_targeted_summary.csv")
  out_exon6_event_def <- file.path(result_dir, "Pakap_exon6_current_event_definition.csv")

  junctions <- read.csv(junction_file, stringsAsFactors = FALSE)
  exons <- read.csv(manual_exon_file, stringsAsFactors = FALSE)

  junctions$SRR <- as.character(junctions$SRR)
  junctions$chrom <- as.character(junctions$chrom)
  junctions$intron_start <- as.integer(junctions$intron_start)
  junctions$intron_end <- as.integer(junctions$intron_end)
  junctions$unique_reads <- as.numeric(junctions$unique_reads)

  exons$chrom <- as.character(exons$chrom)
  exons$start <- as.integer(exons$start)
  exons$end <- as.integer(exons$end)
  exons$gene_name <- as.character(exons$gene_name)
  exons$exon_number <- as.character(exons$exon_number)
  exons$exon_id <- as.character(exons$exon_id)
  exons$manual_structure <- as.character(exons$manual_structure)

  get_junction <- function(chrom, start, end) {
    junctions[junctions$chrom == chrom & junctions$intron_start == start & junctions$intron_end == end, ]
  }

  all_events <- list()
  all_cell_rows <- list()
  all_exon_status <- list()

  for (structure in unique(exons$manual_structure)) {
    structure_exons <- exons[exons$manual_structure == structure, ]
    structure_exons <- structure_exons[order(structure_exons$start, structure_exons$end), ]

    for (i in seq_len(nrow(structure_exons))) {
      exon <- structure_exons[i, ]
      exon_position <- ifelse(i == 1, "first_terminal", ifelse(i == nrow(structure_exons), "last_terminal", "internal"))
      all_exon_status[[length(all_exon_status) + 1]] <- data.frame(
        manual_structure = structure,
        gene_name = exon$gene_name,
        exon_number = exon$exon_number,
        exon_id = exon$exon_id,
        target_start = exon$start,
        target_end = exon$end,
        exon_position = exon_position,
        standard_EIP_calculable = exon_position == "internal",
        stringsAsFactors = FALSE
      )
    }

    if (nrow(structure_exons) < 3) next

    for (i in 2:(nrow(structure_exons) - 1)) {
      upstream <- structure_exons[i - 1, ]
      target <- structure_exons[i, ]
      downstream <- structure_exons[i + 1, ]

      j1_start <- upstream$end + 1
      j1_end <- target$start - 1
      j2_start <- target$end + 1
      j2_end <- downstream$start - 1
      skip_start <- upstream$end + 1
      skip_end <- downstream$start - 1

      inc1 <- if (j1_start <= j1_end) get_junction(target$chrom, j1_start, j1_end) else junctions[FALSE, ]
      inc2 <- if (j2_start <= j2_end) get_junction(target$chrom, j2_start, j2_end) else junctions[FALSE, ]
      exc <- if (skip_start <= skip_end) get_junction(target$chrom, skip_start, skip_end) else junctions[FALSE, ]

      total_inc1 <- ifelse(nrow(inc1) == 0, 0, sum(inc1$unique_reads, na.rm = TRUE))
      total_inc2 <- ifelse(nrow(inc2) == 0, 0, sum(inc2$unique_reads, na.rm = TRUE))
      total_exc <- ifelse(nrow(exc) == 0, 0, sum(exc$unique_reads, na.rm = TRUE))

      all_events[[length(all_events) + 1]] <- data.frame(
        manual_structure = structure, gene_name = target$gene_name,
        exon_number = target$exon_number, exon_id = target$exon_id,
        target_start = target$start, target_end = target$end,
        #target_end = target$end+1,
        upstream_exon_number = upstream$exon_number,
        downstream_exon_number = downstream$exon_number,
        j1_start = j1_start, j1_end = j1_end, j2_start = j2_start, j2_end = j2_end, 
        skip_start = skip_start, skip_end = skip_end,
        total_inc1_unique_reads = total_inc1, total_inc2_unique_reads = total_inc2,
        total_exc_unique_reads = total_exc,
        #I_exon_total = total_inc1, 
        I_exon_total = total_inc1 + total_inc2, E_exon_total = total_exc,
        denominator_total = total_inc1 + total_inc2 + 2 * total_exc,
        EIP_total = safe_divide(total_inc1 + total_inc2, total_inc1 + total_inc2 + 2 * total_exc),
        stringsAsFactors = FALSE
      )

      samples <- unique(c(inc1$SRR, inc2$SRR, exc$SRR))
      samples <- samples[!is.na(samples)]
      if (length(samples) == 0) next

      exon_rows <- data.frame(
        SRR = samples, manual_structure = structure,
        gene_name = target$gene_name, exon_number = target$exon_number,
        exon_id = target$exon_id, target_start = target$start,
        target_end = target$end, stringsAsFactors = FALSE)
      
      exon_rows$I_left <- sum_reads_by_srr(inc1, samples)
      exon_rows$I_right <- sum_reads_by_srr(inc2, samples)
      exon_rows$E_exon <- sum_reads_by_srr(exc, samples)
      exon_rows$I_exon <- exon_rows$I_left + exon_rows$I_right
      exon_rows$denominator <- exon_rows$I_exon + 2*exon_rows$E_exon
      exon_rows$EIP <- safe_divide(exon_rows$I_exon, exon_rows$denominator)

      meta <- unique(junctions[, c("SRR", "sample_name", "class_label", "subclass_label")])
      exon_rows <- merge(exon_rows, meta, by = "SRR", all.x = TRUE)
      all_cell_rows[[length(all_cell_rows) + 1]] <- exon_rows
    }
  }

  event_df <- do.call(rbind, all_events)
  write.csv(event_df, out_event_summary, row.names = FALSE)
  exon_status_df <- do.call(rbind, all_exon_status)
  write.csv(exon_status_df, out_exon_status, row.names = FALSE)

  if (length(all_cell_rows) > 0) {
    cell_df <- do.call(rbind, all_cell_rows)
    write.csv(cell_df, out_by_cell, row.names = FALSE)

    subclass_df <- aggregate(
      cbind(I_left, I_right, I_exon, E_exon, denominator) ~ manual_structure + subclass_label + class_label + gene_name + exon_number + exon_id + target_start + target_end,
      data = cell_df, FUN = sum, na.rm = TRUE)
    
    subclass_df$EIP <- safe_divide(subclass_df$I_exon, subclass_df$denominator)
    write.csv(subclass_df, out_by_subclass, row.names = FALSE)
  }

  chrom <- get_region("exon_5")$chrom
  E5_end <- get_region("exon_5")$end
  E6a_start <- get_region("exon_6a")$start
  E6a_end <- get_region("exon_6a")$end
  E6_terminal_full_start <- get_region("exon_6_terminal")$start
  E6_terminal_full_end <- get_region("exon_6_terminal")$end
  E7_start <- get_region("exon_7")$start

  J_E5_E6a_start <- E5_end + 1
  J_E5_E6a_end <- E6a_start - 1
  J_E6a_E7_start <- E6a_end + 1
  J_E6a_E7_end <- E7_start - 1
  J_E5_E7_start <- E5_end + 1
  #J_E5_E7_end <- E7_start
  J_E5_E7_end <- E7_start - 1

  event_def <- data.frame(
    event = c("J_E5_to_E6a", "J_E6a_to_E7", "J_E5_to_E7_skip_exon6", "terminal_exon6_region"),
    chrom = chrom,
    start = c(J_E5_E6a_start, J_E6a_E7_start, J_E5_E7_start, E6_terminal_full_start),
    end = c(J_E5_E6a_end, J_E6a_E7_end, J_E5_E7_end, E6_terminal_full_end),
    interpretation = c(
      "splice junction from exon 5 into exon 6a; shared by terminal and internal exon 6 transcripts",
      "splice junction from exon 6a to exon 7; proves internal continuation to exon 7",
      "splice junction from exon 5 to exon 7; supports exon 6 skipping",
      "continuous terminal exon 6a/6b/3UTR region; not directly measurable from SJ.out.tab alone"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(event_def, out_exon6_event_def, row.names = FALSE)

  J_E5_E6a <- get_junction(chrom, J_E5_E6a_start, J_E5_E6a_end)
  J_E6a_E7 <- get_junction(chrom, J_E6a_E7_start, J_E6a_E7_end)
  J_E5_E7 <- get_junction(chrom, J_E5_E7_start, J_E5_E7_end)

  samples <- unique(c(J_E5_E6a$SRR, J_E6a_E7$SRR, J_E5_E7$SRR))
  samples <- samples[!is.na(samples)]
  if (length(samples) == 0) next

  exon6_cell <- data.frame(SRR = samples, stringsAsFactors = FALSE)
  exon6_cell$J_E5_E6a <- sum_reads_by_srr(J_E5_E6a, samples)
  exon6_cell$J_E6a_E7 <- sum_reads_by_srr(J_E6a_E7, samples)
  exon6_cell$J_E5_E7_skip <- sum_reads_by_srr(J_E5_E7, samples)

  exon6_cell$I_exon6a_standard <- exon6_cell$J_E5_E6a + exon6_cell$J_E6a_E7
  exon6_cell$E_exon6a_standard <- exon6_cell$J_E5_E7_skip
  exon6_cell$denominator_standard_EIP <- exon6_cell$I_exon6a_standard + 2 * exon6_cell$E_exon6a_standard
  exon6_cell$EIP_exon6a_standard <- safe_divide(exon6_cell$I_exon6a_standard, exon6_cell$denominator_standard_EIP)

  exon6_cell$left_anchor_denominator <- exon6_cell$J_E5_E6a + exon6_cell$J_E5_E7_skip
  exon6_cell$P_left_include_exon6 <- safe_divide(exon6_cell$J_E5_E6a, exon6_cell$left_anchor_denominator)
  exon6_cell$P_left_skip_exon6 <- safe_divide(exon6_cell$J_E5_E7_skip, exon6_cell$left_anchor_denominator)

  exon6_cell$right_internal_denominator <- exon6_cell$J_E6a_E7 + exon6_cell$J_E5_E7_skip
  exon6_cell$P_right_internal_exon6 <- safe_divide(exon6_cell$J_E6a_E7, exon6_cell$right_internal_denominator)
  exon6_cell$P_right_skip_exon6 <- safe_divide(exon6_cell$J_E5_E7_skip, exon6_cell$right_internal_denominator)
  exon6_cell$AltEIP_internal_exon6a_vs_skip <- safe_divide(exon6_cell$J_E6a_E7, exon6_cell$J_E6a_E7 + exon6_cell$J_E5_E7_skip)

  meta <- unique(junctions[, c("SRR", "sample_name", "class_label", "subclass_label")])
  exon6_cell <- merge(exon6_cell, meta, by = "SRR", all.x = TRUE)
  write.csv(exon6_cell, out_exon6_by_cell, row.names = FALSE)

  exon6_summary <- aggregate(
    cbind(J_E5_E6a, J_E6a_E7, J_E5_E7_skip, I_exon6a_standard, E_exon6a_standard, denominator_standard_EIP, left_anchor_denominator, right_internal_denominator) ~ subclass_label + class_label,
    data = exon6_cell, FUN = sum, na.rm = TRUE)

  exon6_summary$EIP_exon6a_standard <- safe_divide(exon6_summary$I_exon6a_standard, exon6_summary$denominator_standard_EIP)
  exon6_summary$P_left_include_exon6 <- safe_divide(exon6_summary$J_E5_E6a, exon6_summary$left_anchor_denominator)
  exon6_summary$P_left_skip_exon6 <- safe_divide(exon6_summary$J_E5_E7_skip, exon6_summary$left_anchor_denominator)
  exon6_summary$P_right_internal_exon6 <- safe_divide(exon6_summary$J_E6a_E7, exon6_summary$right_internal_denominator)
  exon6_summary$P_right_skip_exon6 <- safe_divide(exon6_summary$J_E5_E7_skip, exon6_summary$right_internal_denominator)
  exon6_summary$AltEIP_internal_exon6a_vs_skip <- safe_divide(exon6_summary$J_E6a_E7, exon6_summary$J_E6a_E7 + exon6_summary$J_E5_E7_skip)

  n_cells_exon6 <- aggregate(
    SRR ~ subclass_label + class_label,
    data = exon6_cell[exon6_cell$J_E5_E6a > 0 | exon6_cell$J_E6a_E7 > 0 | exon6_cell$J_E5_E7_skip > 0, ],
    FUN = function(x) length(unique(x))
  )
  colnames(n_cells_exon6)[colnames(n_cells_exon6) == "SRR"] <- "n_cells_with_exon6_junction_evidence"
  exon6_summary <- merge(exon6_summary, n_cells_exon6, by = c("subclass_label", "class_label"), all.x = TRUE)

  write.csv(exon6_summary, out_exon6_summary, row.names = FALSE)
  print(exon6_summary)
}
