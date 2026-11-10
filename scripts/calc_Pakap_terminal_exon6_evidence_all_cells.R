options(stringsAsFactors = FALSE)

personal_dir <- Sys.getenv("PROJECT_DIRECTORY")
if (!nzchar(personal_dir)) stop("Set PROJECT_DIRECTORY")

coverage_file <- file.path(personal_dir, "results", "Pakap_exon6a_6b_3UTR_coverage_by_cell_all.tsv")

run_info <- data.frame(
  run_name = c("L23ITCTX_matched", "Oligo_matched"),
  list_file = c(file.path(personal_dir, "L23ITCTX_Acc_List.txt"), file.path(personal_dir, "Oligo_Acc_List.txt")),
  default_subclass_label = c("L2/3 IT CTX", "Oligo"), default_class_label = c("Glutamatergic", "Non-neuronal"),
  stringsAsFactors = FALSE) # low oligo count - switch to other or more non-neuron cell types later?


#out <- num / den
# to avoid den = 0:
safe_divide <- function(num, den) {
  out <- num / den
  out[den == 0] <- NA
  out
}

coverage <- read.delim(coverage_file, stringsAsFactors = FALSE)
coverage$SRR <- as.character(coverage$SRR)

all_cells_list <- list()
for (i in seq_len(nrow(run_info))) {
  cfg <- run_info[i, ]
  srrs <- trimws(readLines(cfg$list_file))
  srrs <- srrs[nzchar(srrs)]
  all_cells_list[[cfg$run_name]] <- data.frame(
    SRR = srrs, run_name = cfg$run_name, subclass_label = cfg$default_subclass_label,
    class_label = cfg$default_class_label, stringsAsFactors = FALSE)
}
all_cells <- do.call(rbind, all_cells_list)

junction_list <- list()
for (run_name in run_info$run_name) {
  f <- file.path(personal_dir, "results", run_name, "Pakap_exon6_current_targeted_by_cell.csv")
  dat <- read.csv(f, stringsAsFactors = FALSE)
  dat$run_name <- run_name
  junction_list[[run_name]] <- dat
}
junction <- do.call(rbind, junction_list)

junction_keep <- junction[, c(
  "SRR", "run_name", "J_E5_E6a", "J_E6a_E7", "J_E5_E7_skip",
  "I_exon6a_standard", "E_exon6a_standard", "denominator_standard_EIP",
  "EIP_exon6a_standard", "P_left_include_exon6", "P_left_skip_exon6",
  "P_right_internal_exon6", "P_right_skip_exon6", "AltEIP_internal_exon6a_vs_skip"
)]

merged <- merge(all_cells, junction_keep, by = c("SRR", "run_name"), all.x = TRUE)
merged <- merge(merged, coverage, by = "SRR", all.x = TRUE)

junction_zero_cols <- c("J_E5_E6a", "J_E6a_E7", "J_E5_E7_skip", "I_exon6a_standard", "E_exon6a_standard", "denominator_standard_EIP")
for (cc in junction_zero_cols) merged[[cc]][is.na(merged[[cc]])] <- 0

coverage_cols <- c("reads_E6a", "reads_E6b", "reads_3UTR", "reads_E7")
for (cc in coverage_cols) merged[[cc]][is.na(merged[[cc]])] <- 0

merged$EIP_exon6a_standard_recalc <- safe_divide(merged$I_exon6a_standard, merged$denominator_standard_EIP)
merged$has_E5_to_E6a_junction <- merged$J_E5_E6a > 0
merged$has_E6a_to_E7_junction <- merged$J_E6a_E7 > 0
merged$has_E5_to_E7_skip_junction <- merged$J_E5_E7_skip > 0
merged$has_E6a_coverage <- merged$reads_E6a > 0
merged$has_E6b_coverage <- merged$reads_E6b > 0
merged$has_3UTR_coverage <- merged$reads_3UTR > 0
merged$has_E7_coverage <- merged$reads_E7 > 0

merged$has_any_exon6_decision_junction <- (merged$has_E5_to_E6a_junction | merged$has_E6a_to_E7_junction | merged$has_E5_to_E7_skip_junction)

merged$terminal_exon6_strong_evidence <- (
  merged$has_E5_to_E6a_junction &
    merged$has_E6b_coverage &
    merged$has_3UTR_coverage &
    !merged$has_E6a_to_E7_junction
)

merged$terminal_exon6_moderate_evidence <- (
  merged$has_E5_to_E6a_junction &
    (merged$has_E6b_coverage | merged$has_3UTR_coverage) &
    !merged$has_E6a_to_E7_junction &
    !merged$terminal_exon6_strong_evidence
)

merged$terminal_exon6_weak_evidence <- (
  merged$has_E5_to_E6a_junction & !merged$has_E6b_coverage & !merged$has_3UTR_coverage & !merged$has_E6a_to_E7_junction)

merged$internal_exon6_to_E7_evidence <- merged$has_E6a_to_E7_junction
merged$skip_exon6_evidence <- merged$has_E5_to_E7_skip_junction

merged$terminal_evidence_category <- "no_terminal_evidence"
merged$terminal_evidence_category[merged$terminal_exon6_weak_evidence] <- "weak_E5_to_E6a_only_no_E6a_to_E7"
merged$terminal_evidence_category[merged$terminal_exon6_moderate_evidence] <- "moderate_E5_to_E6a_plus_E6b_or_3UTR_no_E6a_to_E7"
merged$terminal_evidence_category[merged$terminal_exon6_strong_evidence] <- "strong_E5_to_E6a_plus_E6b_and_3UTR_no_E6a_to_E7"

out_by_cell <- file.path(personal_dir, "results", "Pakap_terminal_exon6_evidence_by_cell_all_cells.csv")
write.csv(merged, out_by_cell, row.names = FALSE)

summary_df <- aggregate(
  cbind(
    has_any_exon6_decision_junction,
    #has_any_exon6_decision_junction_test,
    has_E5_to_E6a_junction, has_E6a_to_E7_junction, has_E5_to_E7_skip_junction,
    has_E6a_coverage, has_E6b_coverage, has_3UTR_coverage, has_E7_coverage,
    terminal_exon6_strong_evidence, terminal_exon6_moderate_evidence,
    terminal_exon6_weak_evidence, internal_exon6_to_E7_evidence,
    skip_exon6_evidence
  ) ~ run_name + subclass_label + class_label,
  data = merged, FUN = sum, na.rm = TRUE)

n_total <- aggregate(SRR ~ run_name + subclass_label + class_label, data = merged, FUN = function(x) length(unique(x)))
colnames(n_total)[colnames(n_total) == "SRR"] <- "n_total_cells"
summary_df <- merge(summary_df, n_total, by = c("run_name", "subclass_label", "class_label"), all.x = TRUE)

summary_df$prop_strong_terminal_out_of_all_cells <- summary_df$terminal_exon6_strong_evidence / summary_df$n_total_cells
summary_df$prop_moderate_terminal_out_of_all_cells <- summary_df$terminal_exon6_moderate_evidence / summary_df$n_total_cells
summary_df$prop_weak_terminal_out_of_all_cells <- summary_df$terminal_exon6_weak_evidence / summary_df$n_total_cells
summary_df$prop_any_terminal_candidate_out_of_all_cells <- (
  summary_df$terminal_exon6_strong_evidence +
    summary_df$terminal_exon6_moderate_evidence +
    summary_df$terminal_exon6_weak_evidence
) / summary_df$n_total_cells

summary_df$prop_strong_terminal_out_of_exon6_junction_cells <- safe_divide(summary_df$terminal_exon6_strong_evidence, summary_df$has_any_exon6_decision_junction)
summary_df$prop_any_terminal_candidate_out_of_exon6_junction_cells <- safe_divide(
  summary_df$terminal_exon6_strong_evidence + summary_df$terminal_exon6_moderate_evidence + summary_df$terminal_exon6_weak_evidence,
  summary_df$has_any_exon6_decision_junction
)
summary_df$prop_internal_E6a_to_E7_out_of_all_cells <- summary_df$internal_exon6_to_E7_evidence / summary_df$n_total_cells
summary_df$prop_skip_exon6_out_of_all_cells <- summary_df$skip_exon6_evidence / summary_df$n_total_cells

out_summary <- file.path(personal_dir, "results", "Pakap_terminal_exon6_evidence_summary_all_cells.csv")
write.csv(summary_df, out_summary, row.names = FALSE)

category_table <- as.data.frame(table(merged$run_name, merged$subclass_label, merged$terminal_evidence_category, useNA = "ifany"))
colnames(category_table) <- c("run_name", "subclass_label", "terminal_evidence_category", "n_cells")
out_category <- file.path(personal_dir, "results", "Pakap_terminal_exon6_evidence_category_counts_all_cells.csv")
write.csv(category_table, out_category, row.names = FALSE)

cat("cell-level table saved in:\n", out_by_cell, "\n")
cat("Saved summary table:\n", out_summary, "\n")
cat("Saved category table:\n", out_category, "\n")
print(summary_df)
