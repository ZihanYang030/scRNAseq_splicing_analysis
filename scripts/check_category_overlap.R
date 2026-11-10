# script to print category overlap
options(stringsAsFactors = FALSE)

personal_dir <- Sys.getenv("PROJECT_DIRECTORY")
if (!nzchar(personal_dir)) stop("Set PROJECT_DIRECTORY")
cell_file <- file.path(personal_dir, "results", "Pakap_terminal_exon6_evidence_by_cell_all_cells.csv")
df <- read.csv(cell_file, stringsAsFactors = FALSE)

df$any_terminal <- df$terminal_exon6_strong_evidence | df$terminal_exon6_moderate_evidence | df$terminal_exon6_weak_evidence
df$any_exon6_decision <- df$has_any_exon6_decision_junction

for (run in unique(df$run_name)) {
  cat("Run:", run, "\n")

  d <- df[df$run_name == run & df$any_exon6_decision, ]

  cat("Cells with any exon-6 decision junction:", nrow(d), "\n\n")
  cat("Any terminal evidence:", sum(d$any_terminal), "\n")
  cat("strong terminal:", sum(d$terminal_exon6_strong_evidence), "\n")
  cat("Moderate terminal:", sum(d$terminal_exon6_moderate_evidence), "\n")
  cat("Weak terminal:", sum(d$terminal_exon6_weak_evidence), "\n")
  cat("internal E6a->E7:", sum(d$has_E6a_to_E7_junction), "\n")
  cat("skipped E5->E7:", sum(d$has_E5_to_E7_skip_junction), "\n\n")

  cat("Overlap: terminal vs skip\n")
  print(table(terminal = d$any_terminal, skip = d$has_E5_to_E7_skip_junction))

  cat("\nOverlap: terminal vs internal\n")
  print(table(terminal = d$any_terminal, internal = d$has_E6a_to_E7_junction))

  cat("\nOverlap: skip vs internal\n")
  print(table(skip = d$has_E5_to_E7_skip_junction, internal = d$has_E6a_to_E7_junction))
}
