#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
positional_args <- args[!grepl("^--", args)]

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

infer_batch_label <- function(batch_name) {
  sub("_.*$", "", batch_name)
}

parse_num <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", trimws(as.character(x)))))
}

quant <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))

base_dir_default <- if (length(positional_args) >= 1) positional_args[[1]] else getwd()
base_dir <- normalizePath(get_arg("--base_dir", base_dir_default), mustWork = FALSE)
current_batch <- get_arg("--current_batch", basename(dirname(base_dir)))
batch_label <- get_arg("--batch_label", infer_batch_label(current_batch))
out_dir <- normalizePath(get_arg("--out_dir", file.path(base_dir, "boxplot_inputs")), mustWork = FALSE)
feishu_long <- normalizePath(get_arg("--feishu_long", file.path(base_dir, "feishu_abundance_history_long.tsv")), mustWork = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("base_dir: ", base_dir)
message("current_batch: ", current_batch)
message("feishu_long: ", feishu_long)
message("out_dir: ", out_dir)

required_cols <- c("sample", "batch", "taxid", "species", "gc_pct", "expected_pct", "abundance_pct")
if (!file.exists(feishu_long)) stop("Missing Feishu abundance long table: ", feishu_long)

all_abund <- read.delim(feishu_long, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
missing_cols <- setdiff(required_cols, names(all_abund))
if (length(missing_cols) > 0) {
  stop("Feishu abundance long table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

all_abund$sample <- trimws(as.character(all_abund$sample))
all_abund$batch <- trimws(as.character(all_abund$batch))
all_abund$taxid <- trimws(as.character(all_abund$taxid))
all_abund$species <- trimws(as.character(all_abund$species))
all_abund$gc_pct <- parse_num(all_abund$gc_pct)
all_abund$expected_pct <- parse_num(all_abund$expected_pct)
all_abund$abundance_pct <- parse_num(all_abund$abundance_pct)
all_abund$clade_reads <- 0
all_abund$taxon_reads <- 0
all_abund$abundance_source <- "feishu_measured_pct"
all_abund$source_file <- feishu_long

all_abund <- all_abund[
  nzchar(all_abund$sample) &
    nzchar(all_abund$batch) &
    nzchar(all_abund$taxid) &
    nzchar(all_abund$species) &
    is.finite(all_abund$abundance_pct),
]
if (nrow(all_abund) == 0) stop("No valid measured abundance rows parsed from Feishu long table.")

current_targets <- unique(all_abund[all_abund$batch == current_batch, c(
  "sample", "taxid", "species", "gc_pct", "expected_pct"
)])
if (nrow(current_targets) == 0) {
  stop("No current-batch rows found in Feishu long table for current_batch=", current_batch)
}
current_targets$target_source <- basename(feishu_long)
current_targets$target_order <- seq_len(nrow(current_targets))
current_targets$gc_source <- ifelse(is.finite(current_targets$gc_pct), "feishu_sheet", "missing")

species_gc_table <- unique(current_targets[, c("taxid", "species", "gc_pct", "gc_source")])
species_gc_table <- species_gc_table[order(species_gc_table$gc_pct, species_gc_table$species), ]
gc_out <- file.path(out_dir, "species_gc_table.tsv")
write.table(species_gc_table, gc_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote species GC table: ", gc_out)

history_values <- merge(
  unique(current_targets[, c("sample", "taxid", "species", "gc_pct", "gc_source")]),
  all_abund[all_abund$batch != current_batch, c(
    "sample", "batch", "taxid", "abundance_pct", "clade_reads", "taxon_reads", "abundance_source"
  )],
  by = c("sample", "taxid"),
  all.x = FALSE,
  sort = FALSE
)
if (nrow(history_values) == 0) stop("No historical measured abundance rows found for current target species.")

current_values <- all_abund[all_abund$batch == current_batch, c(
  "sample", "taxid", "abundance_pct", "clade_reads", "taxon_reads", "abundance_source", "expected_pct"
)]
names(current_values) <- c(
  "sample", "taxid", "current_abundance_pct", "current_clade_reads",
  "current_taxon_reads", "current_abundance_source", "expected_pct"
)

split_keys <- paste(history_values$sample, history_values$taxid, history_values$species, sep = "||")
summary_list <- lapply(split(history_values, split_keys), function(x) {
  data.frame(
    sample = x$sample[1],
    taxid = x$taxid[1],
    species = x$species[1],
    gc_pct = x$gc_pct[1],
    gc_source = x$gc_source[1],
    n_history = nrow(x),
    history_min = min(x$abundance_pct, na.rm = TRUE),
    history_q1 = quant(x$abundance_pct, 0.25),
    history_median = stats::median(x$abundance_pct, na.rm = TRUE),
    history_q3 = quant(x$abundance_pct, 0.75),
    history_max = max(x$abundance_pct, na.rm = TRUE),
    history_mean = mean(x$abundance_pct, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

summary_tbl <- do.call(rbind, summary_list)
summary_tbl <- merge(summary_tbl, current_values, by = c("sample", "taxid"), all.x = TRUE, sort = FALSE)
summary_tbl$current_abundance_pct[is.na(summary_tbl$current_abundance_pct)] <- 0
summary_tbl$current_clade_reads[is.na(summary_tbl$current_clade_reads)] <- 0
summary_tbl$current_taxon_reads[is.na(summary_tbl$current_taxon_reads)] <- 0
summary_tbl$current_abundance_source[is.na(summary_tbl$current_abundance_source)] <- "missing_current_batch"

summary_tbl$range_status <- ifelse(
  summary_tbl$current_abundance_pct > summary_tbl$history_max,
  "above_history_max",
  ifelse(
    summary_tbl$current_abundance_pct < summary_tbl$history_min,
    "below_history_min",
    "within_history_range"
  )
)
summary_tbl$iq_range_status <- ifelse(
  summary_tbl$current_abundance_pct > summary_tbl$history_q3,
  "above_history_IQR",
  ifelse(
    summary_tbl$current_abundance_pct < summary_tbl$history_q1,
    "below_history_IQR",
    "within_history_IQR"
  )
)

summary_tbl <- summary_tbl[order(summary_tbl$sample, summary_tbl$gc_pct, summary_tbl$species), ]
history_values <- history_values[order(history_values$sample, history_values$gc_pct, history_values$species, history_values$batch), ]
current_targets <- current_targets[order(current_targets$sample, current_targets$target_order), ]

summary_file <- file.path(out_dir, "historical_distribution_summary.tsv")
history_file <- file.path(out_dir, "history_values_for_current_species.tsv")
all_file <- file.path(out_dir, "all_bracken_species_abundance_long.tsv")
targets_file <- file.path(out_dir, "current_targets_for_boxplot.tsv")

write.table(summary_tbl, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(history_values[, c(
  "batch", "sample", "taxid", "species", "gc_pct", "gc_source",
  "abundance_pct", "clade_reads", "taxon_reads", "abundance_source"
)], history_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_abund, all_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(current_targets, targets_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote:")
message(summary_file)
message(history_file)
message(all_file)
message(targets_file)
message(gc_out)
