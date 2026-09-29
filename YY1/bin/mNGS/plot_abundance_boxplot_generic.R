#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
positional_args <- args[!grepl("^--", args)]

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

infer_batch_label <- function(base_dir) {
  batch_name <- basename(dirname(normalizePath(base_dir, mustWork = FALSE)))
  batch_name
}

base_dir_default <- if (length(positional_args) >= 1) positional_args[[1]] else getwd()
base_dir <- normalizePath(get_arg("--base_dir", base_dir_default), mustWork = FALSE)
input_dir <- get_arg("--input_dir", file.path(base_dir, "boxplot_inputs"))
summary_file <- get_arg("--summary", file.path(input_dir, "historical_distribution_summary.tsv"))
history_file <- get_arg("--history", file.path(input_dir, "history_values_for_current_species.tsv"))
out_dir <- get_arg("--out_dir", input_dir)
batch_label <- get_arg("--batch_label", infer_batch_label(base_dir))
sample_filter <- get_arg("--sample", "")
prefix <- get_arg("--prefix", "boxplot_abundance_history")
y_gc_order <- get_arg("--y_gc_order", "bottom_to_top")

if (!file.exists(summary_file)) stop("Missing summary file: ", summary_file)
if (nzchar(history_file) && !file.exists(history_file)) {
  message("History file not found; plotting boxplots without historical batch points: ", history_file)
  history_file <- ""
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("base_dir: ", base_dir)
message("input_dir: ", input_dir)
message("summary: ", summary_file)
message("history: ", ifelse(nzchar(history_file), history_file, "(none)"))
message("out_dir: ", out_dir)
message("batch_label: ", batch_label)

required_summary_cols <- c(
  "sample", "taxid", "species", "gc_pct",
  "history_min", "history_q1", "history_median", "history_q3", "history_max",
  "current_abundance_pct", "range_status"
)
summary_tbl <- read.delim(summary_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
missing_summary <- setdiff(required_summary_cols, names(summary_tbl))
if (length(missing_summary) > 0) {
  stop("Summary file is missing required columns: ", paste(missing_summary, collapse = ", "))
}

num_summary_cols <- c(
  "gc_pct", "history_min", "history_q1", "history_median", "history_q3",
  "history_max", "current_abundance_pct"
)
for (cc in num_summary_cols) summary_tbl[[cc]] <- suppressWarnings(as.numeric(summary_tbl[[cc]]))

if (nzchar(history_file)) {
  history_tbl <- read.delim(history_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_history_cols <- c("sample", "taxid", "species", "abundance_pct")
  missing_history <- setdiff(required_history_cols, names(history_tbl))
  if (length(missing_history) > 0) {
    stop("History file is missing required columns: ", paste(missing_history, collapse = ", "))
  }
  history_tbl$abundance_pct <- suppressWarnings(as.numeric(history_tbl$abundance_pct))
} else {
  history_tbl <- data.frame(sample = character(), taxid = character(), species = character(), abundance_pct = numeric())
}

status_col <- c(
  within_history_range = "#008837",
  above_history_max = "#d7191c",
  below_history_min = "#2c7bb6",
  no_history = "#6b7280"
)
status_label <- c(
  within_history_range = "Current within historical min-max",
  above_history_max = "Current above historical max",
  below_history_min = "Current below historical min",
  no_history = "No historical comparison"
)
status_short <- c(
  within_history_range = "within range",
  above_history_max = "above max",
  below_history_min = "below min",
  no_history = "no history"
)
current_point_col <- "#d7191c"

safe_file <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

draw_one <- function(sid) {
  sm <- summary_tbl[summary_tbl$sample == sid, ]
  sm <- sm[order(ifelse(is.na(sm$gc_pct), Inf, sm$gc_pct), sm$species), ]
  hs <- history_tbl[history_tbl$sample == sid, ]

  species_key <- paste(sm$taxid, sm$species, sep = "||")
  labels <- ifelse(
    is.na(sm$gc_pct),
    sm$species,
    sprintf("%s  | GC %.1f%%", sm$species, sm$gc_pct)
  )
  y_pos <- seq_along(species_key)
  if (identical(y_gc_order, "top_to_bottom")) {
    y_pos <- rev(y_pos)
  }
  names(y_pos) <- species_key
  hs$key <- paste(hs$taxid, hs$species, sep = "||")
  hs$y <- y_pos[hs$key]
  sm$key <- paste(sm$taxid, sm$species, sep = "||")
  sm$y <- y_pos[sm$key]
  hs <- hs[!is.na(hs$y), ]

  n <- nrow(sm)
  xmax <- max(c(hs$abundance_pct, sm$current_abundance_pct, sm$history_max), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  xlim <- c(0, xmax * 1.24)
  height <- max(1800, 120 * n + 580)
  png_file <- file.path(out_dir, paste0(prefix, "_", safe_file(sid), ".png"))

  png(png_file, width = 3800, height = height, res = 180)
  par(mar = c(6, 23.5, 5.5, 14.5), xpd = NA, family = "sans")
  plot(
    NA,
    xlim = xlim,
    ylim = c(0.4, n + 0.85),
    yaxt = "n",
    xlab = "Species abundance (%)",
    ylab = "",
    main = paste0(batch_label, " current abundance over historical boxplot\n", sid),
    cex.main = 1.05
  )
  axis(2, at = sm$y, labels = labels, las = 1, cex.axis = 0.78, tick = FALSE)
  old_xpd <- par("xpd")
  par(xpd = FALSE)
  abline(v = axTicks(1), col = "#eeeeee", lty = 1)
  par(xpd = old_xpd)

  box_h <- 0.18
  whisk_h <- 0.15
  for (i in seq_len(nrow(sm))) {
    row <- sm[i, ]
    yy <- row$y
    vals <- hs$abundance_pct[hs$key == row$key]
    if (length(vals) > 0) {
      points(vals, rep(yy, length(vals)), pch = 16, col = adjustcolor("#374151", alpha.f = 0.24), cex = 0.75)
    }
    segments(row$history_min, yy, row$history_q1, yy, col = "#9ca3af", lwd = 1.5)
    segments(row$history_q3, yy, row$history_max, yy, col = "#9ca3af", lwd = 1.5)
    segments(row$history_min, yy - whisk_h, row$history_min, yy + whisk_h, col = "#9ca3af", lwd = 1.5)
    segments(row$history_max, yy - whisk_h, row$history_max, yy + whisk_h, col = "#9ca3af", lwd = 1.5)
    rect(row$history_q1, yy - box_h, row$history_q3, yy + box_h, col = "#d1d5db", border = "#6b7280", lwd = 0.9)
    segments(row$history_median, yy - box_h, row$history_median, yy + box_h, col = "#111827", lwd = 2.2)
  }

  for (i in seq_len(nrow(sm))) {
    row <- sm[i, ]
    st <- row$range_status
    if (!st %in% names(status_col)) st <- "no_history"
    points(row$current_abundance_pct, row$y, pch = 16, col = current_point_col, lwd = 1.0, cex = 0.85)
    text(
      row$current_abundance_pct,
      row$y + 0.34,
      labels = sprintf("%s %.2f%% %s", batch_label, row$current_abundance_pct, status_short[[st]]),
      cex = 0.64,
      col = status_col[[st]],
      pos = 4,
      offset = 0.35
    )
  }

  legend_x <- xlim[2] + diff(xlim) * 0.035
  legend_y <- n + 0.85
  legend(
    legend_x,
    legend_y,
    legend = c("Historical box: Q1-Q3", "Historical whisker: min-max", "Historical median", "Historical batch point", paste0(batch_label, " current point"), status_label),
    fill = c("#d1d5db", NA, NA, NA, NA, rep(NA, length(status_col))),
    border = c("#6b7280", NA, NA, NA, NA, rep(NA, length(status_col))),
    col = c(NA, "#9ca3af", "#111827", adjustcolor("#374151", alpha.f = 0.45), current_point_col, status_col),
    lwd = c(NA, 1.5, 2.2, NA, NA, rep(NA, length(status_col))),
    pch = c(NA, NA, NA, 16, 16, rep(NA, length(status_col))),
    pt.bg = c(NA, NA, NA, NA, NA, rep(NA, length(status_col))),
    pt.cex = c(NA, NA, NA, 0.9, 1.0, rep(NA, length(status_col))),
    text.col = c(rep("#111827", 5), status_col),
    bty = "n",
    cex = 0.78,
    ncol = 1,
    xjust = 0,
    yjust = 1
  )

  gc_sentence <- if (identical(y_gc_order, "top_to_bottom")) {
    "Species ordered by GC content from low to high, top to bottom."
  } else {
    "Species ordered by GC content from low to high, bottom to top."
  }
  mtext(
    paste("Read guide: gray box = historical Q1-Q3; whiskers = historical min-max; black tick = median; small gray dots = historical batches; red dot = current batch.", gc_sentence),
    side = 1,
    line = 4.6,
    cex = 0.72,
    col = "#374151"
  )
  dev.off()
  png_file
}

samples <- unique(summary_tbl$sample)
if (nzchar(sample_filter)) samples <- intersect(samples, sample_filter)
if (length(samples) == 0) stop("No samples to plot.")

files <- vapply(samples, draw_one, character(1))
writeLines(files)
