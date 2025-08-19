#!/usr/bin/env Rscript

# =============================================================================
# VarClock Simple Robust Plot Generator with Variant Grouping
# =============================================================================
# This script creates reliable plots that properly handle all support types
# and maintain the 72-hour time window while plotting all time points correctly.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(lubridate)
  library(scales)
  library(glue)
  library(paletteer)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input TSV file with variant grouping", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "robust_plots",
              help = "Output directory (default: robust_plots)", metavar = "DIR"),
  make_option(c("-t", "--timing"), type = "character", default = NULL,
              help = "Optional timing file for start time", metavar = "TIMING"),
  make_option(c("--min-reads"), type = "integer", default = 1,
              help = "Minimum reads for plot generation (default: 1)", metavar = "N"),
  make_option(c("--time-window"), type = "integer", default = 72,
              help = "Analysis time window in hours (default: 72)", metavar = "HOURS"),
  make_option(c("--group-multiallelic"), action = "store_true", default = FALSE,
              help = "Group multi-allelic variants together")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: input file is required.", call. = FALSE)
}

input_file <- opt$input
output_dir <- opt$output
timing_file <- opt$timing
min_reads <- opt$`min-reads`
time_window_hours <- opt$`time-window`
group_multiallelic <- opt$`group-multiallelic`

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("VarClock Simple Robust Plot Generator\n")
cat("=====================================\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Minimum reads:", min_reads, "\n")
cat("Time window:", time_window_hours, "hours\n")
cat("Group multi-allelic:", group_multiallelic, "\n\n")

# =============================================================================
# DATA LOADING AND PROCESSING
# =============================================================================
cat("Loading and processing data...\n")

# Load data
df <- read_tsv(input_file, col_types = cols(.default = "c"))
df$num_alts <- as.numeric(df$num_alts)
df$timestamp <- ymd_hms(df$timestamp)

cat("Loaded", nrow(df), "records\n")

# Check required columns
required_cols <- c("variant_group_id", "variant_summary", "allele_match", "timestamp")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Check for gene name column (region_name)
has_gene_names <- "region_name" %in% colnames(df)
if (!has_gene_names) {
  cat("Warning: region_name column not found. Plot titles will not include gene names.\n")
}

# Classify support types more simply and robustly
df <- df %>%
  mutate(
    support_type = case_when(
      str_starts(allele_match, "VARIANT") ~ "variant",
      allele_match == "REFERENCE" ~ "reference",
      str_starts(allele_match, "OTHER") ~ "other",
      TRUE ~ "unknown"
    ),
    # For plotting, combine reference and other as "reference"
    plot_support = case_when(
      str_starts(allele_match, "VARIANT") ~ "variant",
      TRUE ~ "reference"
    ),
    # Use variant group ID for grouping
    plot_group = variant_group_id
  ) %>%
  filter(!is.na(timestamp))

cat("Support type distribution:\n")
print(table(df$support_type))
cat("Plot support distribution:\n")
print(table(df$plot_support))

# =============================================================================
# TIME RANGE SETUP
# =============================================================================
cat("\nSetting up time range...\n")

# Determine start time
if (!is.null(timing_file) && file.exists(timing_file)) {
  cat("Reading timing file:", timing_file, "\n")
  timing_lines <- readLines(timing_file)
  earliest_line <- timing_lines %>%
    str_subset("^Earliest\\s*:") %>%
    str_remove("Earliest\\s*:\\s*")

  if (length(earliest_line) > 0) {
    global_time_min <- ymd_hms(earliest_line)
    cat("Using timing file start time:", as.character(global_time_min), "\n")
  } else {
    global_time_min <- min(df$timestamp, na.rm = TRUE)
    cat("Timing file format issue, using data minimum:", as.character(global_time_min), "\n")
  }
} else {
  global_time_min <- min(df$timestamp, na.rm = TRUE)
  cat("Using data minimum start time:", as.character(global_time_min), "\n")
}

global_time_max <- global_time_min + hours(time_window_hours)

cat("Analysis window:", as.character(global_time_min), "to", as.character(global_time_max), "\n")
cat("Data range:", as.character(min(df$timestamp)), "to", as.character(max(df$timestamp)), "\n")

# =============================================================================
# GROUP FILTERING
# =============================================================================
cat("\nFiltering groups by minimum reads...\n")

group_counts <- df %>%
  group_by(plot_group, variant_summary) %>%
  summarise(
    total_reads = n(),
    variant_reads = sum(plot_support == "variant"),
    reference_reads = sum(plot_support == "reference"),
    region_id = if(has_gene_names) {
      raw_name <- first(region_name)
      # Extract region ID (full region name), excluding coordinate-like patterns
      if (!is.na(raw_name) && raw_name != "" && !str_detect(raw_name, "^chr\\d+:")) {
        raw_name
      } else {
        NA_character_
      }
    } else {
      NA_character_
    },
    variant_type = {
      raw_types <- first(variant_type)
      if (!is.na(raw_types) && raw_types != "") {
        # Split by comma and get unique types
        unique_types <- unique(str_split(raw_types, ",")[[1]])
        paste(unique_types, collapse = ",")
      } else {
        raw_types
      }
    },
    gene_name = if(has_gene_names) {
      raw_name <- first(region_name)
      # Extract gene name from format: 0018_chr1_FUBP1 (take third part)
      if (!is.na(raw_name) && raw_name != "" && !str_detect(raw_name, "^chr\\d+:")) {
        gene_parts <- str_split(raw_name, "_")[[1]]
        if (length(gene_parts) >= 3 && gene_parts[3] != "") {
          gene_parts[3]  # Take third part as gene name (FUBP1)
        } else if (length(gene_parts) > 0 && gene_parts[1] != "") {
          gene_parts[1]  # Fallback to first part for other formats
        } else {
          NA_character_
        }
      } else {
        NA_character_
      }
    } else {
      NA_character_
    },
    .groups = "drop"
  ) %>%
  filter(total_reads >= min_reads) %>%
  arrange(desc(total_reads))

cat("Groups meeting minimum read threshold:", nrow(group_counts), "\n")
print(head(group_counts))

valid_groups <- group_counts$plot_group
df_filtered <- df %>% filter(plot_group %in% valid_groups)

cat("Filtered to", nrow(df_filtered), "records for plotting\n")

# =============================================================================
# TIME SERIES CREATION
# =============================================================================
cat("\nCreating time series...\n")

create_time_series <- function(group_data, group_id, group_summary, region_id = NA, gene_name = NA, variant_type = NA) {

  cat("Processing group:", group_id, "\n")
  cat("  Total reads:", nrow(group_data), "\n")

  # Filter group data to the analysis time window first
  group_data_filtered <- group_data %>%
    filter(timestamp >= global_time_min & timestamp <= global_time_max)

  cat("  Reads in time window:", nrow(group_data_filtered), "\n")

  # Create time sequence with 10-minute intervals for the full 72-hour window
  time_seq <- seq(from = global_time_min, to = global_time_max, by = "10 mins")

  # Round timestamps to 10-minute intervals
  group_data_filtered$timestamp_rounded <- floor_date(group_data_filtered$timestamp, "10 mins")

  # Count reads at each time point for each support type
  time_counts <- group_data_filtered %>%
    count(timestamp_rounded, plot_support, name = "reads_at_time")

  cat("  Time counts summary:\n")
  print(time_counts %>% group_by(plot_support) %>% summarise(total = sum(reads_at_time), .groups = "drop"))

  # Create complete time grid with rounded timestamps to match data
  time_seq_rounded <- floor_date(time_seq, "10 mins")
  support_types <- c("variant", "reference")
  time_grid <- expand_grid(
    timestamp_rounded = unique(time_seq_rounded),
    plot_support = support_types
  )

  # Merge and fill missing values
  complete_data <- time_grid %>%
    left_join(time_counts, by = c("timestamp_rounded", "plot_support")) %>%
    mutate(
      reads_at_time = replace_na(reads_at_time, 0),
      timestamp = timestamp_rounded  # Use rounded timestamp for consistency
    ) %>%
    arrange(plot_support, timestamp)

  cat("  Complete data summary:\n")
  print(complete_data %>% group_by(plot_support) %>% summarise(total = sum(reads_at_time), .groups = "drop"))

  # Calculate cumulative sums
  cumulative_data <- complete_data %>%
    group_by(plot_support) %>%
    mutate(
      cumulative_reads = cumsum(reads_at_time),
      rel_time_hours = as.numeric(difftime(timestamp, global_time_min, units = "hours"))
    ) %>%
    ungroup()

  # Calculate final counts and VAF
  final_variant <- max(cumulative_data[cumulative_data$plot_support == "variant", "cumulative_reads"], na.rm = TRUE)
  final_reference <- max(cumulative_data[cumulative_data$plot_support == "reference", "cumulative_reads"], na.rm = TRUE)

  # Handle edge cases
  if (is.infinite(final_variant)) final_variant <- 0
  if (is.infinite(final_reference)) final_reference <- 0

  total_reads <- final_variant + final_reference
  vaf <- if (total_reads > 0) final_variant / total_reads else 0

  cat("  Final counts - Variant:", final_variant, "Reference:", final_reference, "VAF:", round(vaf, 3), "\n")

  # Create plot
  p <- ggplot(cumulative_data, aes(x = rel_time_hours, y = cumulative_reads, color = plot_support)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    scale_color_manual(
      values = c("variant" = paletteer_d("ggsci::category10_d3")[2],
                 "reference" = paletteer_d("ggsci::category10_d3")[1]),
      labels = c("variant" = "Variant Support", "reference" = "Reference Support"),
      name = "Support Type"
    ) +
    scale_x_continuous(
      name = "Time (hours of sequencing)",
      breaks = seq(0, time_window_hours, by = 12),
      minor_breaks = seq(0, time_window_hours, by = 6),
      limits = c(0, time_window_hours)
    ) +
    scale_y_continuous(
      name = "Cumulative Read Count",
      labels = comma_format()
    ) +
    labs(
      title = {
        clean_summary <- str_replace_all(group_summary, "\\d+-allelic\\s+", "")

        # Add variant type after position if available
        if (!is.na(variant_type) && variant_type != "") {
          # Split clean_summary to insert variant type after position
          parts <- str_split(clean_summary, " ", n = 2)[[1]]
          if (length(parts) == 2) {
            clean_summary <- paste(parts[1], paste("-", variant_type, "-"), parts[2])
          } else {
            clean_summary <- paste(clean_summary, paste0("-", variant_type, "-"))
          }
        }

        if (!is.na(region_id) && region_id != "" && !is.na(gene_name) && gene_name != "") {
          # Both region ID and gene name available: Extract numeric portion - GENEID - variant_info
          numeric_id <- str_extract(region_id, "^\\d+")
          if (!is.na(numeric_id)) {
            title_text <- paste(numeric_id, "-", gene_name, "-", clean_summary)
          } else {
            title_text <- paste(region_id, "-", gene_name, "-", clean_summary)
          }
          str_wrap(title_text, 80)
        } else if (!is.na(gene_name) && gene_name != "" && !str_detect(gene_name, "^chr\\d+:")) {
          # Only gene name available
          str_wrap(paste(gene_name, "-", clean_summary), 75)
        } else {
          # No gene information
          str_wrap(clean_summary, 60)
        }
      },
      subtitle = glue("VAF: {percent(vaf, accuracy = 0.1)} | Total: {comma(total_reads)} | Variant: {comma(final_variant)} | Reference: {comma(final_reference)}")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray50"),
      legend.position = "bottom",
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.3),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5)
    )

  # Add detection thresholds
  thresholds <- c(1, 5, 8, 10, 20, 50)
  palette_colors <- paletteer_d("ggsci::category10_d3")
  colors <- palette_colors[3:8]  # Use colors 3-8 from the palette

  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    if (final_variant >= thresh) {
      # Find first time point where threshold was reached
      first_detection <- cumulative_data %>%
        filter(plot_support == "variant", cumulative_reads >= thresh) %>%
        slice_min(rel_time_hours, n = 1)

      if (nrow(first_detection) > 0) {
        # Calculate precise time elapsed to the minute
        time_hours <- first_detection$rel_time_hours
        hours_part <- floor(time_hours)
        minutes_part <- round((time_hours - hours_part) * 60)

        # Handle cases where rounding gives 60 minutes
        if (minutes_part == 60) {
          hours_part <- hours_part + 1
          minutes_part <- 0
        }

        time_label <- if (hours_part > 0) {
          sprintf("%d reads - %dh %dm", thresh, hours_part, minutes_part)
        } else {
          sprintf("%d reads - %dm", thresh, minutes_part)
        }

        p <- p +
          geom_vline(xintercept = first_detection$rel_time_hours,
                     color = colors[i], linetype = "dashed", alpha = 0.7, size = 0.8) +
          annotate("text", x = first_detection$rel_time_hours,
                   y = max(cumulative_data$cumulative_reads) * 0.95,
                   label = time_label, angle = 90, vjust = -0.5,
                   size = 3, color = colors[i])
      }
    }
  }

  return(list(plot = p, data = cumulative_data, stats = list(
    total_reads = total_reads,
    variant_reads = final_variant,
    reference_reads = final_reference,
    vaf = vaf
  )))
}

# =============================================================================
# GENERATE PLOTS
# =============================================================================
cat("\nGenerating plots...\n")

plot_results <- list()
summary_stats <- data.frame()

for (i in 1:nrow(group_counts)) {
  group_id <- group_counts$plot_group[i]
  group_summary <- group_counts$variant_summary[i]
  region_id <- group_counts$region_id[i]
  gene_name <- group_counts$gene_name[i]
  variant_type <- group_counts$variant_type[i]

  group_data <- df_filtered %>% filter(plot_group == group_id)

  tryCatch({
    result <- create_time_series(group_data, group_id, group_summary, region_id, gene_name, variant_type)

    # Save plot
    safe_filename <- str_replace_all(group_id, "[^A-Za-z0-9_.-]", "_")
    plot_file <- file.path(output_dir, paste0(safe_filename, "_timeseries.png"))

    ggsave(plot_file, result$plot, width = 12, height = 8, dpi = 300)
    cat("Saved plot:", plot_file, "\n")

    # Store results
    plot_results[[group_id]] <- result

    # Add to summary
    summary_stats <- rbind(summary_stats, data.frame(
      group_id = group_id,
      group_summary = str_replace_all(group_summary, "\\d+-allelic\\s+", ""),
      region_id = if(!is.na(region_id)) region_id else "",
      gene_name = if(!is.na(gene_name)) gene_name else "",
      variant_type = if(!is.na(variant_type)) variant_type else "",
      total_reads = result$stats$total_reads,
      variant_reads = result$stats$variant_reads,
      reference_reads = result$stats$reference_reads,
      vaf = result$stats$vaf,
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    cat("ERROR creating plot for", group_id, ":", e$message, "\n")
  })
}

# =============================================================================
# SAVE SUMMARY
# =============================================================================
cat("\nSaving summary statistics...\n")

summary_file <- file.path(output_dir, "plot_summary.tsv")
write_tsv(summary_stats, summary_file)

cat("Summary saved to:", summary_file, "\n")

# =============================================================================
# COMPLETION REPORT
# =============================================================================
cat("\n=== Plot Generation Complete ===\n")
cat("Results:\n")
cat("  - Plots generated:", length(plot_results), "\n")
cat("  - Output directory:", output_dir, "\n")
cat("  - Summary file:", summary_file, "\n")

if (nrow(summary_stats) > 0) {
  cat("\nTop variant groups by total reads:\n")
  print(summary_stats %>% arrange(desc(total_reads)) %>% head())
}

cat("\nConfiguration used:\n")
cat("  - Time window:", time_window_hours, "hours\n")
cat("  - Minimum reads:", min_reads, "\n")
cat("  - Group multi-allelic:", group_multiallelic, "\n")

cat("\nAll plots saved with proper 72-hour time windows and complete time series!\n")
