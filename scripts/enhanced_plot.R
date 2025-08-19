#!/usr/bin/env Rscript

# =============================================================================
# VARCLOCK Enhanced Plot Generator with Multi-Allelic Grouping Support
# =============================================================================
# This script generates cumulative read count plots for variant analysis,
# with enhanced support for multi-allelic variant grouping. It uses the new
# variant_group_id and variant_summary columns to properly group and display
# complex variant patterns.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)    # Command line argument parsing
  library(marquee)     # Enhanced text rendering in plots
  library(lubridate)   # Date/time manipulation
  library(patchwork)   # Plot composition
  library(glue)        # String interpolation
  library(magrittr)    # Pipe operators
  library(dplyr)       # Data manipulation
  library(ggplot2)     # Plotting
  library(readr)       # Reading data files
  library(stringr)     # String manipulation
  library(purrr)       # Functional programming tools
  library(tidyr)       # Data tidying
  library(scales)      # Plot scaling utilities
  library(viridis)     # Color palettes
})

# =============================================================================
# COMMAND LINE ARGUMENT SETUP
# =============================================================================
cat("VarClock Enhanced Plot Generator with Multi-Allelic Grouping Support\n")
cat("====================================================================\n")

# Define options
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to input TSV file (with variant grouping columns)", metavar = "FILE"),
  make_option(c("-t", "--timing"), type = "character", default = NULL,
              help = "Optional file with Earliest/Latest timestamps", metavar = "TIMING"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for plots (default: input_enhanced_plots/)", metavar = "DIR"),
  make_option(c("--group-multiallelic"), action = "store_true", default = FALSE,
              help = "Group multi-allelic variants together in plots"),
  make_option(c("--separate-alleles"), action = "store_true", default = FALSE,
              help = "Create separate plots for each alternate allele"),
  make_option(c("--time-window"), type = "integer", default = 72,
              help = "Analysis time window in hours (default: 72)", metavar = "HOURS"),
  make_option(c("--min-reads"), type = "integer", default = 5,
              help = "Minimum reads required for plot generation (default: 5)", metavar = "N")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: input file is required.", call. = FALSE)
}

varclock_results <- opt$input
timing_file <- opt$timing
output_dir <- opt$output
group_multiallelic <- opt$`group-multiallelic`
separate_alleles <- opt$`separate-alleles`
time_window_hours <- opt$`time-window`
min_reads <- opt$`min-reads`

# Set default output directory
if (is.null(output_dir)) {
  input_basename_noext <- sub("\\.[^.]*$", "", basename(varclock_results))
  output_dir <- paste0(input_basename_noext, "_enhanced_plots")
}

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

cat("Configuration:\n")
cat("  Input file:", varclock_results, "\n")
cat("  Timing file:", ifelse(is.null(timing_file), "None", timing_file), "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Group multi-allelic:", group_multiallelic, "\n")
cat("  Separate alleles:", separate_alleles, "\n")
cat("  Time window:", time_window_hours, "hours\n")
cat("  Minimum reads:", min_reads, "\n")

# =============================================================================
# DATA LOADING AND VALIDATION
# =============================================================================
cat("\n=== Loading and validating data ===\n")

# Load the main varclock results file
cat("Reading TSV file:", varclock_results, "\n")
df <- readr::read_tsv(varclock_results, col_types = readr::cols(.default = "c"))
cat("Initial data dimensions:", nrow(df), "rows x", ncol(df), "columns\n")

# Check for required columns
required_cols <- c("read_id", "timestamp", "allele_match", "variant_group_id",
                  "variant_summary", "region", "region_name", "num_alts")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("This script requires VarClock output with enhanced variant grouping.\n")
  cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
  stop("Missing required columns for enhanced plotting", call. = FALSE)
}

cat("All required columns found!\n")
cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")

# =============================================================================
# DATA PREPROCESSING WITH GROUPING SUPPORT
# =============================================================================
cat("\n=== Processing data with variant grouping ===\n")

# Transform and clean the data
df <- df %>%
  mutate(
    # Parse timestamp and round to nearest second
    timestamp_raw = ymd_hms(timestamp, tz = "UTC"),
    timestamp = round_date(timestamp_raw, unit = "secs"),
    # Parse numeric fields
    num_alts = as.numeric(num_alts),
    # Determine if this read contains a variant
    contains_variant = str_starts(allele_match, "VARIANT"),
    # Extract specific alternate allele information
    alt_allele_num = case_when(
      str_starts(allele_match, "VARIANT:ALT") ~ str_extract(allele_match, "ALT\\d+"),
      TRUE ~ NA_character_
    ),
    # Classify support type with more detail
    support_type = case_when(
      allele_match == "REFERENCE" ~ "reference",
      str_starts(allele_match, "VARIANT:ALT") ~ "variant",
      str_starts(allele_match, "OTHER:") ~ "other",
      TRUE ~ "unknown"
    ),
    # Create simplified support classification that includes 'other'
    support = case_when(
      str_starts(allele_match, "VARIANT") ~ "variant",
      allele_match == "REFERENCE" ~ "reference",
      str_starts(allele_match, "OTHER:") ~ "other",
      TRUE ~ "unknown"
    ),
    # Add multi-allelic flag
    is_multiallelic = num_alts > 1,
    # Create enhanced plot grouping variable
    plot_group = case_when(
      group_multiallelic ~ variant_group_id,
      separate_alleles & !is.na(alt_allele_num) ~ paste0(variant_group_id, ":", alt_allele_num),
      TRUE ~ variant_group_id
    ),
    # Create plot title from variant summary and region
    plot_title = case_when(
      group_multiallelic ~ paste0(variant_summary, " in ", region_name),
      separate_alleles & !is.na(alt_allele_num) ~ paste0(variant_summary, " (", alt_allele_num, ") in ", region_name),
      TRUE ~ paste0(variant_summary, " in ", region_name)
    )
  ) %>%
  # Remove rows with invalid timestamps
  filter(!is.na(timestamp_raw))

cat("After preprocessing:\n")
cat("  - Records with valid timestamps:", nrow(df), "\n")
cat("  - Timestamp range:", as.character(min(df$timestamp, na.rm = TRUE)),
    "to", as.character(max(df$timestamp, na.rm = TRUE)), "\n")
cat("  - Unique variant groups:", length(unique(df$variant_group_id)), "\n")
cat("  - Unique plot groups:", length(unique(df$plot_group)), "\n")
cat("  - Multi-allelic variants:", sum(df$is_multiallelic, na.rm = TRUE), "records\n")

# Support type distribution
cat("  - Support type distribution:\n")
support_summary <- df %>% count(support_type, sort = TRUE)
print(support_summary)

# Multi-allelic summary
cat("  - Multi-allelic breakdown:\n")
multiallelic_summary <- df %>%
  filter(is_multiallelic) %>%
  count(variant_summary, sort = TRUE) %>%
  head(10)
print(multiallelic_summary)

# =============================================================================
# TIME RANGE DETERMINATION
# =============================================================================
cat("\n=== Determining analysis time range ===\n")

# Determine global time minimum from timing file or data
if (!is.null(timing_file) && file.exists(timing_file)) {
  cat("Reading timing file:", timing_file, "\n")
  timing_lines <- readLines(timing_file)

  earliest_line <- timing_lines %>%
    str_subset("^Earliest\\s*:") %>%
    str_remove("Earliest\\s*:\\s*")

  if (length(earliest_line) > 0) {
    global_time_min_raw <- ymd_hms(earliest_line, tz = "UTC")
    global_time_min <- floor_date(global_time_min_raw, unit = "secs")
    cat("Using timing file for start time:", as.character(global_time_min), "\n")
  } else {
    cat("Warning: Could not find 'Earliest:' line in timing file, using data minimum\n")
    global_time_min <- floor_date(min(df$timestamp, na.rm = TRUE), unit = "secs")
  }
} else {
  cat("Using data minimum for start time\n")
  global_time_min <- floor_date(min(df$timestamp, na.rm = TRUE), unit = "secs")
}

# Set global time maximum based on time window
global_time_max <- global_time_min + hours(time_window_hours)

cat("Analysis time window:\n")
cat("  - Start:", as.character(global_time_min), "\n")
cat("  - End:", as.character(global_time_max), "\n")
cat("  - Duration:", time_window_hours, "hours\n")

# =============================================================================
# PLOT GROUP FILTERING AND VALIDATION
# =============================================================================
cat("\n=== Filtering plot groups ===\n")

# Count reads per plot group
group_counts <- df %>%
  group_by(plot_group, plot_title) %>%
  summarise(
    total_reads = n(),
    variant_reads = sum(support == "variant"),
    reference_reads = sum(support == "reference"),
    other_reads = sum(support_type == "other"),
    unique_alt_alleles = length(unique(alt_allele_num[!is.na(alt_allele_num)])),
    is_multiallelic = any(is_multiallelic),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads))

cat("Plot group summary (top 10):\n")
print(head(group_counts, 10))

# Filter groups with sufficient reads
valid_groups <- group_counts %>%
  filter(total_reads >= min_reads) %>%
  pull(plot_group)

cat("\nFiltering results:\n")
cat("  - Total plot groups:", nrow(group_counts), "\n")
cat("  - Groups with >=", min_reads, "reads:", length(valid_groups), "\n")
cat("  - Groups to be plotted:", length(valid_groups), "\n")

if (length(valid_groups) == 0) {
  stop("No plot groups meet the minimum read threshold (", min_reads, " reads)", call. = FALSE)
}

# Filter data to valid groups
df_filtered <- df %>%
  filter(plot_group %in% valid_groups)

cat("Filtered data:", nrow(df_filtered), "records\n")

# =============================================================================
# TIME SERIES PREPARATION
# =============================================================================
cat("\n=== Preparing time series data ===\n")

# Generate time grid for analysis - use more reasonable intervals for plotting
# For 72 hours, use 5-minute intervals to reduce data size while maintaining resolution
time_points <- seq(
  from = global_time_min,
  to = global_time_max,
  by = "5 mins"  # 5-minute intervals for better performance
)

cat("Generated", length(time_points), "time points (5-minute intervals)\n")
cat("Time range:", as.character(min(time_points)), "to", as.character(max(time_points)), "\n")

# Create expanded grid for all combinations
group_vars <- df_filtered %>%
  distinct(plot_group, plot_title, support)

cat("Group combinations:", nrow(group_vars), "\n")

time_grid <- expand_grid(
  timestamp = time_points,
  group_vars
)

cat("Time grid dimensions:", nrow(time_grid), "rows\n")
cat("Expected data points per group:", length(time_points), "\n")

# =============================================================================
# CUMULATIVE READ COUNTING
# =============================================================================
cat("\n=== Calculating cumulative read counts ===\n")

# Round timestamps to match our time grid intervals (5-minute boundaries)
df_filtered_rounded <- df_filtered %>%
  mutate(
    timestamp_rounded = floor_date(timestamp, unit = "5 mins")
  )

# Aggregate reads by rounded time and group
df_agg <- df_filtered_rounded %>%
  group_by(plot_group, plot_title, support, timestamp_rounded) %>%
  summarise(reads_at_time = n(), .groups = "drop") %>%
  rename(timestamp = timestamp_rounded)

cat("Aggregated data points:", nrow(df_agg), "\n")

# Merge with time grid and fill missing values
df_complete <- time_grid %>%
  left_join(df_agg, by = c("plot_group", "plot_title", "support", "timestamp")) %>%
  mutate(reads_at_time = replace_na(reads_at_time, 0)) %>%
  arrange(plot_group, support, timestamp)

cat("Complete time grid after merge:", nrow(df_complete), "rows\n")

# Calculate cumulative sums - ensure they're calculated properly across the full time range
df_cumulative <- df_complete %>%
  group_by(plot_group, plot_title, support) %>%
  arrange(timestamp) %>%
  mutate(
    cumulative_reads = cumsum(reads_at_time),
    rel_time_hours = as.numeric(difftime(timestamp, global_time_min, units = "hours"))
  ) %>%
  ungroup()

cat("Cumulative data prepared with", nrow(df_cumulative), "data points\n")

# Debug: Show sample of cumulative data for verification
cat("Sample cumulative data (first group):\n")
sample_group <- df_cumulative %>%
  filter(plot_group == first(df_cumulative$plot_group)) %>%
  filter(support == "variant") %>%
  filter(reads_at_time > 0 | row_number() %% 100 == 1) %>%
  select(plot_group, support, rel_time_hours, reads_at_time, cumulative_reads) %>%
  head(10)
print(sample_group)

# =============================================================================
# PLOT GENERATION FUNCTIONS
# =============================================================================
cat("\n=== Defining plot generation functions ===\n")

create_variant_plot <- function(plot_data, plot_id, plot_title) {

  # Separate variant, reference, and other data
  variant_data <- plot_data %>% filter(support == "variant")
  reference_data <- plot_data %>% filter(support == "reference")
  other_data <- plot_data %>% filter(support == "other")

  # Calculate key metrics
  final_variant_count <- ifelse(nrow(variant_data) > 0, max(variant_data$cumulative_reads, na.rm = TRUE), 0)
  final_reference_count <- ifelse(nrow(reference_data) > 0, max(reference_data$cumulative_reads, na.rm = TRUE), 0)
  final_other_count <- ifelse(nrow(other_data) > 0, max(other_data$cumulative_reads, na.rm = TRUE), 0)

  # Handle -Inf case when no data exists
  final_variant_count <- ifelse(is.infinite(final_variant_count), 0, final_variant_count)
  final_reference_count <- ifelse(is.infinite(final_reference_count), 0, final_reference_count)
  final_other_count <- ifelse(is.infinite(final_other_count), 0, final_other_count)

  total_reads <- final_variant_count + final_reference_count + final_other_count

  if (total_reads == 0) {
    cat("Warning: No reads found for", plot_id, "\n")
    return(NULL)
  }

  vaf <- final_variant_count / total_reads

  # Debug: Show data range for this plot
  cat("  Plot data range:", min(plot_data$rel_time_hours), "to", max(plot_data$rel_time_hours), "hours\n")
  cat("  Variant data points:", nrow(variant_data), "Reference data points:", nrow(reference_data), "\n")

  # Create the plot with improved line rendering
  p <- ggplot(plot_data, aes(x = rel_time_hours, y = cumulative_reads, color = support)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(
      values = c("variant" = "#E31A1C", "reference" = "#1F78B4", "other" = "#33A02C"),
      labels = c("variant" = "Variant Support", "reference" = "Reference Support", "other" = "Other/Novel")
    ) +
    scale_x_continuous(
      name = "Time (hours)",
      breaks = seq(0, time_window_hours, by = 12),
      limits = c(0, time_window_hours),
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      name = "Cumulative Read Count",
      labels = comma_format(),
      expand = c(0.05, 0.05)
    ) +
    labs(
      title = plot_title,
      subtitle = glue("VAF: {percent(vaf, accuracy = 0.1)} | Total Reads: {comma(total_reads)} | Variant: {comma(final_variant_count)} | Reference: {comma(final_reference_count)} | Other: {comma(final_other_count)}"),
      color = "Support Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.major.y = element_line(colour = "grey90"),
      strip.text = element_text(face = "bold")
    )

  # Add detection thresholds
  thresholds <- c(1, 5, 10, 20, 50)
  threshold_colors <- viridis::viridis(length(thresholds), alpha = 0.3)

  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    if (final_variant_count >= thresh) {
      # Find first time point where threshold was reached
      first_detection <- variant_data %>%
        filter(cumulative_reads >= thresh) %>%
        slice_min(rel_time_hours, n = 1)

      if (nrow(first_detection) > 0) {
        p <- p + geom_vline(
          xintercept = first_detection$rel_time_hours,
          color = threshold_colors[i],
          linetype = "dashed",
          alpha = 0.7
        ) +
        annotate(
          "text",
          x = first_detection$rel_time_hours,
          y = max(plot_data$cumulative_reads) * 0.9,
          label = glue("{thresh} reads"),
          angle = 90,
          vjust = -0.5,
          size = 3,
          color = threshold_colors[i]
        )
      }
    }
  }

  return(p)
}

create_multiallelic_detail_plot <- function(plot_data, plot_id, plot_title) {

  # This function creates detailed plots for multi-allelic variants
  # showing individual alternate allele support

  # Get the original data with allele details and proper time grid alignment
  detailed_data <- df_filtered %>%
    filter(plot_group == plot_id) %>%
    filter(!is.na(alt_allele_num)) %>%
    mutate(timestamp_rounded = floor_date(timestamp, unit = "5 mins")) %>%
    group_by(alt_allele_num, timestamp_rounded) %>%
    summarise(reads_at_time = n(), .groups = "drop") %>%
    rename(timestamp = timestamp_rounded) %>%

    # Create full time grid for each allele
    complete(alt_allele_num, timestamp = time_points, fill = list(reads_at_time = 0)) %>%

    arrange(alt_allele_num, timestamp) %>%
    group_by(alt_allele_num) %>%
    mutate(
      cumulative_reads = cumsum(reads_at_time),
      rel_time_hours = as.numeric(difftime(timestamp, global_time_min, units = "hours"))
    ) %>%
    ungroup()

  if (nrow(detailed_data) == 0) {
    return(NULL)
  }

  # Create the detailed plot
  p <- ggplot(detailed_data, aes(x = rel_time_hours, y = cumulative_reads, color = alt_allele_num)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_viridis_d(name = "Alternate\nAllele") +
    scale_x_continuous(
      name = "Time (hours)",
      breaks = seq(0, time_window_hours, by = 12),
      limits = c(0, time_window_hours)
    ) +
    scale_y_continuous(
      name = "Cumulative Read Count",
      labels = comma_format()
    ) +
    labs(
      title = paste0(plot_title, " - Allele Detail"),
      subtitle = "Individual alternate allele support over time"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray50"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  return(p)
}

# =============================================================================
# PLOT GENERATION AND EXPORT
# =============================================================================
cat("\n=== Generating and saving plots ===\n")

# Split data by plot group
plot_data_list <- split(df_cumulative, df_cumulative$plot_group)

# Generate plots
successful_plots <- 0
failed_plots <- 0

for (plot_id in names(plot_data_list)) {

  plot_data <- plot_data_list[[plot_id]]
  plot_title <- unique(plot_data$plot_title)[1]

  cat("Creating plot for:", plot_id, "\n")

  tryCatch({

    # Main variant vs reference plot
    main_plot <- create_variant_plot(plot_data, plot_id, plot_title)

    if (!is.null(main_plot)) {

      # Clean filename
      safe_filename <- str_replace_all(plot_id, "[^A-Za-z0-9_.-]", "_")
      main_file <- file.path(output_dir, paste0(safe_filename, "_main.png"))

      ggsave(main_file, main_plot, width = 12, height = 8, dpi = 300)
      cat("  Saved main plot:", main_file, "\n")

      # Create detailed multi-allelic plot if applicable
      is_multiallelic_group <- any(df_filtered$plot_group == plot_id & df_filtered$is_multiallelic)

      if (is_multiallelic_group && separate_alleles) {
        detail_plot <- create_multiallelic_detail_plot(plot_data, plot_id, plot_title)

        if (!is.null(detail_plot)) {
          detail_file <- file.path(output_dir, paste0(safe_filename, "_allele_detail.png"))
          ggsave(detail_file, detail_plot, width = 12, height = 8, dpi = 300)
          cat("  Saved allele detail plot:", detail_file, "\n")
        }
      }

      successful_plots <- successful_plots + 1
    }

  }, error = function(e) {
    cat("  ERROR creating plot for", plot_id, ":", e$message, "\n")
    failed_plots <- failed_plots + 1
  })
}

# =============================================================================
# SUMMARY REPORT GENERATION
# =============================================================================
cat("\n=== Generating summary report ===\n")

# Create summary statistics
summary_stats <- df_filtered %>%
  group_by(plot_group, plot_title) %>%
  summarise(
    total_reads = n(),
    variant_reads = sum(support == "variant"),
    reference_reads = sum(support == "reference"),
    other_reads = sum(support_type == "other"),
    vaf = variant_reads / (variant_reads + reference_reads),
    is_multiallelic = any(is_multiallelic),
    unique_alleles = length(unique(alt_allele_num[!is.na(alt_allele_num)])),
    time_span_hours = as.numeric(difftime(max(timestamp), min(timestamp), units = "hours")),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads))

# Save summary report
summary_file <- file.path(output_dir, "plot_summary.tsv")
write_tsv(summary_stats, summary_file)

cat("Summary report saved:", summary_file, "\n")

# =============================================================================
# COMPLETION SUMMARY
# =============================================================================
cat("\n=== Plot Generation Complete ===\n")
cat("Results summary:\n")
cat("  - Successful plots:", successful_plots, "\n")
cat("  - Failed plots:", failed_plots, "\n")
cat("  - Output directory:", output_dir, "\n")
cat("  - Summary report:", summary_file, "\n")

if (successful_plots > 0) {
  cat("\nGenerated plot types:\n")
  cat("  - Main plots (variant vs reference):", successful_plots, "\n")
  if (separate_alleles) {
    allele_detail_files <- list.files(output_dir, pattern = "_allele_detail\\.png$")
    cat("  - Allele detail plots:", length(allele_detail_files), "\n")
  }
}

cat("\nConfiguration used:\n")
cat("  - Group multi-allelic variants:", group_multiallelic, "\n")
cat("  - Separate allele plots:", separate_alleles, "\n")
cat("  - Time window:", time_window_hours, "hours\n")
cat("  - Minimum reads per plot:", min_reads, "\n")

cat("\n=== Analysis Complete ===\n")
