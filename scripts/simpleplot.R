#!/usr/bin/env Rscript

# =============================================================================
# VARCLOCK Simple Plot Generator
# =============================================================================
# This script generates cumulative read count plots for variant analysis,
# showing reference vs variant support over time with detection thresholds.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)    # Command line argument parsing
  library(marquee)     # Enhanced text rendering in plots
  library(lubridate)   # Date/time manipulation
  library(patchwork)   # Plot composition (loaded but not used in current version)
  library(glue)        # String interpolation
  library(magrittr)    # Pipe operators
  library(dplyr)       # Data manipulation
  library(ggplot2)     # Plotting
  library(readr)       # Reading data files
  library(stringr)     # String manipulation
  library(purrr)       # Functional programming tools
  library(tidyr)       # Data tidying
})

# =============================================================================
# COMMAND LINE ARGUMENT SETUP
# =============================================================================
cat("Setting up command line arguments...\n")

# Define options
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to input TSV file", metavar = "FILE"),
  make_option(c("-t", "--timing"), type = "character", default = NULL,
              help = "Optional file with Earliest/Latest timestamps", metavar = "TIMING")
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

# Extract base filename without extension
input_basename <- basename(varclock_results)
input_basename_noext <- sub("\\.[^.]*$", "", input_basename)

cat("Input file:", varclock_results, "\n")
cat("Timing file:", ifelse(is.null(timing_file), "None", timing_file), "\n")

# =============================================================================
# DATA LOADING AND INITIAL PROCESSING
# =============================================================================
cat("\n=== Loading and preparing data ===\n")

# Load the main varclock results file
cat("Reading TSV file:", varclock_results, "\n")
df <- readr::read_tsv(varclock_results, col_types = readr::cols(.default = "c"))
cat("Initial data dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Column names:", paste(colnames(df), collapse = ", "), "\n")

# Transform and clean the data - PRESERVE MINUTE PRECISION
df <- df %>%
  mutate(
    # Parse timestamp and round to nearest minute (preserve precision)
    timestamp_raw = ymd_hms(timestamp, tz = "UTC"),
    timestamp = round_date(timestamp_raw, unit = "secs"),  # Round to minute boundary
    # Determine if this read contains a variant
    contains_variant = stringr::str_starts(allele_match, "VARIANT"),
    # Create a readable variant region identifier
    var_region = paste0(variant_description, " @ ", variant_pos),
    # Classify support type (variant or reference)
    support = ifelse(contains_variant, "variant", "reference")
  )

cat("After initial processing:\n")
cat("  - Raw timestamp range:", as.character(min(df$timestamp_raw, na.rm = TRUE)),
    "to", as.character(max(df$timestamp_raw, na.rm = TRUE)), "\n")
cat("  - Floored timestamp range:", as.character(min(df$timestamp, na.rm = TRUE)),
    "to", as.character(max(df$timestamp, na.rm = TRUE)), "\n")
cat("  - Unique variant regions:", length(unique(df$var_region)), "\n")
cat("  - Support distribution:\n")
print(table(df$support))

# Show rounding effects
cat("  - Rounding analysis (first 10 rows):\n")
rounding_check <- df %>%
  select(timestamp_raw, timestamp) %>%
  mutate(
    time_diff_minutes = as.numeric(difftime(timestamp_raw, timestamp, units = "secs"))
  ) %>%
  head(10)
print(rounding_check)

# =============================================================================
# GROUPING AND UUID ASSIGNMENT
# =============================================================================
cat("\n=== Creating group identifiers ===\n")

# Create unique identifiers for each variant region + region name combination
df <- df %>%
  group_by(var_region, region_name) %>%
  mutate(
    # Create a unique identifier for each group
    uuid = paste0(cur_group_id())
  ) %>%
  ungroup()

cat("Created", length(unique(df$uuid)), "unique variant-region combinations\n")
cat("Group summary:\n")
group_summary <- df %>%
  group_by(uuid, var_region, region_name) %>%
  summarise(read_count = n(), .groups = "drop")
print(head(group_summary))

# =============================================================================
# TIME RANGE DETERMINATION
# =============================================================================
cat("\n=== Determining time range ===\n")

# Determine global time minimum from timing file or data
if (!is.null(timing_file) && file.exists(timing_file)) {
  cat("Reading timing file:", timing_file, "\n")
  timing_lines <- readLines(timing_file)
  cat("Timing file contents:\n")
  cat(paste(timing_lines, collapse = "\n"), "\n")

  earliest_line <- timing_lines %>%
    str_subset("^Earliest\\s*:") %>%
    str_remove("Earliest\\s*:\\s*")

  if (length(earliest_line) > 0) {
    global_time_min_raw <- ymd_hms(earliest_line, tz = "UTC")
    # Floor the timing file time to hour boundary for consistency
    global_time_min <- floor_date(global_time_min_raw, unit = "secs")
    cat("Using timing file for start time (raw):", as.character(global_time_min_raw), "\n")
    cat("Using timing file for start time (floored):", as.character(global_time_min), "\n")
  } else {
    cat("Warning: Could not find 'Earliest:' line in timing file, using data minimum\n")
    global_time_min <- floor_date(min(df$timestamp, na.rm = TRUE), unit = "secs")
  }
} else {
  cat("No timing file provided or file doesn't exist, using data minimum\n")
  global_time_min <- floor_date(min(df$timestamp, na.rm = TRUE), unit = "secs")
}

# Set global time maximum to 72 hours after start, floored to hour
global_time_max <- global_time_min + hours(72)

cat("Time analysis window:\n")
cat("  - Start:", as.character(global_time_min), "\n")
cat("  - End:", as.character(global_time_max), "\n")
cat("  - Duration:", as.numeric(difftime(global_time_max, global_time_min, units = "secs")), "seconds\n")

# =============================================================================
# TIME GRID GENERATION WITH PROPER BOUNDARIES
# =============================================================================
cat("\n=== Generating time grid with proper boundaries ===\n")

# Generate full hourly sequence - ensure we capture full range
full_seconds <- seq(
  from = global_time_min,
  to = global_time_max,
  by = "1 sec"
)

cat("Generated", length(full_seconds), "hourly time points\n")
cat("First time point:", as.character(full_seconds[1]), "\n")
cat("Last time point:", as.character(full_seconds[length(full_seconds)]), "\n")

# Get unique combinations of grouping variables
group_vars <- df %>%
  distinct(uuid, var_region, region_name, support)

cat("Unique group combinations:", nrow(group_vars), "\n")
print(head(group_vars))

# Expand grid: all time points x group combinations
cat("Creating complete time x group grid...\n")
time_grid <- tidyr::expand_grid(
  timestamp = full_seconds,
  group_vars
)
cat("\n=======================================\n")
cat("\n=========  Testing Time Grid  =========\n")
cat("\n=======================================\n")
print(time_grid)

cat("Time grid dimensions:", nrow(time_grid), "rows\n")
cat("This represents", length(full_seconds), "time points ×", nrow(group_vars), "group combinations\n")

# =============================================================================
# HOURLY AGGREGATION (NO ADDITIONAL ROUNDING)
# =============================================================================
cat("\n=== Aggregating data to hourly bins ===\n")

# Aggregate minute-level data to hourly bins
df_hourly <- df %>%
  mutate(second_floor = floor_date(timestamp, unit = "secs")) %>%
  group_by(second_floor, var_region, region_name, support) %>%
  summarise(reads = n(), .groups = "drop") %>%
  rename(timestamp = second_floor)

cat("Hourly aggregation complete:\n")
cat("  - Aggregated data dimensions:", nrow(df_hourly), "rows\n")
cat("  - Hours with data:", length(unique(df_hourly$timestamp)), "\n")
cat("  - Total reads in aggregated data:", sum(df_hourly$reads), "\n")

# Show sample of hourly data
cat("Sample of hourly data:\n")
print(head(df_hourly))

# Check for any timestamps outside our expected range
out_of_range <- df_hourly %>%
  filter(timestamp < global_time_min | timestamp > global_time_max)
if (nrow(out_of_range) > 0) {
  cat("WARNING: Found timestamps outside expected range:\n")
  print(out_of_range)
}

# =============================================================================
# CUMULATIVE CALCULATION WITH PROPER BOUNDARIES
# =============================================================================
cat("\n=== Computing cumulative sums with proper boundaries ===\n")

# Join with grid and compute cumulative sums
agg_cum <- time_grid %>%
  left_join(df_hourly, by = c("timestamp", "var_region", "region_name", "support")) %>%
  mutate(reads = replace_na(reads, 0)) %>%
  arrange(var_region, region_name, support, timestamp) %>%
  group_by(var_region, region_name, support) %>%
  mutate(
    cumulative_reads = cumsum(reads),
    rel_time_seconds = as.numeric(difftime(timestamp, global_time_min, units = "secs"))
  ) %>%
  ungroup()

cat("Cumulative calculation complete:\n")
cat("  - Final data dimensions:", nrow(agg_cum), "rows\n")
cat("  - Zero-read entries:", sum(agg_cum$reads == 0), "\n")
cat("  - Non-zero entries:", sum(agg_cum$reads > 0), "\n")

# Verify boundary conditions
cat("Boundary condition verification:\n")
boundary_check <- agg_cum %>%
  group_by(var_region, region_name, support) %>%
  summarise(
    min_time = min(rel_time_seconds),
    max_time = max(rel_time_seconds),
    min_cumulative = min(cumulative_reads),
    max_cumulative = max(cumulative_reads),
    .groups = "drop"
  )
cat("All series time ranges:\n")
print(boundary_check)

# Check that all series start at time 0
if (any(boundary_check$min_time != 0)) {
  cat("WARNING: Some series don't start at time 0!\n")
}

# Check that all series end at the same time
expected_max_time <- as.numeric(difftime(global_time_max, global_time_min, units = "secs"))
if (any(boundary_check$max_time != expected_max_time)) {
  cat("WARNING: Some series don't end at expected time", expected_max_time, "!\n")
}

# =============================================================================
# PLOT IDENTIFIER GENERATION
# =============================================================================
cat("\n=== Generating plot identifiers ===\n")

# Annotate with index + short ID for plot organization
agg_cum <- agg_cum %>%
  group_by(uuid, region_name, var_region) %>%
  mutate(
    # Create padded variant index for sorting
    variant_index = stringr::str_pad(dense_rank(var_region), width = 4, pad = "0"),
    # Extract position from variant region string
    var_pos_clean = stringr::str_extract(var_region, "(?<=@\\s)\\d+"),
    var_pos_clean = ifelse(is.na(var_pos_clean), "posNA", var_pos_clean),
    # Create clean plot identifier
    plot_id = glue("{uuid}_{region_name}_{variant_index}_@{var_pos_clean}")
  ) %>%
  ungroup()

cat("Plot identifiers created:\n")
unique_plots <- unique(agg_cum$plot_id)
cat("  - Number of unique plots:", length(unique_plots), "\n")
cat("  - Plot IDs:\n")
for (i in seq_along(unique_plots)) {
  cat("    ", i, ":", unique_plots[i], "\n")
}

print(agg_cum)

# =============================================================================
# DATA SPLITTING FOR PLOTTING
# =============================================================================
cat("\n=== Splitting data for individual plots ===\n")

# Split by clean short plot ID
plot_data_list <- split(agg_cum, agg_cum$plot_id)

cat("Data split into", length(plot_data_list), "plot datasets\n")

# Show summary of each plot dataset
for (plot_id in names(plot_data_list)) {
  data <- plot_data_list[[plot_id]]

  # Check boundary conditions for each plot
  variant_data <- data %>% filter(support == "variant")
  reference_data <- data %>% filter(support == "reference")

  cat("Plot:", plot_id, "\n")
  cat("  - Data points per series:", nrow(variant_data), "\n")
  cat("  - Variant series: starts at hour", min(variant_data$rel_time_seconds),
      "with", min(variant_data$cumulative_reads), "reads\n")
  cat("  - Variant series: ends at hour", max(variant_data$rel_time_seconds),
      "with", max(variant_data$cumulative_reads), "reads\n")
  cat("  - Reference series: starts at hour", min(reference_data$rel_time_seconds),
      "with", min(reference_data$cumulative_reads), "reads\n")
  cat("  - Reference series: ends at hour", max(reference_data$rel_time_seconds),
      "with", max(reference_data$cumulative_reads), "reads\n")
}


# =============================================================================
# PLOT GENERATION
# =============================================================================
cat("\n=== Generating plots ===\n")

# Make plots
plot_list <- purrr::imap(plot_data_list, function(region_data, plot_id) {
  cat("Creating plot for:", plot_id, "\n")

  print(region_data)


  # Extract metadata for this plot
  region_name <- unique(region_data$region_name)
  var_region  <- unique(region_data$var_region)

  cat("  - Region name:", region_name, "\n")
  cat("  - Variant region:", var_region, "\n")

  # Verify this plot data has proper boundaries
  time_range <- range(region_data$rel_time_seconds)
  cat("  - Time range: [", time_range[1], ",", time_range[2], "] Seconds\n")

  # Check start and end points
  start_points <- region_data %>%
    filter(rel_time_seconds == min(rel_time_seconds)) %>%
    select(support, cumulative_reads)
  cat("  - Start points:\n")
  print(start_points)

  end_points <- region_data %>%
    filter(rel_time_seconds == max(rel_time_seconds)) %>%
    select(support, cumulative_reads)
  cat("  - End points:\n")
  print(end_points)

  # Format variant region for display (clean up formatting)
  var_region2 <- var_region |>
    gsub('<(.*)>', '\\1', x = _) |>           # Remove < > brackets
    gsub('\\s?@', '\n\n@', x = _) |>          # Add line breaks before @
    gsub('(.{175})', '\\1<br>', x =  _) |>    # Add line breaks every 175 chars
    gsub('@<br>', '@', x =  _) |>             # Fix @ after line breaks
    gsub('<br>([0-9]+)', '\\1', x = _)        # Fix numbers after line breaks

  # Calculate VAF statistics
  variant_max <- max(region_data$cumulative_reads[region_data$support == "variant"], na.rm = TRUE)
  reference_max <- max(region_data$cumulative_reads[region_data$support == "reference"], na.rm = TRUE)
  total_reads <- variant_max + reference_max
  vaf <- round(100 * variant_max / total_reads, 2)

  # Create formatted title with filename and VAF details
  md_title <- glue("### {input_basename_noext} - {region_name} Reference/Variant Support in Reads\n\n<i>{var_region2}</i>\n\nVAF: {vaf}% ({variant_max} variant / {reference_max} reference / {total_reads} total reads)")

  # Create the main plot
  main <- ggplot(
    region_data,
    aes(x = rel_time_seconds/3600, y = cumulative_reads, color = support)) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c("variant" = "orange", "reference" = "gray60"),
      labels = c("reference" = "Reference Reads", "variant" = "Variant Reads"),
      name = NULL
    ) +
   	scale_x_continuous(
       name = "Sequencing Time (H)",
    		breaks = seq(0, 72, by = 12),
    		minor_breaks = seq(0, 72, by = 6)
   	) +
    labs(y = "Cumulative Read Count", color = "Support Type") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.margin = margin(t = -10, b = 10),
      legend.box.margin = margin(b = -10),
      legend.text = element_text(size = 15)
    ) +
    ggtitle(md_title) +
    theme(plot.title = element_marquee(
      size = 10,
      width = 1,
      lineheight = 1.2,
      margin = margin(b = 10)
    ))

  # =============================================================================
  # THRESHOLD ANNOTATIONS
  # =============================================================================
  # Add vertical lines for variant detection thresholds
  variant_thresholds <- c(1, 5, 10)
  threshold_colors <- c("1" = "red", "5" = "green", "10" = "blue")

  cat("  - Adding threshold annotations...\n")

  # Collect all threshold times first to handle x-axis overlaps
  # Use original minute-level data to get precise timestamps
  threshold_times <- list()

  # Get the variant/region info for this plot
  plot_var_region <- unique(region_data$var_region)[1]  # Take first in case of multiple
  plot_region_name <- unique(region_data$region_name)[1]  # Take first in case of multiple

  for (thresh in variant_thresholds) {
    # Find precise minute-level timestamp when threshold was reached
    # First get cumulative counts at minute level
    minute_cumulative <- df %>%
      filter(var_region == plot_var_region, region_name == plot_region_name, support == "variant") %>%
      arrange(timestamp) %>%
      mutate(
        cumulative_reads = row_number(),
        rel_time_seconds = as.numeric(difftime(timestamp, global_time_min, units = "secs"))
      ) %>%
      filter(cumulative_reads >= thresh) %>%
      slice(1)

    if (nrow(minute_cumulative) == 1 && !is.na(minute_cumulative$rel_time_seconds)) {
      time_point <- minute_cumulative$rel_time_seconds
      threshold_times[[as.character(thresh)]] <- time_point
      cat("    - Threshold", thresh, "reached at hour", round(time_point/3600, 3), "\n")
    } else {
      cat("    - Threshold", thresh, "never reached\n")
    }
  }

  # Add annotations with x-nudging to avoid overlap
  if (length(threshold_times) > 0) {
    # Create data frame for easier manipulation
    thresh_df <- data.frame(
      threshold = as.numeric(names(threshold_times)),
      time_point = unlist(threshold_times),
      stringsAsFactors = FALSE
    )

    # Safety check - ensure we have valid data
    if (nrow(thresh_df) > 0 && !any(is.na(thresh_df$time_point))) {
      thresh_df <- thresh_df %>% arrange(time_point, threshold)
    } else {
      # Create empty data frame with proper structure if no valid thresholds
      thresh_df <- data.frame(threshold=numeric(0), time_point=numeric(0), x_nudge=numeric(0))
    }

    # Calculate x-nudges to avoid overlap when times are too close
    thresh_df$x_nudge <- 0
    if (nrow(thresh_df) > 1) {
      for (i in 2:nrow(thresh_df)) {
        time_diff <- thresh_df$time_point[i] - thresh_df$time_point[i-1]
        if (!is.na(time_diff) && time_diff < 3600) {  # If within 1 hour and not NA, nudge
          thresh_df$x_nudge[i] <- 1.5  # Nudge right
        }
      }
    }

    # Add vertical lines and annotations
    for (i in 1:nrow(thresh_df)) {
      thresh <- thresh_df$threshold[i]
      time_point <- thresh_df$time_point[i]
      x_nudge <- thresh_df$x_nudge[i]

      main <- main +
        geom_vline(
          xintercept = as.numeric(time_point/3600),
          linetype = "dotted",
          linewidth = 0.6,
          color = threshold_colors[as.character(thresh)]
        ) +
        annotate(
          "text",
          x = time_point/3600 + x_nudge,
          y = Inf,
          label =
          glue(
            "Time {sprintf('%02d:%02d',
                            floor(time_point / 3600),
                            ceiling((time_point %% 3600) / 60)
             )} — ≥{thresh} variant read support"
          ),
          angle = 90,
          vjust = -0.5,
          hjust = 1,
          size = 6,
          color = threshold_colors[as.character(thresh)]
        )
    }
  }

  cat("  - Plot creation complete for:", plot_id, "\n")
  main
})

# =============================================================================
# PLOT SAVING
# =============================================================================
cat("\n=== Saving plots ===\n")

# Save plots
purrr::walk2(
  names(plot_list), plot_list,
  function(plot_name, plot_obj) {
    filename <- glue("cumulative_reads_{plot_name}.png")
    cat("Saving plot:", filename, "\n")
    ggsave(filename, plot_obj, width = 15, height = 9, dpi = 600)
    cat("  - Saved successfully:", filename, "\n")
  }
)

cat("\n=== SCRIPT COMPLETE ===\n")
cat("Generated", length(plot_list), "plots successfully\n")
cat("All plots saved to current directory with pattern: cumulative_reads_*.png\n")
