#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(dplyr)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Input TSV file path [required]", metavar = "FILE"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "variant_plot.png",
    help = "Output plot file path [default: %default]", metavar = "FILE"
  ),
  make_option(c("-t", "--title"),
    type = "character", default = "Variant-supporting Reads Over Time",
    help = "Plot title [default: %default]", metavar = "STRING"
  ),
  make_option(c("-w", "--width"),
    type = "integer", default = 10,
    help = "Plot width in inches [default: %default]", metavar = "INTEGER"
  ),
  make_option(c("-v", "--height"),
    type = "integer", default = 6,
    help = "Plot height in inches [default: %default]", metavar = "INTEGER"
  ),
  make_option(c("-b", "--binwidth"),
    type = "numeric", default = 60,
    help = "Histogram bin width in seconds [default: %default]", metavar = "NUMERIC"
  ),
  make_option("--dpi",
    type = "integer", default = 300,
    help = "Plot resolution in DPI [default: %default]", metavar = "INTEGER"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(
  usage = "Usage: %prog -i INPUT_FILE [OPTIONS]",
  description = "Create a histogram plot of variant-supporting reads over time from VarClock output",
  option_list = option_list,
  epilogue = "Example: Rscript simpleplot.R -i output.tsv -o my_plot.png -t 'My Variants' -w 12 -h 8"
)

opt <- parse_args(opt_parser)

# Check if input file is specified
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file (-i/--input) is required", call. = FALSE)
}

# Check if input file exists
if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call. = FALSE)
}

cat("Reading data from:", opt$input, "\n")

# Load data
df <- read_tsv(opt$input)

# Check if required columns exist
required_cols <- c("timestamp", "contains_variant", "variant_type")
missing_cols <- required_cols[!required_cols %in% colnames(df)]
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

# Convert timestamp to POSIXct (handle different timestamp formats)
df <- df %>%
  mutate(timestamp = case_when(
    timestamp == "NA" ~ as.POSIXct(NA),
    TRUE ~ as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  ))

# Filter for reads that contain the variant
df_variant <- df %>% filter(contains_variant == TRUE | contains_variant == "true")

# Check if we have any variant-containing reads
if (nrow(df_variant) == 0) {
  stop("No variant-supporting reads found in the data", call. = FALSE)
}

# Remove rows with NA timestamps for plotting
df_variant_clean <- df_variant %>% filter(!is.na(timestamp))

if (nrow(df_variant_clean) == 0) {
  warning("No valid timestamps found, cannot create time-based plot")
  stop("All timestamps are NA or invalid", call. = FALSE)
}

cat("Found", nrow(df_variant), "variant-supporting reads\n")
cat("Found", nrow(df_variant_clean), "reads with valid timestamps\n")

# Create plot
p <- ggplot(df_variant_clean, aes(x = timestamp, fill = variant_type)) +
  geom_histogram(binwidth = opt$binwidth, position = "stack") +
  labs(
    title = opt$title,
    x = "Time",
    y = "Read Count",
    fill = "Variant Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
cat("Saving plot to:", opt$output, "\n")
ggsave(opt$output, plot = p, width = opt$width, height = opt$height, dpi = opt$dpi)

cat("Plot saved successfully!\n")
