#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(marquee)
  library(lubridate)
  library(patchwork)
  library(glue)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

# Define options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to input TSV file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: input file is required.", call. = FALSE)
}

varclock_results <- opt$input

# Load and prepare data
df <- readr::read_tsv(varclock_results, col_types = readr::cols(.default = "c")) %>%
  mutate(
    timestamp = ymd_hms(timestamp, tz = "UTC"),
    contains_variant = stringr::str_starts(allele_match, "VARIANT"),
    var_region = paste0(variant_description, " @ ", variant_pos),
    support = ifelse(contains_variant, "variant", "reference")
  )

# Aggregate + cumulative sum
agg_cum <- df %>%
  group_by(timestamp, var_region, region_name, support) %>%
  summarise(reads = n(), .groups = "drop") %>%
  arrange(var_region, region_name, support, timestamp) %>%
  group_by(var_region, region_name, support) %>%
  mutate(cumulative_reads = cumsum(reads)) %>%
  ungroup()

# Compute global time range
global_time_min <- min(agg_cum$timestamp, na.rm = TRUE)
global_time_max <- max(agg_cum$timestamp, na.rm = TRUE)

# Add a relative time column in hours
agg_cum <- agg_cum %>%
  mutate(rel_time_hours = as.numeric(difftime(timestamp, global_time_min, units = "hours")))

# Annotate with index + short ID
agg_cum <- agg_cum %>%
  group_by(region_name, var_region) %>%
  mutate(
    variant_index = stringr::str_pad(dense_rank(var_region), width = 4, pad = "0"),
    var_pos_clean = stringr::str_extract(var_region, "(?<=@\\s)\\d+"),
    var_pos_clean = ifelse(is.na(var_pos_clean), "posNA", var_pos_clean),
    plot_id = glue("{region_name}_{variant_index}_@{var_pos_clean}")
  ) %>%
  ungroup()

# Split by clean short plot ID
plot_data_list <- split(agg_cum, agg_cum$plot_id)

# Make plots
plot_list <- purrr::imap(plot_data_list, function(region_data, plot_id) {
  region_name <- unique(region_data$region_name)
  var_region  <- unique(region_data$var_region)

  var_region2 <- var_region |>
    gsub('<(.*)>', '\\1', x = _) |>
    gsub('\\s?@', '\n\n@', x = _) |>
    gsub('(.{175})', '\\1<br>', x =  _) |>
    gsub('@<br>', '@', x =  _) |>
    gsub('<br>([0-9]+)', '\\1', x = _)

  md_title <- glue("### {region_name} Reference/Variant Support in Reads\n\n<i>{var_region2}</i>")

  main <- ggplot(region_data, aes(x = rel_time_hours, y = cumulative_reads, color = support)) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c("variant" = "red", "reference" = "gray60"),
      labels = c("reference" = "Reference Reads", "variant" = "Variant Reads"),
      name = NULL
    ) +
    scale_x_continuous(
      name = "Sequencing Time (H)",
      breaks = seq(0, ceiling(as.numeric(difftime(global_time_max, global_time_min, units = "hours"))), by = 6),
      minor_breaks = seq(0, ceiling(as.numeric(difftime(global_time_max, global_time_min, units = "hours"))), by = 1)
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

  variant_thresholds <- c(1, 5, 10)
  threshold_colors <- c("1" = "red", "5" = "green", "10" = "blue")

  for (thresh in variant_thresholds) {
    time_point <- region_data %>%
      filter(support == "variant", cumulative_reads >= thresh) %>%
      arrange(timestamp) %>%
      slice(1) %>%
      mutate(rel_time_hours = as.numeric(difftime(timestamp, global_time_min, units = "hours"))) %>%
      pull(rel_time_hours)

    if (length(time_point) == 1) {
      main <- main +
        geom_vline(
          xintercept = as.numeric(time_point),
          linetype = "dotted",
          linewidth = 0.6,
          color = threshold_colors[as.character(thresh)]
        ) +
        annotate(
          "text",
          x = time_point,
          y = Inf,
          label = glue(
            "Time {floor(time_point)}:{stringr::str_pad(round((time_point %% 1) * 60), 2, pad = '0')} — ≥{thresh} variant read support"
          ),
          angle = 90,
          vjust = -0.5,
          hjust = 1,
          size = 6,
          color = threshold_colors[as.character(thresh)]
        )
    }
  }

  main
})

purrr::walk2(
  names(plot_list), plot_list,
  ~ggsave(glue("cumulative_reads_{.x}.png"), .y, width = 10, height = 6, dpi = 300)
)
