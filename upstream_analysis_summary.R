#!/usr/bin/env Rscript

# RNA-seq Pipeline Upstream Summary Report Generator
# Fixed version - no markdown strings, direct file writing

# ==================== ARGUMENT PARSING ====================
args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("\nUsage: Rscript generate_upstream_summary.R --metadata <file> --sample <file> [options]\n\n")
  cat("Required Arguments:\n")
  cat("  --metadata FILE    Metadata CSV file with 'sample' and 'group' columns\n")
  cat("  --sample FILE      Sample CSV file with sample-barcode mapping\n\n")
  cat("Optional Arguments:\n")
  cat("  --base-dir DIR     Base working directory (default: current directory)\n")
  cat("  --output-dir DIR   Output directory (default: output/upstream_analysis_summary)\n")
  cat("  --help             Show this help message\n\n")
  quit(status = 1)
}

parse_args <- function(args) {
  parsed <- list(
    metadata = NULL,
    sample = NULL,
    base_dir = ".",
    output_dir = "output/upstream_analysis_summary"
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--help") usage()
    
    if (arg %in% c("--metadata", "--sample", "--base-dir", "--output-dir")) {
      if (i == length(args) || grepl("^--", args[i+1])) {
        stop("Missing value for ", arg)
      }
      key <- gsub("-", "_", gsub("^--", "", arg))
      parsed[[key]] <- args[i + 1]
      i <- i + 2
    } else {
      stop("Unknown argument '", arg, "'")
    }
  }
  
  if (is.null(parsed$metadata) || is.null(parsed$sample)) {
    stop("Both --metadata and --sample arguments are required")
  }
  
  return(parsed)
}

params <- parse_args(args)
message("Parsed arguments: metadata=", params$metadata, ", sample=", params$sample)

# ==================== LIBRARY LOADING ====================
required_packages <- c("dplyr", "ggplot2", "knitr", "rmarkdown", "DT", "data.table", "stringr", "tidyr", "readr", "pagedown")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})

# ==================== PATH SETUP ====================
BASE_DIR <- normalizePath(params$base_dir)
METADATA_FILE <- file.path(BASE_DIR, params$metadata)
SAMPLE_FILE <- file.path(BASE_DIR, params$sample)
OUTPUT_DIR <- file.path(BASE_DIR, params$output_dir)
PIPELINE_OUTPUT <- file.path(BASE_DIR, "output")

# Validate
check_file <- function(path, name) {
  if (!file.exists(path)) stop(name, " file not found: ", path)
}
check_dir <- function(path, name) {
  if (!dir.exists(path)) stop(name, " directory not found: ", path)
}

check_file(METADATA_FILE, "Metadata")
check_file(SAMPLE_FILE, "Sample")
check_dir(PIPELINE_OUTPUT, "Pipeline output")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)

# ==================== HELPER FUNCTIONS ====================
count_reads_safe <- function(file_path) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NA)
  }
  
  cmd <- if (grepl("\\.gz$", file_path)) {
    paste("zcat", shQuote(file_path), "| wc -l")
  } else {
    paste("wc -l", shQuote(file_path))
  }
  
  tryCatch({
    lines <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)))
    return(ifelse(is.na(lines) || lines < 4, NA, lines / 4))
  }, error = function(e) {
    warning("Failed to count reads: ", e$message)
    return(NA)
  })
}

parse_star_log_safe <- function(log_file) {
  if (!file.exists(log_file)) {
    warning("STAR log not found: ", log_file)
    return(NULL)
  }
  
  tryCatch({
    content <- readLines(log_file)
    
    extract_value <- function(pattern, content, numeric = TRUE) {
      line <- grep(pattern, content, value = TRUE)
      if (length(line) == 0) return(NA)
      value <- str_extract(line, "\\d+\\.?\\d*")
      return(if (numeric) as.numeric(value) else value)
    }
    
    data.frame(
      input_reads = extract_value("Number of input reads", content),
      avg_input_length = extract_value("Average input read length", content),
      uniquely_mapped = extract_value("Uniquely mapped reads number", content),
      uniquely_mapped_pct = extract_value("Uniquely mapped reads %", content),
      multi_mapped = extract_value("Number of reads mapped to multiple loci", content),
      multi_mapped_pct = extract_value("% of reads mapped to multiple loci", content),
      avg_mapped_length = extract_value("Average mapped length", content),
      mismatch_rate = extract_value("Mismatch rate per base", content),
      deletion_rate = extract_value("Deletion rate per base", content),
      insertion_rate = extract_value("Insertion rate per base", content),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Failed to parse STAR log: ", e$message)
    return(NULL)
  })
}

parse_umi_log_safe <- function(log_file) {
  if (!file.exists(log_file)) {
    warning("UMI log not found: ", log_file)
    return(NULL)
  }
  
  tryCatch({
    content <- readLines(log_file)
    input_line <- grep("INFO Input Reads:", content, value = TRUE)
    output_line <- grep("INFO (Reads output|Output Reads):", content, value = TRUE)
    
    if (length(input_line) == 0 || length(output_line) == 0) {
      warning("Could not find read counts in UMI log: ", log_file)
      return(NULL)
    }
    
    input_reads <- as.numeric(str_extract(input_line, "\\d+"))
    output_reads <- as.numeric(str_extract(output_line, "\\d+"))
    
    if (is.na(input_reads) || is.na(output_reads) || input_reads == 0) {
      warning("Invalid read counts in UMI log: ", log_file)
      return(NULL)
    }
    
    data.frame(
      umi_input_reads = input_reads,
      umi_output_reads = output_reads,
      umi_retention_pct = (output_reads / input_reads) * 100,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Failed to parse UMI log: ", e$message)
    return(NULL)
  })
}

count_bam_reads <- function(bam_file) {
  if (!file.exists(bam_file)) {
    warning("BAM file not found: ", bam_file)
    return(NA)
  }
  
  tryCatch({
    cmd <- paste("samtools view -c", shQuote(bam_file))
    as.numeric(system(cmd, intern = TRUE))
  }, error = function(e) {
    warning("Failed to count BAM reads: ", e$message)
    return(NA)
  })
}

# ==================== EXTRACT METRICS ====================
message("=== Extracting pipeline metrics ===")

metadata <- tryCatch({
  fread(METADATA_FILE)
}, error = function(e) stop("Failed to load metadata: ", e$message))

if (!"sample" %in% names(metadata)) stop("Metadata must contain 'sample' column")
if (!"group" %in% names(metadata)) {
  message("Warning: Metadata missing 'group' column. Using 'undefined'.")
  metadata$group <- "undefined"
}

samples <- metadata$sample
summary_table <- data.frame(sample = samples, stringsAsFactors = FALSE)
summary_table <- left_join(summary_table, metadata, by = "sample")

# Parse metrics
demux_dir <- file.path(PIPELINE_OUTPUT, "demuxed")
if (dir.exists(demux_dir)) {
  for (s in samples) {
    file <- file.path(demux_dir, paste0(s, "_R1.fastq.gz"))
    summary_table[summary_table$sample == s, "demux_reads"] <- count_reads_safe(file)
  }
}

umi_dir <- file.path(PIPELINE_OUTPUT, "umi_extracted")
if (dir.exists(umi_dir)) {
  for (s in samples) {
    r1_file <- file.path(umi_dir, paste0(s, "_R1.fastq.gz"))
    log_file <- file.path(umi_dir, paste0(s, "_umi.log"))
    
    summary_table[summary_table$sample == s, "umi_reads"] <- count_reads_safe(r1_file)
    
    log_data <- parse_umi_log_safe(log_file)
    if (!is.null(log_data)) {
      idx <- summary_table$sample == s
      summary_table[idx, names(log_data)] <- log_data[1, ]
    }
  }
}

trim_dir <- file.path(PIPELINE_OUTPUT, "trimmed")
if (dir.exists(trim_dir)) {
  for (s in samples) {
    file <- file.path(trim_dir, paste0(s, "_R1.fastq"))
    summary_table[summary_table$sample == s, "trimmed_reads"] <- count_reads_safe(file)
  }
}

aligned_dir <- file.path(PIPELINE_OUTPUT, "aligned")
alignment_list <- list()
for (s in samples) {
  log_file <- file.path(aligned_dir, s, paste0(s, "_Log.final.out"))
  star_data <- parse_star_log_safe(log_file)
  if (!is.null(star_data)) {
    star_data$sample <- s
    alignment_list[[s]] <- star_data
  }
}
if (length(alignment_list) > 0) {
  alignment_df <- do.call(rbind, alignment_list)
  summary_table <- left_join(summary_table, alignment_df, by = "sample")
}

# UMI deduplication
for (s in samples) {
  bam_file <- file.path(aligned_dir, s, paste0(s, "_dedup.bam"))
  dedup_reads <- count_bam_reads(bam_file)
  summary_table[summary_table$sample == s, "dedup_reads"] <- dedup_reads
  
  aligned_reads <- summary_table[summary_table$sample == s, "uniquely_mapped"]
  if (!is.na(dedup_reads) && !is.na(aligned_reads) && aligned_reads > 0) {
  # dedup_reads should be <= aligned_reads, but if not, cap it
  dedup_reads_fixed <- min(dedup_reads, aligned_reads)
  dedup_rate <- max(0, min(100, ((aligned_reads - dedup_reads_fixed) / aligned_reads) * 100))
  summary_table[summary_table$sample == s, "dedup_rate_pct"] <- dedup_rate
  }
}

# Final gene counts
counts_file <- file.path(PIPELINE_OUTPUT, "analysis", "counts", "counts.matrix.txt")
if (file.exists(counts_file)) {
  tryCatch({
    counts_mat <- fread(counts_file, data.table = FALSE)
    rownames(counts_mat) <- counts_mat[, 1]
    counts_mat <- counts_mat[, -1]
    genes_detected <- apply(counts_mat, 2, function(x) sum(x > 0, na.rm = TRUE))
    
    for (s in samples) {
      if (s %in% names(genes_detected)) {
        summary_table[summary_table$sample == s, "genes_detected"] <- genes_detected[s]
      }
    }
  }, error = function(e) warning("Failed to parse counts: ", e$message))
}

# Calculate retention rates with proper NA handling
calculate_retention <- function(from, to, name) {
  pct <- (to / from) * 100
  pct[is.infinite(pct) | is.nan(pct) | is.na(pct)] <- 0  # Set invalid to 0%
  pct[pct > 100] <- 100  # Cap at 100%
  return(pct)
}

summary_table <- summary_table %>%
  mutate(
    demux_to_umi_pct = calculate_retention(demux_reads, umi_reads, "demux_to_umi"),
    umi_to_trimmed_pct = calculate_retention(umi_reads, trimmed_reads, "umi_to_trimmed"),
    trimmed_to_aligned_pct = calculate_retention(trimmed_reads, uniquely_mapped, "trimmed_to_aligned"),
    aligned_to_dedup_pct = calculate_retention(uniquely_mapped, dedup_reads, "aligned_to_dedup"),
    overall_retention_pct = calculate_retention(demux_reads, dedup_reads, "overall_retention")
  )

saveRDS(summary_table, file.path(OUTPUT_DIR, "upstream_summary_table.RDS"))
write.csv(summary_table, file.path(OUTPUT_DIR, "upstream_summary_table.csv"), row.names = FALSE)

message("✓ Extracted metrics for ", nrow(summary_table), " samples")

# ==================== CREATE VISUALIZATIONS ====================
message("=== Creating visualizations ===")

unique_groups <- unique(na.omit(summary_table$group))
group_colors <- setNames(
  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")[1:length(unique_groups)],
  unique_groups
)

# Read counts plot - linear scale with value labels
if (!all(is.na(summary_table$demux_reads))) {
  read_counts_long <- summary_table %>%
    select(sample, group, demux_reads, umi_reads, trimmed_reads, uniquely_mapped, dedup_reads) %>%
    pivot_longer(cols = -c(sample, group), names_to = "step", values_to = "reads") %>%
    filter(!is.na(reads)) %>%
    mutate(step = factor(step, levels = c("demux_reads", "umi_reads", "trimmed_reads", "uniquely_mapped", "dedup_reads"),
                         labels = c("Demultiplexed", "UMI Extracted", "Trimmed", "Aligned (Unique)", "Deduplicated")))
  
  # Calculate labels for top of bars (format as M for millions)
  read_counts_long <- read_counts_long %>%
    group_by(step) %>%
    mutate(label = ifelse(reads > 1e6, sprintf("%.1fM", reads / 1e6), 
                          ifelse(reads > 1e3, sprintf("%.0fK", reads / 1e3), 
                                 sprintf("%d", reads))))
  
  p_read_counts <- ggplot(read_counts_long, aes(x = sample, y = reads, fill = step)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    geom_text(aes(label = label, y = reads + max(reads, na.rm = TRUE) * 0.02), 
              position = position_dodge(width = 0.8), size = 3, angle = 90, hjust = 0) +
    scale_y_continuous(labels = scales::comma_format()) +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")) +
    scale_x_discrete(labels = ~ gsub("(_1|_2|_3)$", "", .x)) +
    labs(title = "Read Counts Through Pipeline Steps",
         x = "Sample", y = "Number of Reads", fill = "Processing Step") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  
  ggsave(file.path(OUTPUT_DIR, "figures", "read_counts.png"), p_read_counts, width = 12, height = 6, dpi = 300)
}

# Alignment rate plot
if (!all(is.na(summary_table$uniquely_mapped_pct))) {
  p_alignment <- ggplot(summary_table, aes(x = sample, y = uniquely_mapped_pct, fill = group)) +
    geom_bar(stat = "identity", show.legend = length(unique_groups) > 1) +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = ~ gsub("(_1|_2|_3)$", "", .x)) +
    labs(title = "STAR Alignment Rate by Sample",
         x = "Sample", y = "Uniquely Mapped Reads (%)",
         caption = "Percentage of trimmed reads uniquely mapped to genome") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "figures", "alignment_rate.png"), p_alignment, width = 10, height = 5, dpi = 300)
}

# Deduplication rate plot
if (!all(is.na(summary_table$dedup_rate_pct))) {
  p_dedup <- ggplot(summary_table, aes(x = sample, y = dedup_rate_pct, fill = group)) +
    geom_bar(stat = "identity", show.legend = length(unique_groups) > 1) +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = ~ gsub("(_1|_2|_3)$", "", .x)) +
    labs(title = "UMI Deduplication Rate",
         x = "Sample", y = "Duplicate Removal (%)",
         caption = "Percentage of aligned reads removed as PCR duplicates") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "figures", "dedup_rate.png"), p_dedup, width = 10, height = 5, dpi = 300)
}

# Retention flow plot - show actual median read counts for clarity
if (any(!is.na(summary_table[, c("demux_reads", "umi_reads", "trimmed_reads", "uniquely_mapped", "dedup_reads")]))) {
  # Calculate median reads for each step
  median_flow <- summary_table %>%
    summarise(
      Demultiplexed = median(demux_reads, na.rm = TRUE),
      `UMI Extracted` = median(umi_reads, na.rm = TRUE),
      Trimmed = median(trimmed_reads, na.rm = TRUE),
      `Aligned (Unique)` = median(uniquely_mapped, na.rm = TRUE),
      Deduplicated = median(dedup_reads, na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "Step", values_to = "Median_Reads")
  
  # Add filtering after pivot_longer:
  median_flow <- median_flow %>% 
    filter(!is.na(Median_Reads) & Median_Reads > 0) %>%
    mutate(Step = factor(Step, levels = c("Demultiplexed", "UMI Extracted", "Trimmed", "Aligned (Unique)", "Deduplicated")))
  
  # Remove NA values
  median_flow <- median_flow %>% filter(!is.na(Median_Reads) & Median_Reads > 0)
  
  p_retention <- ggplot(median_flow, aes(x = Step, y = Median_Reads, fill = Step)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1fM", Median_Reads / 1e6)), vjust = -0.5, size = 4) +
    scale_fill_brewer(palette = "Pastel1") +
    labs(title = "Median Read Counts Across Pipeline Steps",
         x = "Processing Step", y = "Median Reads (Millions)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "figures", "retention_flow.png"), p_retention, width = 10, height = 6, dpi = 300)
}

# ==================== GENERATE RMD REPORT ====================
message("=== Generating RMarkdown report ===")

# Open file connection for writing Rmd
rmd_file <- file.path(OUTPUT_DIR, "upstream_analysis_report.Rmd")
con <- file(rmd_file, "w", encoding = "UTF-8")

# Write YAML header
writeLines("---", con)
writeLines("title: \"RNA-seq Upstream Analysis Summary Report\"", con)
writeLines("author: \"3' RNA-seq Pipeline\"", con)
writeLines("date: \"`r Sys.Date()`\"", con)
writeLines("output:", con)
writeLines("  html_document:", con)
writeLines("    toc: true", con)
writeLines("    toc_depth: 3", con)
writeLines("    theme: cosmo", con)
writeLines("    highlight: tango", con)
writeLines("    code_folding: hide", con)
writeLines("    fig_width: 10", con)
writeLines("    fig_height: 6", con)
writeLines("---", con)
writeLines("", con)

# Write setup chunk
writeLines("```{r setup, include=FALSE}", con)
writeLines("knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.retina = 2)", con)
writeLines("library(DT)", con)
writeLines("library(dplyr)", con)
writeLines("library(knitr)", con)
writeLines("library(ggplot2)", con)
writeLines("", con)
writeLines(sprintf("OUTPUT_DIR <- \"%s\"", OUTPUT_DIR), con)
writeLines("summary_table <- readRDS(file.path(OUTPUT_DIR, \"upstream_summary_table.RDS\"))", con)
writeLines("", con)
writeLines(sprintf("# Define colors for %d groups", length(unique_groups)), con)
writeLines(sprintf("group_colors <- setNames(c(%s), c(%s))", 
                   toString(shQuote(group_colors, type = "cmd")),
                   toString(shQuote(names(group_colors), type = "cmd"))), con)
writeLines("```", con)
writeLines("", con)

# Write CSS
writeLines("<style>", con)
writeLines("body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }", con)
writeLines("h1, h2, h3 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 5px; }", con)
writeLines(".table th { background-color: #34495e; color: white; }", con)
writeLines(".alert { padding: 10px; margin: 10px 0; border-radius: 5px; border-left: 5px solid; }", con)
writeLines(".alert-success { background-color: #d4edda; border-color: #28a745; }", con)
writeLines(".alert-warning { background-color: #fff3cd; border-color: #ffc107; }", con)
writeLines("</style>", con)
writeLines("", con)

# Create content sections as separate lines
content_sections <- list(
  "# Pipeline Overview",
  "",
  "This report summarizes the **upstream analysis** performance for the RNA-seq dataset.",
  "",
  "## Experimental Design",
  "",
  "```{r metadata-tab, echo=FALSE}",
  "datatable(summary_table[, c(\"sample\", \"group\")],",
  "          options = list(dom = \"t\", ordering = FALSE, pageLength = 10),",
  "          rownames = FALSE,",
  "          caption = \"Sample Metadata\")",
  "```",
  "",
  "```{r qc-checks, results=\"asis\"}",
  "qc_alerts <- c()",
  "if (mean(summary_table$uniquely_mapped_pct, na.rm = TRUE) < 70) {",
  "  qc_alerts <- c(qc_alerts, \"WARNING: Alignment Rate Low - Average uniquely mapped reads < 70%\")",
  "}",
  "if (mean(summary_table$dedup_rate_pct, na.rm = TRUE) > 50) {",
  "  qc_alerts <- c(qc_alerts, \"WARNING: High Deduplication - >50% PCR duplication detected\")",
  "}",
  "if (mean(summary_table$genes_detected, na.rm = TRUE) < 5000) {",
  "  qc_alerts <- c(qc_alerts, \"WARNING: Low Gene Detection - < 5,000 genes detected on average\")",
  "}",
  "",
  "if (length(qc_alerts) > 0) {",
  "  cat(\"<div class=\\\"alert alert-warning\\\"><h4>Quality Alerts</h4>\")",
  "  cat(paste(qc_alerts, collapse = \"<br>\"))",
  "  cat(\"</div>\")",
  "} else {",
  "  cat(\"<div class=\\\"alert alert-success\\\"><h4>All QC Metrics Passed</h4>",
  "      All samples meet quality thresholds for downstream analysis.</div>\")",
  "}",
  "```",
  "",
  "## Key Metrics Summary",
  "",
  "```{r key-metrics, echo=FALSE}",
  "library(dplyr)",
  "total_reads <- sum(summary_table$demux_reads, na.rm = TRUE)",
  "avg_alignment <- mean(summary_table$uniquely_mapped_pct, na.rm = TRUE)",
  "avg_dedup <- mean(summary_table$dedup_rate_pct, na.rm = TRUE)",
  "avg_genes <- mean(summary_table$genes_detected, na.rm = TRUE)",
  "",
  "metrics_df <- data.frame(",
  "  Metric = c(\"Total Reads Processed\", \"Avg. Alignment Rate\", \"Avg. Duplicate Removal\", \"Avg. Genes Detected\"),",
  "  Value = c(",
  "    sprintf(\"%.2f M\", total_reads / 1e6),",
  "    sprintf(\"%.1f%%\", avg_alignment),",
  "    sprintf(\"%.1f%%\", avg_dedup),",
  "    sprintf(\"%.0f\", avg_genes)",
  "  )",
  ")",
  "",
  "knitr::kable(metrics_df, align = \"lrr\", caption = \"Pipeline Performance Summary\")",
  "```",
  "",
  "## Read Processing",
  "",
  "### Read Counts Through Pipeline Steps",
  "",
  "```{r fig1, fig.cap=\"Figure 1: Read counts at each processing step. Values are not log10 scaled for clarity.\", echo=FALSE}",
  "include_graphics(file.path(OUTPUT_DIR, \"figures\", \"read_counts.png\"))",
  "```",
  "",
  "### Retention Rates Between Steps",
  "",
  "```{r fig-retention, fig.cap=\"Figure 2: Average percentage of reads retained between pipeline steps.\", echo=FALSE}",
  "include_graphics(file.path(OUTPUT_DIR, \"figures\", \"retention_flow.png\"))",
  "```",
  "",
  "## Alignment Quality",
  "",
  "```{r fig2, fig.cap=\"Figure 3: Percentage of trimmed reads that mapped uniquely to the genome.\", echo=FALSE}",
  "if (file.exists(file.path(OUTPUT_DIR, \"figures\", \"alignment_rate.png\"))) {",
  "  include_graphics(file.path(OUTPUT_DIR, \"figures\", \"alignment_rate.png\"))",
  "}",
  "```",
  "",
  "### STAR Alignment Details",
  "",
  "```{r star-table, echo=FALSE}",
  "required_cols <- c(\"sample\", \"uniquely_mapped_pct\", \"avg_mapped_length\", \"mismatch_rate\")",
  "available_cols <- intersect(required_cols, names(summary_table))",
  "",
  "if (length(available_cols) > 1) {",
  "  knitr::kable(",
  "    summary_table[, available_cols] %>% arrange(desc(uniquely_mapped_pct)),",
  "    digits = 2,",
  "    caption = \"STAR Alignment Quality Metrics\"",
  "  )",
  "} else {",
  "  cat(\"No alignment metrics available.\")",
  "}",
  "```",
  "",
  "## UMI Deduplication",
  "",
  "```{r fig3, fig.cap=\"Figure 4: Percentage of duplicate reads removed by UMI deduplication.\", echo=FALSE}",
  "if (file.exists(file.path(OUTPUT_DIR, \"figures\", \"dedup_rate.png\"))) {",
  "  include_graphics(file.path(OUTPUT_DIR, \"figures\", \"dedup_rate.png\"))",
  "}",
  "```",
  "",
  "### Deduplication Summary",
  "",
  "```{r dedup-table, echo=FALSE}",
  "if (\"dedup_rate_pct\" %in% names(summary_table)) {",
  "  dedup_summary <- summary_table %>%",
  "    select(sample, group, uniquely_mapped, dedup_reads, dedup_rate_pct) %>%",
  "    arrange(desc(dedup_rate_pct))",
  "",
  "  knitr::kable(dedup_summary, caption = \"UMI Deduplication Summary\", digits = 0)",
  "} else {",
  "  cat(\"No deduplication metrics available.\")",
  "}",
  "```",
  "",
  "## Quantification",
  "",
"### Genes Detected Per Sample",
"",
"```{r barplot-genes, echo=FALSE}",
"# Filter out NA values before plotting",
"plot_data <- summary_table %>% filter(!is.na(genes_detected))",
"if (nrow(plot_data) > 0) {",
"  p_genes <- ggplot(plot_data, aes(x = sample, y = genes_detected, fill = group)) +",
"    geom_bar(stat = \"identity\", show.legend = length(unique(plot_data$group)) > 1) +",
"    scale_fill_manual(values = group_colors) +",
"    scale_x_discrete(labels = function(x) gsub(\"(_1|_2|_3)$\", \"\", x)) +",
"    labs(title = \"Genes Detected Per Sample\",",
"         x = \"Sample\", y = \"Number of Genes (count > 0)\",",
"         caption = \"Total genes with non-zero counts after quantification\") +",
"    theme_minimal(base_size = 12) +",
"    theme(axis.text.x = element_text(angle = 45, hjust = 1))",
"",
"  print(p_genes)",
"} else {",
"  cat(\"No gene detection data available for plotting.\")",
"}",
"```",
  "",
  "### Gene Detection Statistics",
  "",
  "```{r gene-stats, echo=FALSE}",
  "if (!all(is.na(summary_table$genes_detected))) {",
  "  gene_summary <- summary_table %>%",
  "    filter(!is.na(genes_detected)) %>%",
  "    group_by(group) %>%",
  "    summarise(",
  "      Count = n(),",
  "      Avg_Genes = mean(genes_detected),",
  "      Median_Genes = median(genes_detected),",
  "      Min_Genes = min(genes_detected),",
  "      Max_Genes = max(genes_detected)",
  "    )",
  "",
  "  knitr::kable(gene_summary, caption = \"Gene Detection by Group\", digits = 0)",
  "}",
  "```",
  "",
  "## Raw Data",
  "",
  "### Complete Metrics Table",
  "",
  "```{r complete-table, echo=FALSE}",
  "library(DT)",
  "datatable(",
  "  summary_table,",
  "  options = list(",
  "    pageLength = 10,",
  "    dom = \"Bfrtip\",",
  "    buttons = c(\"copy\", \"csv\", \"excel\"),",
  "    scrollX = TRUE",
  "  ),",
  "  rownames = FALSE,",
  "  caption = \"Complete Upstream Analysis Metrics\"",
  ")",
  "```",
  "",
  "### Data Files",
  "",
  "Full summary table is available as:",
  "- **CSV**: `upstream_summary_table.csv` (for Excel/R)",
  "- **RDS**: `upstream_summary_table.RDS` (for R only)",
  "",
  "---",
  "",
  sprintf("**Report Generated:** `r Sys.Date()`  "),
  sprintf("**Base Directory:** '%s'  ", BASE_DIR),
  sprintf("**Metadata:** '%s'  ", basename(METADATA_FILE)),
  sprintf("**Sample File:** '%s'", basename(SAMPLE_FILE))
)

# Write content sections
for (line in content_sections) {
  writeLines(line, con)
}

close(con)

message("✓ RMarkdown file written to: ", rmd_file)

# ==================== RENDER REPORTS ====================
message("=== Rendering HTML report ===")
html_file <- file.path(OUTPUT_DIR, "upstream_analysis_report.html")

tryCatch({
  render(rmd_file, output_file = html_file, quiet = TRUE)
  message("✅ HTML report: ", html_file)
}, error = function(e) {
  stop("Failed to render HTML report: ", e$message)
})

message("=== Converting to PDF ===")
pdf_file <- file.path(OUTPUT_DIR, "upstream_analysis_report.pdf")

if (suppressWarnings(require("pagedown", quietly = TRUE))) {
  tryCatch({
    chrome_print(html_file, output = pdf_file, timeout = 300)
    message("✅ PDF report: ", pdf_file)
  }, error = function(e) {
    warning("Pagedown PDF conversion failed: ", e$message)
  })
} else {
  message("Note: Install 'pagedown' for better PDF rendering: install.packages('pagedown')")
}

message("=== Complete! Reports saved to: ", OUTPUT_DIR)