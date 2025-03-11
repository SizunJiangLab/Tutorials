#!/usr/bin/env Rscript
#
# Rscript render_CosMxQC_report.R "MyQC.html"
# Rscript render_CosMxQC_report.R "output/Report_2024.html" "CosMxQC_report.Rmd"

args <- commandArgs(trailingOnly = TRUE)

# parse parameters
if(length(args) < 1) {
  stop("Usage: Rscript render_report.R <output_filename> [input_file.Rmd]")
}

output_file <- args[1]
input_file <- ifelse(length(args) > 1, args[2], "CosMxQC_report.Rmd")

output_dir <- dirname(output_file)
if(!dir.exists(output_dir) && nchar(output_dir) > 0) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# render
rmarkdown::render(
  input = input_file,
  output_file = output_file,
  output_format = "html_document"
)

message("\n\033[32mRender success:\033[0m ", normalizePath(output_file))
