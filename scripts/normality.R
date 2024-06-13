#!/usr/bin/env Rscript
# Example usage:
# Rscript path-to/normality.R -i path-to/DEG_limma_(date).csv -p path-to/batch_included.csv -o path-to/

# Load required libraries
suppressPackageStartupMessages(library("optparse"))

# Specify command line options
option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character", default="./",
              help="The output directory [default: ./]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-p", "--passbatch"), type="character", default=NULL,
              help="List of passed batches. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE))

# Show help message and exit if required options are not provided
if (is.null(opt$input) || is.null(opt$passbatch)) {
  print_help(opt_parser)
  stop("Input file and list of passed batches must be supplied.", call.=FALSE)
}

# Define output directory
out_dir <- file.path(gsub("/$", "", opt$out_dir), "/")

# Load input files
DEGs <- read.csv(opt$input, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
if (nrow(DEGs) == 0) {
  stop("ERROR: input expression profile is not loaded!")
} else {
  message("Input expression profile is loaded.")
}

passbatch <- read.csv(opt$passbatch, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
passbatch <- as.character(passbatch[,1])
if (length(passbatch) == 0) {
  stop("ERROR: list of passed batches is not loaded!")
} else {
  message("List of passed batches is loaded.")
}

# Filter DEGs for passed batches and create gene comparison column
DEGs_filtered <- DEGs[DEGs$batch %in% passbatch, ]
DEGs_filtered$gene_compare <- paste(DEGs_filtered$gene, DEGs_filtered$compare)

# Initialize normality matrix
normality_results <- data.frame()

# Perform Shapiro-Wilk normality test for each gene comparison
unique_gene_comparisons <- unique(DEGs_filtered$gene_compare)
for (i in seq_along(unique_gene_comparisons)) {
  gene_compare <- unique_gene_comparisons[i]
  subset_data <- DEGs_filtered[DEGs_filtered$gene_compare == gene_compare, ]
  if (nrow(subset_data) > 3) {
    p_value <- shapiro.test(subset_data$logfc)$p.value
    normality_results <- rbind(normality_results, c(gene_compare, p_value))
  }
  if (i %% 1000 == 0) {
    print(paste(i, "/", length(unique_gene_comparisons)))
  }
}

# Process normality results
colnames(normality_results) <- c("gene_compare", "log2FC_p")
normality_results$log2FC_p <- as.numeric(as.character(normality_results$log2FC_p))
normality_results$gene <- sapply(strsplit(as.character(normality_results$gene_compare), " "), function(x) x[1])
normality_results$compare <- sapply(strsplit(as.character(normality_results$gene_compare), " "), function(x) x[2])

# Calculate and print the percentage of genes passing the Shapiro-Wilk test
passed_percentage <- round(length(which(normality_results$log2FC_p > 0.05)) / nrow(normality_results) * 100, 2)
print(paste(passed_percentage, "% of genes passed the Shapiro-Wilk test (p > 0.05). The expression profiles across batches are normally distributed."))

# Ensure column names are lowercase with underscores
colnames(normality_results) <- tolower(colnames(normality_results))

# Save results to a CSV file
output_file <- paste0(out_dir, "normality_genelist_", Sys.Date(), ".csv")
write.csv(normality_results, output_file, row.names=FALSE)

# Print completion message
message("Finished!")
message(paste("Output file: \n", output_file))