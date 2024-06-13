#!/usr/bin/env Rscript
# Example usage:
# Rscript path-to/uchar.R -i path-to/DEG_limma_(date).csv -p path-to/passbatch.csv -r path-to/ref_expr_(date).csv -o path-to/

suppressPackageStartupMessages(library("optparse"))

# Define command-line options
option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character", default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-p", "--passbatch"), type="character", default=NULL,
              help="List of passbatch. Required!"),
  make_option(c("-r", "--refdata"), type="character", default=NULL,
              help="File name of reference data in ratio-scale. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit")
)

# Parse command-line options
opt <- parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

# Pre-analysis checks
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

# Set output directory
out_dir <- paste(gsub("/$", "", opt$out_dir), "/", sep="")

# Load expression profile
DEGs <- read.csv(opt$input, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
if (nrow(DEGs) > 0) {
  message("Input expression profile is loaded.")
} else {
  stop("ERROR: input expression profile is not loaded!")
}

# Load list of passed batches
passbatch <- read.csv(opt$passbatch, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
passbatch <- as.character(passbatch[,1])
if (length(passbatch) > 0) {
  message("List of passed batches is loaded.")
} else {
  stop("ERROR: list of passed batches is not loaded!")
}

# Load reference dataset
refdata <- read.csv(opt$refdata, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(refdata) > 0) {
  message("Input reference dataset is loaded.")
} else {
  stop("ERROR: input reference dataset is not loaded!")
}

# Filter DEGs for passed batches
DEGs_p <- DEGs[DEGs$batch %in% passbatch,]
DEGs_p$gene_compare <- paste(DEGs_p$gene, DEGs_p$compare)

# Calculate UCHAR
rse <- function(x) {
  sd(x, na.rm=TRUE) / (mean(x, na.rm=TRUE) * sqrt(length(x[!is.na(x)])))
}
gene_compare <- levels(as.factor(DEGs_p$gene_compare))

uchar_mat <- data.frame(
  uchar = tapply(DEGs_p$logfc, as.factor(as.character(DEGs_p$gene_compare)), function(x) { rse(2^x) * 100 }),
  gene_compare = gene_compare,
  gene = sapply(strsplit(gene_compare, " "), function(x) { x[1] }),
  compare = sapply(strsplit(gene_compare, " "), function(x) { x[2] })
)

refdata_uchar <- refdata
refdata_uchar$uchar <- uchar_mat$uchar[match(refdata_uchar$gene_compare, uchar_mat$gene_compare)]

# Prepare output data
refdata_uchar_out <- refdata_uchar[, c("gene_compare", "gene", "compare", "uchar")]

ucharfile <- paste(out_dir, "Refdata_uchar_", Sys.Date(), ".csv", sep="")

# Print summary
print("Summary")
print(tapply(refdata_uchar$uchar, refdata_uchar$compare, summary))

# Save output
write.csv(refdata_uchar_out, ucharfile, row.names=FALSE)

message("Finished!")
message(paste("Output file: \n", ucharfile, sep=" "))