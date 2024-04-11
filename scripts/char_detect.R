#!/usr/bin/env Rscript
# char_detect.R
# example:
# Rscript path-to/char_detect.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -o path-to/

# Set CRAN mirror for the script
options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

# Function to install a package if it's not already installed
install_if_not_installed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Check if required R packages are installed and prompt the user to install if not
required_packages <- c("optparse", "reshape2")

# Function to check and install packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Package", pkg, "is not installed.\n")
      cat("Please run the following command to install it:\n")
      cat("install.packages(\"", pkg, "\")\n")
      # Install and load the required packages
      lapply(required_packages, install_if_not_installed)
      stop("Exiting script. Please install the required packages and re-run the script.\n", call.=FALSE)
    }
  }
}

# Run the function to check for required packages
check_and_install_packages(required_packages)

# Now load the required R packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(reshape2))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript char_detect.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/char_detect.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -o path-to/\n")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()


# Define command-line option specifications
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input file, including the expression matrix and metadata."),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to the metadata file of the input expression. Note: The metadata file must include three columns: bacth, sample and type. Required!"),
  make_option(c("-o", "--out_dir"), type="character", default="./",
              help="Path to the output directory [default current directory]."),
  make_option(c("-h", "--help"), action="store_true", default=FALSE,
              help="Show this help message and exit")
)

# Define OptionParser object
opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)

# Parse command-line options
opt <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE))

# Show help message and exit
if (opt$help) {
  print_help(opt_parser)
  quit(save="no", status=0)
}

# Check for required command-line options
if (is.null(opt$input) || is.null(opt$metadata) || is.null(opt$out_dir)) {
  print_help(opt_parser)
  stop("Input file, metadata file, and output directory parameters must be provided.", call.=FALSE)
}

# Read input files and metadata
exprMat_count <- read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(exprMat_count)>0){
  message("Input expression profile is loaded.")
  }else{
stop("ERROR: input expression profile is not loaded!")
}

meta <- read.csv(opt$metadata,header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Process metadata
meta$batch_sample <- paste(meta$batch, meta$sample, sep="_")
meta_1 <- meta

# Process expression matrix
exprMat_count_M <- melt(as.matrix(exprMat_count))
colnames(exprMat_count_M)[1:3] <- c("gene", "library", "count")

# Merge metadata
exprMat_count_M$batch_sample <- meta_1$batch_sample[match(exprMat_count_M$library, meta_1$library)]

# Filter genes with counts greater than or equal to 3
expd_M <- exprMat_count_M[which(exprMat_count_M$count >= 3),]
expd_M <- expd_M[expd_M$library %in% meta_1$library,]

# Calculate the number of genes for each sample
expd_pg <- data.frame(tapply(expd_M$gene, as.factor(paste(expd_M$gene, expd_M$batch_sample)), length))

# Process output data frame
colnames(expd_pg) <- c("Num")
expd_pg$gene <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[1])
expd_pg$batch_sample <- sapply(strsplit(as.character(rownames(expd_pg)), " "), function(x) x[2])
expd_pg$batch <- meta_1$batch[match(expd_pg$batch_sample, meta_1$batch_sample)]
expd_pg$sample <- meta_1$sample[match(expd_pg$batch_sample, meta_1$batch_sample)]

# Filter genes detected in multiple samples
expd_pg_f <- expd_pg[expd_pg$Num > 1,]

# # Save results to csv file
# output_file <- paste0(opt$out_dir, "/detect_genelist_acrossSample_", Sys.Date(), ".csv")
# write.csv(expd_pg_f, output_file)

# Ensure metadata columns are character type
for (i in 1:ncol(meta)) {
  meta[,i] <- as.character(meta[,i])
}

# Read the list of detected genes
expdgene1 <- expd_pg_f

# Sample combinations
sample_combn <- data.frame(
  sampleA=rep("D6", 3),
  sampleB=c("D5", "F7", "M8")
)
sample_combn$compare <- paste(sample_combn$sampleB, "/", sample_combn$sampleA, sep="")

# Detect genes in all batches
expdgene1_p <- expdgene1

# Process samples
usample <- unique(as.character(meta$sample))

expd_refs <- c()
for (i in 1:length(usample)) {
  x <- expdgene1_p[expdgene1_p$sample == usample[i],]
  expd_sample <- data.frame(table(x$gene))
  expd_sample <- cbind(expd_sample, usample[i])
  expd_refs <- rbind(expd_refs, expd_sample)
}

expd_refs <- data.frame(expd_refs)
colnames(expd_refs) <- c("gene", "freq", "sample")
expd_refs$gene <- as.character(expd_refs$gene)
expd_refs$sample <- as.character(expd_refs$sample)

# Filter genes detected in all batches
expd_refs_f <- expd_refs[expd_refs$freq == length(unique(meta$batch)),]

# Save results to csv file
output_file2 <- paste0(opt$out_dir, "/detect_genelist_acrossBatch_", Sys.Date(), ".csv")
write.csv(expd_refs_f, output_file2)

# # Calculate the number of genes detected in each sample
# alldet <- names(which(table(expd_refs_f$gene) == 4))

# # Output results to CSV file
# output_csv <- paste0(opt$out_dir, "/char_detect_summary_", Sys.Date(), ".csv")
# write.csv(data.frame(sample = unique(expd_refs_f$sample),
#                      percent = tapply(expd_refs_f$gene, expd_refs_f$sample, length) / nrow(exprMat_count) * 100), 
#           file = output_csv, row.names = FALSE)

# Print completion message
message("Analysis complete! Results saved to: ")
# message(paste("-detect_genelist_acrossSample.csv: ", output_file))
message(paste("-detect_genelist_acrossBatch.csv: ", output_file2))
# message(paste("-char_detect_summary.csv: ", output_csv))