#!/usr/bin/env Rscript
# example:
# Rscript path-to/homogeneity_assess.R -i path-to/example_expr_homo_log2.csv -m metadata_homo.csv -o path-to/

suppressPackageStartupMessages(library("optparse"))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript homogeneity_assess.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/homogeneity_assess.R -i path-to/example_expr_homo_log2.csv -m metadata_homo.csv -o path-to/\n")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
option_list <- list(
  make_option(c("-o", "--out_dir"), type="character", default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Metadata of the input expression file. Note: metadata file shall include two columns: sample and type (intra- or between-unit). Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE,
              help="Show this help message and exit")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))

#test
#opt$input<-"expr_mat/example_expr_homo_log2.csv"
#opt$metadata<-"expr_mat/metadata_homo.csv"


# Pre-analysis
if (is.null(opt$input)) {
  print_help(OptionParser(option_list=option_list))
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

# Import expression file
out_dir <- paste0(gsub("/$", "", opt$out_dir), "/")
logexpr <- read.csv(opt$input, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

if (nrow(logexpr) > 0) {
  message("Input expression profile is loaded.")
} else {
  stop("ERROR: input expression profile is not loaded!")
}

metadata <- read.csv(opt$metadata, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
usample <- unique(metadata$sample)

if (nrow(metadata) > 0) {
  message("Metadata is loaded.")
} else {
  stop("ERROR: metadata is not loaded!")
}

# Analysis
homoge_p <- data.frame()
homoge <- data.frame()

for (j in seq_along(usample)) {
  for (i in seq_len(nrow(logexpr))) {
    g <- rownames(logexpr)[i]
    df <- metadata[metadata$sample == usample[j],]
    df$value <- as.numeric(logexpr[g, match(rownames(df), colnames(logexpr))])
    
    p <- summary(aov(df$value ~ df$type))[[1]][["Pr(>F)"]][1]
    means <- tapply(df$value, df$type, mean)
    
    homoge_p <- rbind(homoge_p, c(usample[j], g, p, means))
  }
  
  colnames(homoge_p) <- c("sample", "gene", "anova_p", "mean_between", "mean_intra")
  homoge_p <- as.data.frame(homoge_p)
  homoge_p$gene <- as.character(homoge_p$gene)
  
  homoge_p[, 3:ncol(homoge_p)] <- lapply(homoge_p[, 3:ncol(homoge_p)], function(x) round(as.numeric(as.character(x)), 4))
  
  homoge_p$fdr <- round(p.adjust(homoge_p$anova_p, method="fdr"), 4)
  homoge <- rbind(homoge, homoge_p)
  
  homoge_p <- data.frame()
}

sums <- data.frame(
  sample=levels(as.factor(homoge$sample)),
  length_sig=tapply(homoge$fdr, homoge$sample, function(x) length(which(x < 0.05))),
  length_total=tapply(homoge$fdr, homoge$sample, length),
  length_ratio=tapply(homoge$fdr, homoge$sample, function(x) length(which(x < 0.05)) / length(x))
)
rownames(sums) <- seq_len(nrow(sums))

# Output
homogefile <- paste0(out_dir, "homogeneity_Ftest_detail_", Sys.Date(), ".csv")
sumfile <- paste0(out_dir, "homogeneity_Ftest_summary_", Sys.Date(), ".csv")

write.csv(homoge, homogefile, row.names=FALSE)
write.csv(sums, sumfile, row.names=FALSE)

message("Finished!")
message(paste("Output file: \n", homogefile, "\nand\n", sumfile, sep=" "))