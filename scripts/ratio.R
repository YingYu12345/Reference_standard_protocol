#!/usr/bin/env Rscript
# ratio.R
# example:
# Rscript path-to/ratio.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/

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
required_packages <- c("optparse")

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

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript ratio.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/ratio.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/\n")
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
logFPKM <- read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(logFPKM)>0){
  message("Input expression profile is loaded.")
  }else{
stop("ERROR: input expression profile is not loaded!")
}

meta <- read.csv(opt$metadata,header=T,stringsAsFactors=F,row.names=1,check.names=F)

# Process metadata
meta$rep<-sapply(strsplit(as.character(meta$library),"_"),function(x){x[6]})
meta$samplerep<-paste(meta$sample,meta$rep,sep="_")
meta$library<-as.character(meta$library)
meta$sample<-as.character(meta$sample)

ubatch<-unique(as.character(meta$batch))

# Process expression matrix
ratio_D6<-matrix(0,ncol=ncol(logFPKM),nrow=nrow(logFPKM))
rownames(ratio_D6)<-rownames(logFPKM)
colnames(ratio_D6)<-colnames(logFPKM)  

for ( i in 1:length(ubatch)){
  name<-as.character(meta$library[meta$batch==ubatch[i]])
  logexpr_batch<-logFPKM[,colnames(logFPKM)%in% name]
  m<-rowMeans(logexpr_batch[,grep("D6",colnames(logexpr_batch))])
  
  mat<-apply(logexpr_batch,2,function(x){x-m})
  
  ratio_D6[rownames(mat),colnames(mat)]<-mat
  print(i)
}


# Output results to CSV file
output_csv <- paste0(opt$out_dir, "/RNA_ratio3D6_", Sys.Date(), ".csv")
write.csv(ratio_D6, output_csv)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-RNA_ratio3D6_.csv: ", output_csv))