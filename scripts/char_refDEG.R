#!/usr/bin/env Rscript
# char_refDEG.R
# example:
# Rscript path-to/char_refDEG.R -d path-to/DEG_limma_2024-04-02.csv -r ref_expr_2024-04-02.csv -o path-to/

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
  cat("Actual Usage: Rscript char_refDEG.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/char_refDEG.R -d path-to/DEG_limma_(date).csv -r path-to/ref_expr_(date).csv -o path-to/\n")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()


# Define command-line option specifications
option_list <- list(
  make_option(c("-d", "--deg"), type="character", default=NULL,
              help="Path to the deg file from the output of deg_limma.R."),
  make_option(c("-r", "--refFC"), type="character", default=NULL,
              help="Path to the ref expr file from the output of ratio_expr.R."),
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
if (is.null(opt$deg) || is.null(opt$refFC) || is.null(opt$out_dir)) {
  print_help(opt_parser)
  stop("DEG file from the output of deg_limma.R, and output directory parameters must be provided.", call.=FALSE)
}

# Read input files
DEGs <- read.csv(opt$deg,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(DEGs)>0){
  message("Input DEG is loaded.")
  }else{
stop("ERROR: input DEG is not loaded!")
}

## preprocess
DEGs_p<-DEGs
DEGs_p$protocol<-sapply(strsplit(as.character(DEGs_p$batch),"_"),function(x){x[1]})
DEGs_p$gene_compare<-paste(DEGs_p$gene,DEGs_p$compare)

DEGs_p$DEG_type<-"non-DEG"
DEGs_p$DEG_type[intersect(which(DEGs_p$P.Value<0.05),which(DEGs_p$logFC>=1))]<-"up-regulate"
DEGs_p$DEG_type[intersect(which(DEGs_p$P.Value<0.05),which(DEGs_p$logFC<=(-1)))]<-"down-regulate"

########DEGs_cal 
DEG_ref_cal<-data.frame(
  N_up=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="up-regulate"))}),
  N_non=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="non-DEG"))}),
  N_down=tapply(DEGs_p$DEG_type,as.factor(DEGs_p$gene_compare),function(x){length(which(x=="down-regulate"))}))

DEG_ref_cal$Final<-"non-DEG"
DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_up>=2),which(DEG_ref_cal$N_down>=2))]<-"conflicting"

DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_up>=4),which(DEG_ref_cal$N_down==0))]<-"up-regulate"
DEG_ref_cal$Final[intersect(which(DEG_ref_cal$N_down>=4),which(DEG_ref_cal$N_up==0))]<-"down-regulate"
DEG_ref_cal$gene<-sapply(strsplit(rownames(DEG_ref_cal)," "),function(x){x[1]})
DEG_ref_cal$compare<-sapply(strsplit(rownames(DEG_ref_cal)," "),function(x){x[2]})

DEG_ref_cal$Var1<-rownames(DEG_ref_cal)


######ref_FC
refFC_final <- read.csv(opt$refFC,header=T,stringsAsFactors=F,row.names=1,check.names=F)
DEG_ref_cal_f<-DEG_ref_cal[DEG_ref_cal$Var1 %in% refFC_final$gene_compare, ]

refFC_final$N_up<-DEG_ref_cal$N_up[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$N_non<-DEG_ref_cal$N_non[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$N_down<-DEG_ref_cal$N_down[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]
refFC_final$Final<-DEG_ref_cal$Final[match(refFC_final$gene_compare,DEG_ref_cal$Var1)]

table(refFC_final$Final,refFC_final$compare)


DEG_ref<-refFC_final[grep("regulate",refFC_final$Final),]

# Ensure column names are lowercase with underscores
colnames(refFC_final) <- tolower(gsub(".", "_", colnames(refFC_final), fixed = TRUE))


# Output results to CSV file
output_file2 <- paste0(opt$out_dir, "/RefData_DEGs_", Sys.Date(), ".csv")
write.csv(refFC_final, output_file2)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-RefData_DEGs.csv: ", output_file2))