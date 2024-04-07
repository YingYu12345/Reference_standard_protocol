#!/usr/bin/env Rscript
# pvca.R
# example:
# Rscript path-to/pvca.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/

# pre-install
# conda create -n pvca_env r bioconductor-pvca=1.42.0
# conda activate pvca_env
# Rscript pvca.R -i example_expr_multibatch_log2.csv -m metadata.csv -o .

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
suppressPackageStartupMessages(library(pvca))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(ggplot2))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript pvca.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/pvca.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/\n")
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
meta$library<-NULL

row.names(meta)<-paste0(meta$batch,"_",meta$sample_rep)
meta$batch<-NULL
meta$sample_rep<-NULL

# PVCA
RNA_data<-ExpressionSet(as.matrix(logFPKM),
                        phenoData=AnnotatedDataFrame(meta))

pct_threshold <- 0.6
batch.factors <- c("sample","lab", "platform", "protocol","rep" )

pvcaObj <- pvca::pvcaBatchAssess(RNA_data, batch.factors, pct_threshold)
RNA_bats<-data.frame(value=as.numeric(pvcaObj$dat),
                     name=pvcaObj$label)
RNA_bats$value2<-round(RNA_bats$value,3)
RNA_bats_f<-RNA_bats[RNA_bats$value2>=0.001,]
r1 <- RNA_bats[which(RNA_bats$name=="platform:protocol"),]
RNA_bats_f<-rbind(RNA_bats_f,r1)
RNA_bats_f$type<-ifelse(RNA_bats_f$name=="sample","biological",ifelse(grepl("sample",RNA_bats_f$name),"mixed","techinical"))
RNA_bats_f$type[RNA_bats_f$name=="resid"]<-"other"

# Output results to CSV file
output_csv <- paste0(opt$out_dir, "/RNA_PVCA_", Sys.Date(), ".csv")
write.csv(RNA_bats_f, output_csv)

# prepare for plot
RNA_bats_f$type <- factor(RNA_bats_f$type,levels=c('biological','techinical','mixed','other'))
# plot
p1 <- ggplot(RNA_bats_f,aes(x=reorder(name,-value),y=value))+
  geom_bar(aes(fill=type),stat = 'identity')+
  geom_text(aes(label=sprintf("%4.3f",as.numeric(value))),size=2.5,hjust=0.5,vjust=-0.5)+
  labs(x = "Effect",y = "Weighted average proportion variance")+
  scale_fill_manual(values = c("#CD534CFF","#DF8F44FF","#5A9599FF","#868686FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color='black'),
        axis.text.x = element_text(color='black', angle = 30, hjust = 1),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_blank(),
        legend.text = element_text(color='black'),
        legend.position = c(0.5,0.98),
        legend.direction = "horizontal",
        legend.key.size = unit(0.3,'cm'),
        legend.margin = margin(0,0,0,0));p1

# Output results to pdf file
output_pdf <- paste0(opt$out_dir, "/RNA_PVCA_", Sys.Date(), ".pdf")  
ggsave(output_pdf, plot = p1,
       width = 5, height = 4, units = "in", useDingbats=FALSE)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-RNA_PVCA.csv: ", output_csv))
message(paste("-RNA_PVCA.pdf: ", output_pdf))
