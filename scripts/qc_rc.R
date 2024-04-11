# RC
#!/usr/bin/env Rscript
# example:
# Rscript path-to/qc_rc.R -i path-to/DEG_limma_2024-04-02.csv -r path-to/example_ref_data.csv -o path-to/

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

suppressPackageStartupMessages(library("optparse"))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript qc_rc.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/qc_rc.R -i path-to/DEG_limma_2024-04-02.csv -r path-to/example_ref_data.csv -o path-to/\n")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# Define command-line option specifications
option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"),type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-r", "--refdata"),type="character",  default=NULL,
              help="Reference datases. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 

##import file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

Ref_RE<-read.csv(opt$refdata)
Ref_RE$gene_group<-paste(Ref_RE$gene,Ref_RE$compare)

logFC_raw<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)
logFC_raw<-data.frame(logFC_raw)
logFC_raw$gene_compare<-paste(logFC_raw$gene,logFC_raw$compare)
ubatch<-unique(as.character(logFC_raw$batch))

log2FC_QC_cor<-c()  

for ( i in 1:length(ubatch)){
  testlab_logFC<-logFC_raw[logFC_raw$batch==ubatch[i],]
  testlab_logFC$gene_group<-paste(testlab_logFC$gene,testlab_logFC$compare)
  int<-intersect(Ref_RE$gene_group,testlab_logFC$gene_group)
  testlab_logFC_f<-testlab_logFC[match(int,testlab_logFC$gene_group),]
  Ref_RE_f<-Ref_RE[match(int,Ref_RE$gene_group),]
  
  colnames(Ref_RE_f)<-paste("ref_",colnames(Ref_RE_f),sep="")
  colnames(testlab_logFC_f)<-paste("testlab_",colnames(testlab_logFC_f),sep="")
  refvstest<-merge(Ref_RE_f[,c("ref_gene_group","ref_compare","ref_meanlogFC")],testlab_logFC_f[,c("testlab_gene_group","testlab_logFC")],by.x="ref_gene_group",by.y="testlab_gene_group")
  refvstest<-refvstest[apply(refvstest,1,function(x){length(which(is.na(x)))==0}),]
  
  c1<-cor(refvstest$ref_meanlogFC,refvstest$testlab_logFC)
  log2FC_QC_cor<-rbind(log2FC_QC_cor,c(c1,ubatch[i]))
  print(i)
}


log2FC_QC_cor<-data.frame(log2FC_QC_cor)
colnames(log2FC_QC_cor)<-c("RC","batch")
log2FC_QC_cor$RC<-as.numeric(as.character(log2FC_QC_cor$RC))

##output
RCfile<-paste(out_dir,"RC_each_batch_",gsub("-","",Sys.Date()),".csv",sep="")
write.csv(log2FC_QC_cor,RCfile)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-RC_each_batch.csv: ", RCfile))