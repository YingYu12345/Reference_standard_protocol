# Root Mean Square Error (RMSE) is the standard deviation of the residuals (prediction errors). 
#!/usr/bin/env Rscript
# example:
# Rscript path-to/qc_rmse.R -i path-to/DEG_limma_(date).csv -r path-to/example_ref_data.csv -o path-to/

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
required_packages <- c("optparse", "Metrics")

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
suppressPackageStartupMessages(library("Metrics"))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript qc_rmse.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/qc_rmse.R -i path-to/DEG_limma_(date).csv -r path-to/example_ref_data.csv -o path-to/\n")
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

##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

##import exp file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

#metadata<-read.csv(opt$metadata,header=T,stringsAsFactors=F,row.names=1,check.names=F)
Ref_RE<-read.csv(opt$refdata)

logFC_raw<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)
ubatch<-unique(as.character(logFC_raw$batch))
#logexpr
if(nrow(logFC_raw)>0){
  message("Input expression profile is loaded.")
}else{
  stop("ERROR: input expression profile is not loaded!")
}

#for ( i in 1:ncol(meta)){
#  meta[,i]<-as.character(meta[,i])
#}


log2FC_QC_rmse<-c()  
Ref_RE$gene_group<-paste(Ref_RE$gene,Ref_RE$compare)

for ( i in 1:length(ubatch)){
  testlab_logFC<-logFC_raw[logFC_raw$batch==ubatch[i],]
  testlab_logFC$gene_group<-paste(testlab_logFC$gene,testlab_logFC$compare)
  int<-intersect(Ref_RE$gene_group,testlab_logFC$gene_group)
  testlab_logFC_f<-testlab_logFC[match(int,testlab_logFC$gene_group),]
  Ref_RE_f<-Ref_RE[match(int,Ref_RE$gene_group),]
  
  colnames(Ref_RE_f)<-paste("ref_",colnames(Ref_RE_f),sep="")
  colnames(testlab_logFC_f)<-paste("testlab_",colnames(testlab_logFC_f),sep="")
  refvstest<-merge(Ref_RE_f[,c("ref_gene_group","ref_compare","ref_meanlogFC")],testlab_logFC_f[,c("testlab_gene_group","testlab_logfc")],by.x="ref_gene_group",by.y="testlab_gene_group")
  refvstest<-refvstest[apply(refvstest,1,function(x){length(which(is.na(x)))==0}),]
  
  rmse_value<-rmse(refvstest$testlab_logfc,refvstest$ref_meanlogFC)
  
  log2FC_QC_rmse<-rbind(log2FC_QC_rmse,c(rmse_value,ubatch[i]))
  
  print(i)
}

log2FC_QC_rmse<-data.frame(log2FC_QC_rmse)
colnames(log2FC_QC_rmse)<-c("rmse","batch")
log2FC_QC_rmse$rmse<-as.numeric(as.character(log2FC_QC_rmse$rmse))

#DEG_QC21$RMSE<-log2FC_QC_rmse$RMSE[match(rownames(DEG_QC21),log2FC_QC_rmse$batch)]

##output
RMSEfile<-paste(out_dir,"RMSE_each_batch_",gsub("-","",Sys.Date()),".csv",sep="")
write.csv(log2FC_QC_rmse,RMSEfile)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-RMSE_each_batch.csv: ", RMSEfile))