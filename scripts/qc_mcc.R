# MCC based on both detected DEG and nonDEG,
#!/usr/bin/env Rscript
# example:
# Rscript path-to/qc_mcc.R -i path-to/DEG_limma_(date).csv -r path-to/example_ref_data.csv -o path-to/

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
  cat("Actual Usage: Rscript qc_mcc.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/qc_mcc.R -i path-to/DEG_limma_(date).csv -r path-to/example_ref_data.csv -o path-to/\n")
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
              help="File name of Fold change file in log2 scale. Required!"),
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

##import file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

Ref_RE<-read.csv(opt$refdata)
logFC_raw<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)
ubatch<-unique(as.character(logFC_raw$batch))
#logexpr
if(nrow(logFC_raw)>0){
  message("Input expression profile is loaded.")
}else{
  stop("ERROR: input expression profile is not loaded!")
}
#--------
#logFC_raw<-readRDS("../RNA_6/refdata_biaowu_20220426/expr_mat_refdata/RNA_DEGanalysis_limma_20220427.rds")
logFC_raw$gene_compare<-paste(logFC_raw$gene,logFC_raw$compare)

#ref_DEG<-readRDS("../RNA_6/refdata_biaowu_20220426/expr_mat_refdata/RefData_DEGs_r5036_20220721.rds")
#ref_nonDEG<-readRDS("../RNA_6/refdata_biaowu_20220426/expr_mat_refdata/RefData_nonDEG_20220724.rds")
unique(Ref_RE[,"Final"])
ref_nonDEG<-Ref_RE[which(Ref_RE$Final=="non-DEG"),]
ref_DEG<-Ref_RE[-which(Ref_RE$Final=="non-DEG"),]

sample_combn<- unique(ref_DEG$compare)
sample_combn<-as.data.frame(sample_combn)
colnames(sample_combn)<-"compare"
DEG_QC<-c()

#########refs
for ( i in 1:nrow(sample_combn)){
  ref_DEG_com<-ref_DEG[ref_DEG$compare==sample_combn$compare[i],]
  ref_nonDEG_com<-ref_nonDEG[ref_nonDEG$compare==sample_combn$compare[i],]
  
  ref_nd<-as.character(ref_DEG_com[which(ref_DEG_com$Final=="down-regulate"),"gene_compare"])
  ref_nu<-as.character(ref_DEG_com[which(ref_DEG_com$Final=="up-regulate"),"gene_compare"])
  ref_no<-as.character(ref_nonDEG_com$gene_compare)
  refgenes<-c(ref_nd,ref_nu,ref_no)
  
  #########
  logFC_raw_pairs<-logFC_raw[logFC_raw$gene_compare %in% refgenes,]
  
  ###########batch 
  for (k in 1:length(ubatch)){
    test_lab<-logFC_raw_pairs[logFC_raw_pairs$batch==ubatch[k],]
    
    bothgene<-intersect(test_lab$gene_compare,refgenes)
    
    ref_no_f<-intersect(ref_no,bothgene)
    ref_nu_f<-intersect(ref_nu,bothgene)
    ref_nd_f<-intersect(ref_nd,bothgene)
    
    test_lab_non<-test_lab$gene_compare[test_lab$type=="non-DEG"]
    test_lab_up<-test_lab$gene_compare[test_lab$type=="up-regulate"]
    test_lab_down<-test_lab$gene_compare[test_lab$type=="down-regulate"]
    
    test_lab_non<-intersect(test_lab_non,bothgene)
    test_lab_up<-intersect(test_lab_up,bothgene)
    test_lab_down<-intersect(test_lab_down,bothgene)
    
    
    TP=as.numeric(length(intersect(ref_nd_f,test_lab_down)))+as.numeric(length(intersect(ref_nu_f,test_lab_up)))
    TN=as.numeric(length(intersect(ref_no_f,test_lab_non)))
    FN=as.numeric(length(setdiff(test_lab_non,ref_no_f)))
    FP=as.numeric(length(intersect(c(test_lab_down,test_lab_up),ref_no_f)))
    
    Precision=TP/(TP+FP)
    Sensitivity=Recall=TPR=TP/(TP+FN)
    Specificity=TNR=TN/(FP+TN)
    F1= 2*(Recall * Precision) / (Recall + Precision)
    MCC= (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    
    DEG_QC<-rbind(DEG_QC,c(TP,TN,FN,FP,Precision,Sensitivity,Specificity,F1,MCC,ubatch[k],sample_combn$compare[i]))
    
  }
  
}

colnames(DEG_QC)<-c("tp","tn","fn","fp","precision","sensitivity","specificity","f1","mcc","batch","compare")

DEG_QC<-data.frame(DEG_QC)
for ( i in 1:9){
  DEG_QC[,i]<-as.numeric(as.character(DEG_QC[,i]))
}

DEG_QC21<-apply(DEG_QC[,1:9],2,function(x){tapply(x,DEG_QC$batch,mean)})

DEG_QC21<-data.frame(DEG_QC21)

DEG_QC21$batch<-rownames(DEG_QC21)

##output

MCCfile<-paste(out_dir,"MCC_each_batch_",gsub("-","",Sys.Date()),".csv",sep="")
write.csv(DEG_QC21,MCCfile)

message("Analysis complete! Results saved to:")
message(paste("-MCC_each_batch.csv:",MCCfile,sep=" "))
