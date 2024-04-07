#!/usr/bin/env Rscript
# deg_limma.R
# example:
# Rscript path-to/deg_limma.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -o path-to/

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
required_packages <- c("optparse", "limma", "edgeR", "data.table")

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
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript deg_limma.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/deg_limma.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -o path-to/\n")
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

# Define DEG function
DEGanalysis <- function(exprMat, group){
  
  library(edgeR)
  library(limma)
  
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)
  
  keep <- filterByExpr(dge, design,min.count = 0)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  
  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$gene = rownames(result)
  result$groupA =  levels(group)[1]
  result$groupB =  levels(group)[2]
  return(as.data.table(result))
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
usample<-c("D5","D6","F7","M8")
sample_combn<-data.frame(
  sampleA=rep("D6",3),
  sampleB= c("D5","F7","M8"))

colnames(sample_combn)<-c("sampleA","sampleB")
sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")

ubatch<-sort(unique(as.character(meta$batch)))
for( i in 1:ncol(sample_combn)){
  sample_combn[,i]<-as.character(sample_combn[,i])
}

# DEG for D5 vs D6
DEG.cal.D56<-c()
for ( i in 1:length(ubatch)){
  meta_dat<-meta[meta$batch==ubatch[i],]
  expr<-exprMat_count[,match(meta_dat$library,colnames(exprMat_count))]
  
  j=1
    lst.library.sampleA=meta_dat$library[meta_dat$sample==sample_combn$sampleA[j]]
    lst.library.sampleB=meta_dat$library[meta_dat$sample==sample_combn$sampleB[j]]
    lst.library.forDEG <- c(lst.library.sampleA, lst.library.sampleB)
    
    y.DEG <- DEGanalysis(expr[,lst.library.forDEG],
                         factor(as.character(meta_dat$sample[match(lst.library.forDEG,meta_dat$library)])))
    
    y.DEG$batch<-ubatch[i]
    y.DEG<-data.frame(y.DEG)
    y.DEG$compare<-paste(y.DEG$groupB,"/",y.DEG$groupA,sep="")
    
    y.DEG_rev<-data.frame(
      logFC= -y.DEG$logFC,
      AveExpr=y.DEG$AveExpr,
      t=y.DEG$t,
      P.Value=y.DEG$P.Value,
      adj.P.Val=y.DEG$adj.P.Val,
      B=y.DEG$B,
      gene=y.DEG$gene,
      groupA=y.DEG$groupB,
      groupB=y.DEG$groupA,
      # When limma does a variance analysis, the name of the numerator group of the FC must be sorted after the character order of the name of the denominator group.
      batch=y.DEG$batch,
      compare=rep(sample_combn$compare[j],nrow(y.DEG))
    )
    
    DEG.cal.D56<-rbind(DEG.cal.D56,y.DEG_rev)
  
  print(i)
}

DEG.cal.D56<-data.frame(DEG.cal.D56)
DEG.cal.D56$Type<-"non-DEG"
DEG.cal.D56$Type[intersect(which(DEG.cal.D56$P.Value<0.05),which(DEG.cal.D56$logFC>=1))]<-"up-regulate"
DEG.cal.D56$Type[intersect(which(DEG.cal.D56$P.Value<0.05),which(DEG.cal.D56$logFC<=(-1)))]<-"down-regulate"

# DEG for F7/D6, M8/D6
DEG.cal.other<-c()
for ( i in 1:length(ubatch)){
  meta_dat<-meta[meta$batch==ubatch[i],]
  expr<-exprMat_count[,match(meta_dat$library,colnames(exprMat_count))]
  
  for (j in 2:nrow(sample_combn)){
    lst.library.sampleA=meta_dat$library[meta_dat$sample==sample_combn$sampleA[j]]
    lst.library.sampleB=meta_dat$library[meta_dat$sample==sample_combn$sampleB[j]]
    lst.library.forDEG <- c(lst.library.sampleA, lst.library.sampleB)
     
    y.DEG <- DEGanalysis(expr[,lst.library.forDEG],
                         factor(as.character(meta_dat$sample[match(lst.library.forDEG,meta_dat$library)])))
    
    y.DEG$batch<-ubatch[i]
    y.DEG$compare<-sample_combn$compare[j]
    
    DEG.cal.other<-rbind(DEG.cal.other,y.DEG)
  }
  print(i)
}

DEG.cal.other<-data.frame(DEG.cal.other)
DEG.cal.other$Type<-"non-DEG"
DEG.cal.other$Type[intersect(which(DEG.cal.other$P.Value<0.05),which(DEG.cal.other$logFC>=1))]<-"up-regulate"
DEG.cal.other$Type[intersect(which(DEG.cal.other$P.Value<0.05),which(DEG.cal.other$logFC<=(-1)))]<-"down-regulate"

DEG.cal3<-rbind(DEG.cal.D56,DEG.cal.other)

# Save results to csv file
output_csv <- paste0(opt$out_dir, "/DEG_limma_", Sys.Date(), ".csv")
fwrite(DEG.cal3, output_csv,row.names = TRUE)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-DEG_limma.csv: ", output_csv))
#message(paste("-deg_limma_summary.csv: ", output_csv))
