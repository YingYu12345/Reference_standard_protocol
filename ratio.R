#!/usr/bin/env Rscript
# ratio_expr.R
# example:
# Rscript path-to/ratio_expr.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -g path-to/detect_genelist_acrossBatch_2024-04-02.csv -d path-to/DEG_limma_2024-04-02.csv -o path-to/

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
  cat("Actual Usage: Rscript ratio_expr.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/ratio_expr.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -g path-to/detect_genelist_acrossBatch_2024-04-02.csv -d path-to/DEG_limma_2024-04-02.csv -o path-to/\n")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()


# Define command-line option specifications
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input file, including the expression matrix and metadata."),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to the metadata file of the input expression. Note: The metadata file must include three columns: bacth, sample and type. Required!"),
  make_option(c("-g", "--genelist"), type="character", default=NULL,
              help="Path to the gene list file from the output of char_detect.R."),
  make_option(c("-d", "--deg"), type="character", default=NULL,
              help="Path to the deg file from the output of deg_limma.R."),
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

expd_refs_f <- read.csv(opt$genelist,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(expd_refs_f)>0){
  message("Input detected gene list is loaded.")
  }else{
stop("ERROR: input detected gene list is not loaded!")
}

DEGs <- read.csv(opt$deg,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(DEGs)>0){
  message("Input DEG is loaded.")
  }else{
stop("ERROR: input DEG is not loaded!")
}

# metainfo 
meta <- read.csv(opt$metadata,header=T,stringsAsFactors=F,row.names=1,check.names=F)
sample_combn <- data.frame(
  sampleA=rep("D6",3),
  sampleB=c("D5","F7","M8")
)

sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")

# define sample pair
detect_genes_pairs<-c()
for (i in 1:nrow(sample_combn)){
  s<-intersect(expd_refs_f$gene[expd_refs_f$sample==sample_combn$sampleA[i]],expd_refs_f$gene[expd_refs_f$sample==sample_combn$sampleB[i]])
  detect_genes_pairs<-rbind(detect_genes_pairs,cbind(s,as.character(sample_combn$compare[i])))
}
colnames(detect_genes_pairs)<-c("gene","compare")
detect_genes_pairs<-data.frame(detect_genes_pairs)

detect_genes_pairs$gene_compare<-paste(detect_genes_pairs$gene,detect_genes_pairs$compare)

# Output results to CSV file
output_file2 <- paste0(opt$out_dir, "/detected_gene_pairs_", Sys.Date(), ".csv")
write.csv(detect_genes_pairs, output_file2)

# FC and P
DEGs_p <- DEGs
DEGs_p$protocol <- sapply(strsplit(as.character(DEGs_p$batch),"_"),function(x){x[1]})
DEGs_p$gene_compare <- paste(DEGs_p$gene,DEGs_p$compare)

DEGs_p_dec <- DEGs_p[DEGs_p$gene_compare %in%  detect_genes_pairs$gene_compare,]

# #filter p<0.05
DEGs_p_f <- DEGs_p_dec[DEGs_p_dec$`P.Value`<0.05,]

DEGs_p_cal <- data.frame(table(paste(DEGs_p_f$gene,DEGs_p_f$compare)))
DEGs_p_cal$gene <- sapply(strsplit(as.character(DEGs_p_cal$Var1)," "),function(x){x[1]})
DEGs_p_cal$compare <- sapply(strsplit(as.character(DEGs_p_cal$Var1)," "),function(x){x[2]})

DEGs_p_cal_f <- DEGs_p_cal[DEGs_p_cal$Freq>=4,]
DEGs_p_cal_f_fd <- DEGs_p_cal_f[DEGs_p_cal_f$Var1 %in% detect_genes_pairs$gene_compare,]

DEGs_p_cal_f_fd <- DEGs_p_cal_f_fd[order(DEGs_p_cal_f_fd$compare),]
DEGs_p_cal_f_fd$gene_compare <- paste(DEGs_p_cal_f_fd$gene,DEGs_p_cal_f_fd$compare)

##filtering protocol-dependent genes
### Deine the pvalue-function
pvalue <- function(x,group){
      obj <- try(t.test(x[group==levels(group)[1]], x[group==levels(group)[2]],var.equal = TRUE), silent=TRUE)
    if (is(obj, "try-error")){p.value =1}	else {p.value = obj$p.value}
 
    return(p.value)   
}

DEGs_p_cal_f_fd2 <- DEGs_p_cal_f_fd
prot_fc<-c()
for ( i in 1:nrow(DEGs_p_cal_f_fd2)){
  
  c=as.character(DEGs_p_cal_f_fd2$compare)[i]
  g=as.character(DEGs_p_cal_f_fd2$gene)[i]
  m<-DEGs_p[intersect(which(DEGs_p$gene==g),which(DEGs_p$compare==c)),]
  
  if(length(which(grepl("P",m$protocol)))>=3 && length(which(grepl("R",m$protocol)))>=3){
    p<-pvalue(m$logFC,as.factor(m$protocol))
    f<-tapply(m$logFC,as.factor(m$protocol),mean)
    fc<-f[1]-f[2]
    prot_fc<-rbind(prot_fc,c(g,c,p,f,fc))
  }
    
  if(i %in% seq(0,nrow(DEGs_p_cal_f_fd2),by=200)){
print(i)
  }
  
}

colnames(prot_fc)<-c("gene","compare","pvalue","mean_polyA","mean_riboZ","FC_P2R")
prot_fc<-data.frame(prot_fc)
for ( i in 1:ncol(prot_fc)){
  prot_fc[,i]<-as.character(prot_fc[,i])
  
}

for ( i in 3:ncol(prot_fc)){
  prot_fc[,i]<-as.numeric(prot_fc[,i])
  
}

prot_fc$gene_compare<-paste(prot_fc$gene,prot_fc$compare)
ff<-prot_fc[intersect(which(prot_fc$pvalue<0.05),which(abs(prot_fc$FC_P2R)>=1)),]

#####filter 
ref_FC_f2<-DEGs_p_cal_f_fd2[!as.character(DEGs_p_cal_f_fd2$gene_compare) %in% as.character(ff$gene_compare),]


###########################
####meanFC medianp,seFC
DEGs_p$gene_compare <- paste(DEGs_p$gene,DEGs_p$compare)

#function #
gene_compare=levels(as.factor(DEGs_p$gene_compare))

mm.DEG<-data.frame(
  meanlogFC=tapply(DEGs_p$logFC,as.factor(DEGs_p$gene_compare),mean),
  medianp=tapply(DEGs_p$P.Value,as.factor(DEGs_p$gene_compare),median),
  gene_compare=levels(as.factor(DEGs_p$gene_compare))
)

ref_FC_f2$FC <- 2^(mm.DEG$meanlogFC[match(ref_FC_f2$gene_compare,mm.DEG$gene_compare)])
ref_FC_f2$medianp <- mm.DEG$medianp[match(ref_FC_f2$gene_compare,mm.DEG$gene_compare)]

# Output results to CSV file
output_file4 <- paste0(opt$out_dir, "/ref_expr_", Sys.Date(), ".csv")
write.csv(ref_FC_f2, output_file4)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-detected_gene_pairs.csv: ", output_file2))
message(paste("-ref_expr.csv: ", output_file4))