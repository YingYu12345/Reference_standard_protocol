#!/usr/bin/env Rscript
# SNR_qc.R
# example:
# Rscript path-to/SNR_qc.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -s sample_info(eg. "D5,D6,F7,M8") -c config_file -o path-to/

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
required_packages <- c("optparse", "data.table", "ggplot2")

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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# Function to print usage information
print_usage <- function() {
  cat("Actual Usage: Rscript SNR_qc.R [options]\n\n")
  cat("Example:\n")
  cat("  Rscript path-to/SNR_qc.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -s \"D5,D6,F7,M8\"\n -c config_file -o path-to/")
  cat("Note: Make sure to replace path-to/ with your actual file paths\n\n")
}
print_usage()

# Define command-line option specifications
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input file, including the expression matrix and metadata."),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to the metadata file of the input expression. Note: The metadata file must include three columns: batch, sample and type. Required!"),
  make_option(c("-o", "--out_dir"), type="character", default="./",
              help="Path to the output directory [default current directory]."),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Comma-separated list of sample names (e.g., \"D5,D6,F7,M8\"). Required!"),
  make_option(c("-c", "--config"), type="character", default="./scripts/config.R",
            help="Path to the configuration file [default ./scripts/config.R]"),
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
if (is.null(opt$input) || is.null(opt$metadata) || is.null(opt$out_dir) || is.null(opt$samples) || is.null(opt$config)) {
  print_help(opt_parser)
  stop("Input file, metadata file, output directory, samples parameters, and config file must be provided.", call.=FALSE)
}

# Parse the samples and configparameter
samples <- unlist(strsplit(opt$samples, ","))
source(opt$config)

# define SNR calculation function
calSNR_pca <- function(exprMat, group) {
  library(data.table)
  
  IDs <- colnames(exprMat)
  IDs.group.mat <- data.table(IDs = IDs, group = group)
  
  # Remove features which variance is zero to ensure the PCA
  exprMat <- exprMat[which(apply(exprMat, 1, var) != 0), ]
  
  # PCA
  # pca_prcomp <- prcomp(t(exprMat), retx = TRUE, center = FALSE, scale. = FALSE)
  pca_prcomp <- prcomp(t(exprMat),scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  # Percent: Proportion of Variance, AccumPercent: Cumulative Proportion
  dt.perc.pcs <- data.table(PCX = 1:length(summary(pca_prcomp)$importance[2,]),
                            Percent = summary(pca_prcomp)$importance[2,],
                            AccumPercent = summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs, each = length(IDs)),
                        ID.B = rep(IDs, time = length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A, IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B, IDs.group.mat$IDs)]$group
  
  dt.dist[, Type := ifelse(ID.A == ID.B, 'Same',
                           ifelse(group.A == group.B, 'Intra', 'Inter'))]
  
  dt.dist[, Dist := (dt.perc.pcs[1]$Percent * (pcs[ID.A, 1] - pcs[ID.B, 1])^2 +
                       dt.perc.pcs[2]$Percent * (pcs[ID.A, 2] - pcs[ID.B, 2])^2)]
  
  dt.dist.stats <- dt.dist[, .(Avg.Dist = mean(Dist)), by = .(Type)]
  setkey(dt.dist.stats, Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist / dt.dist.stats['Intra']$Avg.Dist
  
  signoise_db <- 10 * log10(signoise)
  
  return(list(signoise_db = signoise_db, pca_prcomp = pca_prcomp))
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
nbatch <- unique(meta$batch)
SNR_n <- data.frame("batch_id"=character(), "SNR"=numeric(), stringsAsFactors=FALSE)

type <- function(x, samples) {
  type <- rep(NA, length(x)) # Initialize the type vector with NA values
  for (sample in samples) {
    type[grepl(sample, x)] <- sample
  }
  return(type)
}

# Calculate the SNR for each batch
for (i in seq_along(nbatch)) {
  logFPKM_nb <- logFPKM[, grep(nbatch[i], colnames(logFPKM)), drop=FALSE]
  ntype <- type(colnames(logFPKM_nb), samples)
  # print(ntype) # Debugging line
  SNR_pca <- calSNR_pca(logFPKM_nb, ntype)
  SNR_n <- rbind(SNR_n, data.frame(batch_id=nbatch[i], SNR=round(SNR_pca$signoise_db, 1)))
}

SNR_n$SNR <- as.numeric(SNR_n$SNR)
SNR_n <- SNR_n[order(SNR_n$SNR, decreasing = TRUE), ]
colnames(SNR_n) <- c("batch_id","snr")

# SNR Output results to CSV file
output_csv <- paste0(opt$out_dir, "/SNR_perBatch_", Sys.Date(), ".csv")
fwrite(SNR_n, output_csv)

## plot PCA with SNR
pca_n <- list()

plot_pca <- function(SNR,image_title,expr,image_name, samples) {
  SNR_n<- round(SNR$signoise_db,1)
  pca.all <- SNR$pca_prcomp
  pcs <- pca.all$x
  pcs <- as.data.frame(pcs)
  pcs <- pcs[,1:3]
  pcs$library <- rownames(pcs)
  
  # Merge with metadata
  pcs <- merge(pcs, meta, by="library")
  pcs$sample <- factor(pcs$sample, levels=samples)
  
  # Create plot
  p <-  ggplot(pcs, aes(x = PC1, y = PC2, fill = sample, color = sample, shape = batch)) +
    geom_point(size = 3.5) +
    theme_bw() +
    theme(
      axis.text = element_text(),
      legend.position = "none",
      legend.background = element_blank(),
      plot.title = element_text(hjust=0.5, color="black", face="bold"),
      plot.subtitle = element_text(hjust=0.5)
    ) +
    # guides(
    #   fill = guide_legend("Type", override.aes = list(shape=colors, fill=shapes))
    # ) +
    scale_fill_manual("Type", values = colors) +  # Ensure `colors` is a named vector
    scale_color_manual("Type", values = colors) +  # Ensure `colors` is a named vector
    scale_shape_manual("Batch", values = shapes) +
    labs(
      x=sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
      y=sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
      title=paste0(image_title," ", nrow(expr)),
      subtitle=paste0("SNR=", SNR_n)
    ) +
    theme(aspect.ratio=1/1)
  
  # Save plot to file
  # output_pdf1 <- paste0(opt$out_dir, "/", image_name, "_", Sys.Date(), ".pdf")
  # ggsave(output_pdf1, p, width=5, height=5, dpi=300)
  return(p)
}

## sort_batch by SNR RANK 
sbatch <- SNR_n$batch_id
SNR_s <- data.frame("batch_id"=1,"snr"=1)
SNR_s <- SNR_s[-1,]

for (i in 1:length(sbatch))
{
  print(i);print(sbatch[i])
  logFPKM_nb <- logFPKM[,grep(sbatch[i],colnames(logFPKM))]
  
  ntype <- type(colnames(logFPKM_nb), samples)
  SNR_pca <- calSNR_pca(logFPKM_nb,ntype)
  temp2<- round(SNR_pca$signoise_db,1)
  SNR_s <- rbind(SNR_s,c(sbatch[i],temp2))
  p <- plot_pca(SNR_pca,paste0("SNR","_",sbatch[i]),logFPKM_nb,paste0("SNR","_",sbatch[i]))
  pca_n[[i]]=p
}

# # export the all pdf files to one pdf file
output_pdf2 <- paste0(opt$out_dir, "/SNR_perBatch_", Sys.Date(), ".pdf")
pdf(file = output_pdf2, width = 5, height = 5)
for (i in 1:length(sbatch))
{
  print(pca_n[[i]])
}
dev.off()

## all batches
ntype <- type(colnames(logFPKM), samples)
SNR_pca <- calSNR_pca(logFPKM,ntype)
p <- plot_pca(SNR_pca,"SNR_allBatch",logFPKM,"SNR_allBatch_")
p <- p +
    theme(
      legend.position = "right",
      legend.background = element_blank()
    )
output_pdf3 <-  paste0(opt$out_dir, "/", "SNR_allBatch", "_", Sys.Date(), ".pdf")
ggsave(output_pdf3, p, width=5, height=5, dpi=300)

# Print completion message
message("Analysis complete! Results saved to: ")
message(paste("-SNR_perBatch.csv: ", output_csv))
# message(paste("-PCA_withSNR_perBatch.pdf: ", paste0(opt$out_dir, "/")))
message(paste("-PCA_withSNR_perBatch.pdf: ", output_pdf2))
message(paste("-PCA_withSNR_allBatch.png: ", output_pdf3))