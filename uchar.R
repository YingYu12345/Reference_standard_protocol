#!/usr/bin/env Rscript
# example:

#Rscript path-to/uchar.R -i path-to/DEG_limma.csv -m path-to/passbatch.csv -r path-to/ref_expr_(date).csv -o path-to/

suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 


option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"),type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-p", "--passbatch"),type="character",  default=NULL,
              help="List of passbatch. Required!"),
  make_option(c("-r", "--refdata"),type="character",  default=NULL,
              help="File name of reference data in ratio-scale. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

#test
#opt$input<-"data/DEG_limma_2024-04-02.csv"
#opt$passbatch<-"expr_mat/passbatch.csv"
#opt$refdata<-"expr_mat/ref_expr_example.csv"

##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}



##import exp file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")



DEGs<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(DEGs)>0){
  message("Input expression profile is loaded.")
}else{
  stop("ERROR: input expression profile is not loaded!")
}



passbatch<-read.csv(opt$passbatch,header=T,stringsAsFactors=F,check.names=F)

passbatch<-as.character(passbatch[,1])

if(length(passbatch)>0){
  message("List of passed batches is loaded.")
}else{
  stop("ERROR: list of passed batches is not loaded!")
}


refdata<-read.csv(opt$refdata,header=T,stringsAsFactors=F,check.names=F)


if(nrow(refdata)>0){
  message("Input reference dataset is loaded.")
}else{
  stop("ERROR: input reference dataset is not loaded!")
}

##analysis

DEGs_p<-DEGs[DEGs$batch %in% passbatch,]
DEGs_p$gene_compare<-paste(DEGs_p$gene,DEGs_p$compare)



rse <- function(x){sd(x, na.rm = T)/(mean(x, na.rm=T)*sqrt(length(x[!is.na(x)])))}
gene_compare=levels(as.factor(DEGs_p$gene_compare))


uchar_mat<-data.frame(
  uchar=tapply(DEGs_p$logFC,as.factor(as.character(DEGs_p$gene_compare)),function(x){rse(2^x)*100}),
  gene_compare=gene_compare,
  gene=sapply(strsplit(gene_compare," "),function(x){x[1]}),
  compare=sapply(strsplit(gene_compare," "),function(x){x[2]}))

refdata_uchar<-refdata

refdata_uchar$uchar<-uchar_mat$uchar[match(refdata_uchar$gene_compare,uchar_mat$gene_compare)]


###out
refdata_uchar_out<-refdata_uchar[,c("gene_compare","gene","compare","uchar")]


ucharfile<-paste(out_dir,"Refdata_uchar_",Sys.Date(),".csv",sep="")

####summary
#
print("Summary")
print(tapply(refdata_uchar$uchar,refdata_uchar$compare,summary))


write.csv(refdata_uchar_out,ucharfile,row.names=F)

message("Finished!")
message(paste("Output file: \n",ucharfile,sep=" "))




