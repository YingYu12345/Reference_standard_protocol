#!/usr/bin/env Rscript
# example:

#Rscript path-to/uc.R -c Refdata_uchar_(date).csv -b Refdata_ubb_(date).csv -s Refdata_us_(date).csv -r path-to/ref_expr_(date).csv -o path-to/

suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 


option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-c", "--uchar"),type="character", default=NULL,
              help="File name of characterization uncertainty. Required!"),
  make_option(c("-b", "--ubb"),type="character",  default=NULL,
              help="File name of inhomogeneity uncertainty. Required!"),
  make_option(c("-s", "--us"),type="character",  default=NULL,
              help="File name of instability  uncertainty. Required!"),
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
opt$uchar<-"expr_mat/Refdata_uchar_example.csv"
opt$ubb<-"expr_mat/Refdata_ubb_example.csv"
opt$us<-"expr_mat/Refdata_us_example.csv"
opt$refdata<-"expr_mat/ref_expr_example.csv"

##pre analysis


##import uncertainty files
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

uchar<-read.csv(opt$uchar,header=T,stringsAsFactors=F,check.names=F)
if(nrow(uchar)>0){
  message("Characterization uncertainty is loaded.")
}else{
  stop("ERROR: characterization uncertainty is not loaded!")
}


ubb<-read.csv(opt$ubb,header=T,stringsAsFactors=F,check.names=F)
if(nrow(ubb)>0){
  message("Inhomogeneity uncertainty is loaded.")
}else{
  stop("ERROR: inhomogeneity uncertainty is not loaded!")
}


us<-read.csv(opt$us,header=T,stringsAsFactors=F,check.names=F)
if(nrow(us)>0){
  message("Instability uncertainty is loaded.")
}else{
  stop("ERROR: Instability uncertainty is not loaded!")
}


refdata<-read.csv(opt$refdata,header=T,stringsAsFactors=F,check.names=F)

if(nrow(refdata)>0){
  message("Input reference dataset is loaded.")
}else{
  stop("ERROR: input reference dataset is not loaded!")
}

##analysis

refdata_un<-refdata

refdata_un$uchar<-uchar$uchar[match(refdata_un$gene_compare,uchar$gene_compare)]
refdata_un$ubb<-ubb$ubb[match(refdata_un$gene_compare,ubb$gene_compare)]
refdata_un$us<-us$us[match(refdata_un$gene_compare,us$gene_compare)]

refdata_un$uc<-sqrt(refdata_un$uchar^2+refdata_un$ubb^2+refdata_un$us^2)

ucombinefile<-paste(out_dir,"Refdata_uc_",Sys.Date(),".csv",sep="")

####summary
#
print("Summary")
print(tapply(refdata_un$uc,refdata_un$compare,summary))


write.csv(refdata_un,ucombinefile,row.names=F)

message("Finished!")
message(paste("Output file: \n",ucombinefile,sep=" "))

