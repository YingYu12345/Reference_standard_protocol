#!/usr/bin/env Rscript
# example:

#Rscript path-to/U.R -d Refdata_uc_(date).csv -r path-to/ref_expr_(date).csv -k 2 -o path-to/

suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 


option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-d", "--uc"),type="character", default=NULL,
              help="File name of combined uncertainty. Required!"),
  make_option(c("-r", "--refdata"),type="character",  default=NULL,
              help="File name of reference data in ratio-scale. Required!"),
  make_option(c("-k", "--k"),type="numeric",  default=NULL,
              help="Coverage factor. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit")
  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

#test
# opt$uc<-"expr_mat/Refdata_uc_example.csv"
# opt$refdata<-"expr_mat/ref_expr_example.csv"
# opt$k<-2

##pre analysis


##import uncertainty files
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

uc<-read.csv(opt$uc,header=T,stringsAsFactors=F,check.names=F)
if(nrow(uc)>0){
  message("Combined uncertainty is loaded.")
}else{
  stop("ERROR: combined uncertainty is not loaded!")
}


refdata<-read.csv(opt$refdata,header=T,stringsAsFactors=F,check.names=F)

if(nrow(refdata)>0){
  message("Input reference dataset is loaded.")
}else{
  stop("ERROR: input reference dataset is not loaded!")
}

k<-opt$k
##analysis

refdata_un<-uc
refdata_un$U<-refdata_un$uc*k


uextendfile<-paste(out_dir,"Refdata_U_",Sys.Date(),".csv",sep="")

####summary
#
print("Summary")
print(tapply(refdata_un$U,refdata_un$compare,summary))


write.csv(refdata_un,uextendfile,row.names=F)

message("Finished!")
message(paste("Output file: \n",uextendfile,sep=" "))