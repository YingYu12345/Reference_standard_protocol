#!/usr/bin/env Rscript
# example:
#  Rscript path-to/us.R -i path-to/example_expr_stability_RIN.csv -m path-to/metadata_stability_RIN.csv -r path-to/ref_expr.csv -t R -o path-to/

suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library(""))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 


option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"),type="character", default=NULL,
              help="File name of gene expression file in log2 scale, or file of RIN values. Required!"),
  make_option(c("-m", "--metadata"),type="character",  default=NULL,
              help="Metadata of the input expression file. Note: metadata file shall include at least three columns: sample (library), replicates, and sample group. Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
              help="Show this help message and exit"),
  make_option(c("-r", "--refdata"),type="character",  default=NULL,
              help="File name of reference data in ratio-scale. Required!"),
  make_option(c("-t", "--type"),type="character",  default=NULL,
              help="R or G. R for RIN values, and G for calculation based of genes. Required!")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

#test
#opt$input<-"expr_mat/example_RINvalue_ratio.csv"
#opt$metadata<-"expr_mat/metadata_stability_RIN.csv"
#opt$refdata<-"expr_mat/ref_expr_example.csv"
#opt$type<-"R"

##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}



##import exp file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

ratioexpr<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(ratioexpr)>0){
  message("Input expression profile is loaded.")
}else{
  stop("ERROR: input expression profile is not loaded!")
}

metadata<-read.csv(opt$metadata,header=T,stringsAsFactors=F,row.names=1,check.names=F)

usample<-unique(metadata$sample)


if(nrow(metadata)>0){
  message("Metadata is loaded.")
}else{
  stop("ERROR: metadata is not loaded!")
}


refdata<-read.csv(opt$refdata,header=T,stringsAsFactors=F,check.names=F)
refdata$usample<-sapply(strsplit(as.character(refdata$compare),"/"),function(x){x[1]})

if(nrow(refdata)>0){
  message("Input reference dataset is loaded.")
}else{
  stop("ERROR: input reference dataset is not loaded!")
}


t95<-data.frame(
  N=c(1:25),
  t95=c(12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,2.201,2.179,2.16,2.145,2.131,2.12,2.11,2.101,2.093,2.086,2.08,2.074,2.069,2.064,2.06))


##analysis

if (opt$type=="R"){
######use RIN values 
  
  ratioexpr$library<-rownames(ratioexpr)
  
stability<-c()
for ( i in 1:length(usample)){
  lib<-metadata$library[metadata$sample==usample[i]]
  
  m<-ratioexpr[rownames(ratioexpr)%in% lib,]
  m$month<-metadata$month[match(rownames(m),metadata$library)]
  
  mm<-data.frame(tapply(m$value,m$month,mean))
  colnames(mm)<-"value"
  mm$month<-as.numeric(levels(as.factor(m$month)))
  
  mm<-mm[order(mm$month),]
  
  m1<-mm[,c("month","value")]
  colnames(m1)<-c("X","Y")
  
  b1<-sum((m1$X-mean(m1$X))*(m1$Y-mean(m1$Y)))/sum(((m1$X-mean(m1$X))^2))
  
  b0<-mean(m1$Y)-b1*mean(m1$X)
  
  n<-nrow(m1)
  
  t95_n<-t95$t95[t95$N==(n-2)]
  
  s<-sqrt(sum((m1$Y-b0-b1*m1$X)^2)/(n-2))
  
  s_b1<-s/sqrt(sum((m1$X-mean(m1$X))^2))
  
  if(abs(b1)<(t95_n*s_b1)){
    sig<-F
  }else{
    sig<-T
  }
  stability<-rbind(stability,c(usample[i],mean(m1$Y),b1,b0,n,s,s_b1,t95_n,sig))
  
}

colnames(stability)<-c("usample","meanFC","b1_slope","b0_intercept","n","s","sb1","t95_n","sig")

stability<-data.frame(stability)

for (i in c(1,9)){
  stability[,i]<-as.character( stability[,i])
}

for (i in 2:8){
  stability[,i]<-as.numeric(as.character( stability[,i]))
}

stability$sb1_t095<-stability$sb1*stability$t95_n
stability$ub<-stability$sb1_t095*20*100

#write.csv(stability,"refdata_biaowu_20220426/expr_mat_refdata/wending_RIN_20220724.csv")

refdata$us<-stability$ub[match(refdata$usample,stability$usample)]
  
  
}


###out

refdata_us<-refdata

refdata_us_out<-refdata_us[,c("gene_compare","gene","compare","us")]


usfile<-paste(out_dir,"Refdata_us_",Sys.Date(),".csv",sep="")


write.csv(refdata_us_out,usfile,row.names=F)

message("Finished!")
message(paste("Output file: \n",usfile,sep=" "))






