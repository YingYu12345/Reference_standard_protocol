#!/usr/bin/env Rscript
# example:
#  Rscript path-to/stability_assess.R -i path-to/example_expr_stability_RIN.csv -m path-to/metadata_stability_RIN.txt -o path-to/

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
              help="Show this help message and exit")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

#test
#opt$input<-"expr_mat/example_expr_stability_RIN.csv"
#opt$metadata<-"expr_mat/metadata_stability_RIN.csv"


##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}



##import exp file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

expr<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(expr)>0){
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


t95<-data.frame(
  N=c(1:25),
  t95=c(12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,2.201,2.179,2.16,2.145,2.131,2.12,2.11,2.101,2.093,2.086,2.08,2.074,2.069,2.064,2.06))

if(length(intersect(rownames(expr),metadata$library))==nrow(expr)){
  
  RINvalue_f<-metadata
  RINvalue_f$RIN<-expr$RIN[match(RINvalue_f$library,rownames(expr))]
  
}else{
  stop("ERROR: check metadata and expression!")
}


RINvalue_m<-data.frame(tapply(RINvalue_f$RIN,paste(RINvalue_f$month,RINvalue_f$sample),mean))
colnames(RINvalue_m)<-"RIN"
RINvalue_m$time<-sapply(strsplit(rownames(RINvalue_m)," "),function(x){x[1]})
RINvalue_m$sample<-sapply(strsplit(rownames(RINvalue_m)," "),function(x){x[2]})
RINvalue_m$time_sample<-rownames(RINvalue_m)


########
### relative value

usample<-as.character(unique(RINvalue_m$sample))
udenomsample<-usample[1]
unumsample<-usample[2:length(usample)]


RINvalue_ratio<-c()
utime<-unique(RINvalue_m$time)

for ( i in 1:length(utime)){
  m<-RINvalue_m[RINvalue_m$time==utime[i],]
  m$RINvalue_ratio<-m$RIN/m$RIN[m$sample==udenomsample]
  
  RINvalue_ratio<-rbind(RINvalue_ratio,m)
}


stability<-c()
for ( i in 1:length(unumsample)){
 m1<-RINvalue_ratio[RINvalue_ratio$sample==unumsample[i],c("time","RINvalue_ratio")]
 m1$time<-as.numeric(as.character(m1$time))
 m1<-m1[order(m1$time),]
 
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
stability<-rbind(stability,c(paste(unumsample[i],"/",udenomsample,sep=""),mean(m1$Y),b1,b0,n,s,s_b1,t95_n,sig))

}


colnames(stability)<-c("compare","meanFC","b1_slope","b0_intercept","n","s","sb1","t95_n","sig")

stability<-data.frame(stability)

for (i in c(1,9)){
  stability[,i]<-as.character( stability[,i])
}

for (i in 2:8){
  stability[,i]<-as.numeric(as.character( stability[,i]))
}


##output


stabilityfile<-paste(out_dir,"stability_lm_detail_",Sys.Date(),".csv",sep="")

write.csv(stability,stabilityfile)

message("Finished!")
message(paste("Output file: \n",stabilityfile,sep=" "))


