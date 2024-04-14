#!/usr/bin/env Rscript
# example:
# Rscript path-to/homogeneity_assess.R -i path-to/example_expr_homo_log2.csv -m metadata_homo.csv -o path-to/

suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 


option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"),type="character", default=NULL,
              help="File name of expression file in log2 scale. Required!"),
  make_option(c("-m", "--metadata"),type="character",  default=NULL,
              help="Metadata of the input expression file. Note: metadata file shall include two columns: sample and type (intra- or between-unit). Required!"),
  make_option(c("-h", "--help"), action="store_true", default=FALSE, 
            help="Show this help message and exit")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))
#

#test
#opt$input<-"expr_mat/example_expr_homo_log2.csv"
#opt$metadata<-"expr_mat/metadata_homo.csv"


##pre analysis
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}



##import exp file
out_dir<-paste(gsub("/$","",opt$out_dir),"/",sep="")

logexpr<-read.csv(opt$input,header=T,stringsAsFactors=F,row.names=1,check.names=F)

if(nrow(logexpr)>0){
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

##analysis

homoge_p<-c()
homoge<-c()
for(j in 1:length(usample)){
  
  for (i in 1:nrow(logexpr)){
    g<-rownames(logexpr)[i]
    
    df<-metadata[metadata$sample==usample[j],]
    
    df$value<-as.numeric(logexpr[rownames(logexpr) ==g,match(rownames(df), colnames(logexpr))])
    
    p<- summary(aov(df$value~df$type))[[1]][["Pr(>F)"]][1]  
    means<-tapply(df$value,df$type,mean)
    
    homoge_p<-rbind(homoge_p,c(usample[j],g,p,means)) 
    
  }
  
  
  colnames(homoge_p)<-c("sample","gene","anova_p","mean_between","mean_intra")
  homoge_p<-data.frame(homoge_p)
  homoge_p$gene<-as.character(homoge_p$gene)
  
  for ( k in 3:ncol(homoge_p)){
    homoge_p[,k]<-round(as.numeric(as.character( homoge_p[,k])),4)
  }
  
  homoge_p$fdr<-round(p.adjust(homoge_p$anova_p,method="fdr"),4)
  
  homoge<-rbind(homoge,homoge_p)
  
   homoge_p<-c()
}

sums<-data.frame(
  sample=levels(as.factor(homoge$sample)),
  length_sig=tapply(homoge$fdr,homoge$sample,function(x){length(which(x<0.05))}),
  length_total=tapply(homoge$fdr,homoge$sample,function(x){length(x)}),
  length_ratio=tapply(homoge$fdr,homoge$sample,function(x){length(which(x<0.05))/length(x)})
)
rownames(sums)<-seq(1:nrow(sums))

##output

homogefile<-paste(out_dir,"homogeneity_Ftest_detail_",Sys.Date(),".csv",sep="")
sumfile<-paste(out_dir,"homogeneity_Ftest_summary_",Sys.Date(),".csv",sep="")

write.csv(homoge,homogefile)
write.csv(sums,sumfile)

message("Finished!")
message(paste("Output file: \n",homogefile,"\n and \n",sumfile,sep=" "))
