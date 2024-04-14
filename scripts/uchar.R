#!/usr/bin/env Rscript
# example:
#Rscript path-to/uchar.R -i path-to/DEG_limma.csv -m path-to/passbatch.csv -o path-to/

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



#########summary
cutref<-c(0,15,30,50,80,100)

max=max(uchar_mat$uchar)

cut<-c()
for ( i in 1:length(cutref)){
  if(cutref[i]<max){
    
    cut<-c(cut,cutref[i])
  }else(
    cut<-c(cut,max)
  )
  
}
cut<-unique(cut)

cutout<-round(cut,0)
tot<-data.frame(table(uchar_mat$compare))
uchar_summary<-c()
for (i in 2:length(cut)){
  m<-uchar_mat[intersect(which(uchar_mat$uchar>cut[i-1]),which(uchar_mat$uchar<cut[i])),]
  
  U<-data.frame( table(m$compare))
  U$cut0<-cutout[i-1]
  U$cut1<-cutout[i]
  
  
  U$total<-tot$Freq[match(U$Var1,tot$Var1)]
  U$percent<-paste(U$Freq/U$total*100,"%")
  uchar_summary<-rbind(uchar_summary,U)
  
}
uchar_summary$range<-paste(uchar_summary$cut0,"~",uchar_summary$cut1)

uchar_summary_out<-uchar_summary[,c("range","Var1","Freq","percent")]

print(uchar_summary_out)

ucharfile<-paste(out_dir,"Refdata_uchar_",Sys.Date(),".csv",sep="")

write.csv(uchar_mat,ucharfile)

message("Finished!")
message(paste("Output file: \n",normalityfile,sep=" "))




