#!/usr/bin/env Rscript
# example:
#Rscript path-to/normality.R -i path-to/DEG_limma.csv -m path-to/passbatch.csv -o path-to/


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


gg<-unique(DEGs_p$gene_compare)

normality_mat<-c()

for (i in 1:length(gg)){
  mm<-DEGs_p[DEGs_p$gene_compare==gg[i],]
 if(nrow(mm)>3){
   
 p2<-shapiro.test(mm$logFC)$p.value
 
 normality_mat<-rbind(normality_mat,c(gg[i],p2))
 }
 if(i %in% seq(0,length(gg),1000)){
   
   print(paste(i,"/",length(gg)))
 }
 
}

normality_mat<-data.frame(normality_mat)
colnames(normality_mat)<-c("gene_compare","log2FC_p")
normality_mat$log2FC_p<-as.numeric(as.character(normality_mat$log2FC_p))
normality_mat$gene<-sapply(strsplit(as.character(normality_mat$gene_compare)," "),function(x){x[1]})
normality_mat$compare<-sapply(strsplit(as.character(normality_mat$gene_compare)," "),function(x){x[2]})


s<- round(length(which(normality_mat$log2FC_p>0.05))/nrow(normality_mat)*100,2)
print( paste( s,"% of genes are passed Shapiro-Wilk test (p>0.05). The expresion profiles acorss batches are normally distributed."))


normalityfile<-paste(out_dir,"normality_genelist_",Sys.Date(),".csv",sep="")

write.csv(normality_mat,normalityfile)

message("Finished!")
message(paste("Output file: \n",normalityfile,sep=" "))



