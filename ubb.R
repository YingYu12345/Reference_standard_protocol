#!/usr/bin/env Rscript
# example:


#Rscript path-to/ubb.R -i path-to/example_expr_homo_log2.csv -m path-to/metadata_homo.csv -r path-to/ref_expr_(date).csv -o path-to/
  
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 

option_list <- list( 
  make_option(c("-o", "--out_dir"), type="character",default="./",
              help="The output directory [default ./]"),
  make_option(c("-i", "--input"),type="character", default=NULL,
              help="File name of expression file in ratio-scale. Required!"),
  make_option(c("-m", "--metadata"),type="character",  default=NULL,
              help="Metadata of the input expression file. Note: metadata file shall include two columns: sample and type (intra- or between-unit). Required!"),
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

#opt$input<-"expr_mat/example_expr_homo_ratio.csv"
#opt$metadata<-"expr_mat/metadata_homo.csv"
#opt$refdata<-"expr_mat/ref_expr_example.csv"



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

##analysis


sbb<-function(intra_list,between_list){
  s2<-var(intra_list)
  s1<-var(between_list)
  
  sbb<-sqrt((s1-s2)/length(intra_list))
  return(sbb)
}





df<-metadata
df$library<-rownames(df)
df1<-df[match(colnames(ratio_hom1),df$library),]

usample3<-as.character(unique(df$sample))

ubb_all<-c()
for ( i in 1:nrow(refdata)){
  g<-refdata$gene[i]
  s<-refdata$usample[i]
  
  df2<-df1[df1$sample==s,]
  
  df2$value<-ratio_hom1[rownames(ratio_hom1)==g,match(df2$library,colnames(ratio_hom1))]
  
  ############log2FC ubb_1
  intra_list<-df2$value[df2$type=="intra-unit"]
  between_list<-df2$value[df2$type=="between-unit"]
  
  s2_1<-var(intra_list)
  s1_1<-var(between_list)
  
  ubb_1_logFC<-sbb(intra_list,between_list)*100
  
  #############FC ubb_1
  intra_list<-2^(df2$value[df2$type=="intra-unit"])
  between_list<-2^(df2$value[df2$type=="between-unit"])
  
  s2_2<-var(intra_list)
  s1_2<-var(between_list)
  
  ubb_1_FC<-sbb(intra_list,between_list)*100
  
  ############logFC ubb_2
  intra_list<-df2$value[df2$type=="intra-unit"]
  between_list<-df2$value[df2$type=="between-unit"]
  
  s2_3<-var(intra_list) #### s2^2
  s1_3<-var(between_list) ### s1^2
  
  n<-length(intra_list)
  v<-n-1
  
  ubb_2_logFC<-sqrt(s2_3/n)*((2/v)^0.25)*100
  
  
  ############FC ubb_2
  intra_list<-2^(df2$value[df2$type=="intra-unit"])
  between_list<-2^(df2$value[df2$type=="between-unit"])
  
  s2_4<-var(intra_list) ####jiushi s2^2
  s1_4<-var(between_list) ###jiushi s1^2
  
  n<-length(intra_list)
  v<-n-1
  
  ubb_2_FC<-sqrt(s2_4/n)*((2/v)^0.25)*100
  
  ubb_all<-rbind(ubb_all,c(g,refdata$compare[i],refdata$gene_compare[i],s2_1,s1_1,ubb_1_logFC,s2_2,s1_2,ubb_1_FC,s2_3,s1_3,ubb_2_logFC,s2_4,s1_4,ubb_2_FC))
  
  if(i %in% seq(0,nrow(refdata),by=100)){
    print(i)
  }
}


colnames(ubb_all)<-c("gene","compare","gene_compare","S2_1","S1_1","ubb_1_logFC","S2_2","S1_2","ubb_1_FC","S2_3","S1_3","ubb_2_logFC","S2_4","S1_4","ubb_2_FC")
ubb_all<-data.frame(ubb_all)

for(i in 4:ncol(ubb_all)){
  ubb_all[,i]<-as.numeric(as.character(ubb_all[,i]))
}


ubb_all_f<-ubb_all[,c("gene","compare","gene_compare","S2_2","S1_2","ubb_1_FC","ubb_2_FC")]

ubb_all_f$ubb_2_FC[!is.na(ubb_all_f$ubb_1_FC)]<-NA

ubb_all_f$ubb_all<-ubb_all_f$ubb_1_FC
ubb_all_f$ubb_all[is.na(ubb_all_f$ubb_1_FC)]<-ubb_all_f$ubb_2_FC[is.na(ubb_all_f$ubb_1_FC)]


refdata_Ubb<-merge(refdata,ubb_all_f,by.x="gene_compare",by.y="gene_compare")

refdata_Ubb_out<-refdata_Ubb[,c("gene.x","compare.x","gene_compare","ubb_all")]
colnames(refdata_Ubb_out)<-c("gene","compare","gene_compare","ubb")

ubbfile<-paste(out_dir,"Refdata_ubb_",Sys.Date(),".csv",sep="")

write.csv(refdata_Ubb_out,ubbfile,row.names=F)

####summary
#
print("Summary")
print(tapply(refdata_Ubb_out$ubb,refdata_Ubb_out$compare,summary))


message("Finished!")
message(paste("Output file: \n",ubbfile,sep=" "))




