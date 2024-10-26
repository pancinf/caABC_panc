#This script takes caABC predictions and filters them down to
#only replicated ones:
#Open promoter in at least 2 samples
#Enhancer rediscovered in at least one other sample (has to overlap with >=X bp)

###
###
###Libraries + Parameters
set.seed(1234)
options(scipen=999)
library(GenomicRanges)
library(data.table)
library(optparse)

###
###
###External arguments
option_list = list(
  make_option("--caABCRepliList", type="character", default=NULL, 
              help="Path to file containing promoter, enhancer and gene paths to use for replication filtering"),
  
  make_option("--minOvBpEnh", type="numeric", default=NULL, 
              help="Minimum overlap to consider an enhancer replicated (bp)"),
  
  make_option("--outPath", type="character", default=NULL, 
              help="pathToOutputDirectory"),
  
  make_option("--outTag", type="character", default=NULL, 
              help="outTag")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$caABCRepliList) | is.null(opt$minOvBpEnh) | is.null(opt$outPath) | is.null(opt$outTag)) {
  print("Pleas provide all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Main

##
##Functions
repliD <- function(regType,bpOv){
  filePaths <- repliPaths[repliPaths$regType == regType,]
  for(i in 1:nrow(filePaths)){
    subReg <- fread(filePaths$Path[i],data.table = FALSE, header = FALSE)
    subReg <- subReg[subReg$V4 %in% genes$V4,]
    subReg <- makeGRangesFromDataFrame(subReg,keep.extra.columns = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3",starts.in.df.are.0based = TRUE)
    if(i == 1){
      regs <- subReg
    }else{
      regs <- c(regs,subReg)
    }
  }
  regsBp <- sum(width(GenomicRanges::reduce(regs)))
  regs$ovNum <- 0
  
  for(i in 1:nrow(filePaths)){
    subReg <- fread(filePaths$Path[i],data.table = FALSE, header = FALSE)
    subReg <- subReg[subReg$V4 %in% genes$V4,]
    subReg <- makeGRangesFromDataFrame(subReg,keep.extra.columns = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3",starts.in.df.are.0based = TRUE)
    hits <- findOverlaps(regs,subReg,minoverlap = as.integer(bpOv))
    hits <- data.frame(hits@from,hits@to)
    hits$geneFrom <- regs$V4[hits$hits.from]
    hits$geneTo <- subReg$V4[hits$hits.to]
    hits <- hits[hits$geneFrom == hits$geneTo,]
    regs$ovNum[hits$hits.from] <- regs$ovNum[hits$hits.from] +1
  }
  regs <- regs[regs$ovNum > 1,]
  regsBp <- c(regsBp,sum(width(GenomicRanges::reduce(regs))))
  regs <- data.frame(regs)
  regs$start <- regs$start -1
  regs$width <- NULL
  regs$strand <- NULL
  regs$ovNum <- NULL
  regs <- data.table(regs)
  regs <- regs[, as.data.table(reduce(IRanges(start, end))), by = .(seqnames, V4)]
  regs$width <- NULL
  regs <- data.frame(regs[,c(1,3,4,2)])
  return(list(regs,regsBp))
}

##
##Process genes

#Subset to gene paths
repliPaths <- fread(opt$caABCRepliList,data.table = FALSE, header = TRUE)
genePaths <- repliPaths[repliPaths$regType == "Gene",]

#Keep only genes appearing at least twice (--> promoter open)
genes <- fread(genePaths$Path[1],data.table = FALSE, header = FALSE)
for(i in 2:nrow(genePaths)){
  genesSub <- fread(genePaths$Path[i],data.table = FALSE, header = FALSE)
  genes <- rbind(genes,genesSub)
}
genes <- unique(genes[genes$V4 %in% genes$V4[duplicated(genes)],])

##
##Process promoters & enhancers

#Keep only expressed (>= 2x) genes and keep only promoter regions 
#with (> X% overlap with promoters in other samples / matched for target gene)
promFlt <- repliD(regType = "Promoter",bpOv = 1L)
promFltPred <- data.frame(promFlt[1])
#Keep only expressed (>= 2x) genes and keep only enhancer regions 
#with (> X% overlap with enhancer in other samples / matched for target gene)
enhFlt <- repliD(regType = "Enhancer", bpOv = opt$minOvBpEnh)
enhFltPred <- data.frame(enhFlt[1])
enhFltregsBp <- unlist(enhFlt[2])
enhFltregsBp <- data.frame(c("bpBefore","bpAfter","reduction"),c(enhFltregsBp,1 - (enhFltregsBp[2]/enhFltregsBp[1])))

##
##Write outputs
write.table(genes[order(genes$V1,genes$V2),],file = paste0(opt$outPath,opt$outTag,"_genesRepli.bed"),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(promFltPred[order(promFltPred$seqnames,promFltPred$start),],file = paste0(opt$outPath,opt$outTag,"_","minOvBp1","_promotersRepli.bed"),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(enhFltPred[order(enhFltPred$seqnames,enhFltPred$start),],file = paste0(opt$outPath,opt$outTag,"_","minOvBp",opt$minOvBpEnh,"_enhancersRepli.bed"),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(enhFltregsBp,file = paste0(opt$outPath,opt$outTag,"_","minOvBp",opt$minOvBpEnh,"_enhFltregsBpRepli.txt"),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
