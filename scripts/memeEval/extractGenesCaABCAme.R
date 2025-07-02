#This script extracts regulatory sequences for defined genes from caABC results.
###
###
###Libraries
set.seed(1234)
options(scipen = 999)
library(data.table)
library(optparse)

###
###
###External arguments
option_list = list(
  make_option("--caABCName", type="character", default=NULL, 
              help="name of the caABC file (before .txt.gz ending and without path)"),
  
  make_option("--geneList", type="character", default=NULL, 
              help="One-column list of genes (Ensembl-IDs) to investigate"),
  
  make_option("--out", type="character", default=NULL, 
              help="path to output file containing short IDs for regulatory regions")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###
###
###Main
caABC <- fread(paste0("../../data/GTEx/caAbcRepli/caAbcRepliOut/",opt$caABCName,".txt"), data.table = FALSE, header = TRUE)

#Resize regions
midPoints <- round(rowMeans(caABC[,c(3,4)]))
caABC$Start <- midPoints - 250
caABC$End <- midPoints + 250
geneList <- fread(paste0("../../data/memeEval/geneLists/",opt$geneList,".txt"),data.table = FALSE,header = FALSE)

#Cases
caABCCases <- caABC[caABC$Symbol %in% geneList$V1 & caABC$featureType %in% c("Promoter", "Enhancer"),]
caABCCases <- paste0(caABCCases$Chr,":",caABCCases$Start,"-",caABCCases$End)
caABCCases <- unique(caABCCases)
write.table(x = caABCCases, file = paste0("../../data/memeEval/ame/",opt$caABCName,"_",opt$geneList,"/","regRegionsCase.txt"),quote = FALSE,append = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

#Controls
caABCControls <- caABC[!(caABC$Symbol %in% geneList$V1) & caABC$featureType %in% c("Promoter", "Enhancer"),]
controlGenes <- sample(unique(caABCControls$Symbol), 5000)
caABCControls <- caABCControls[caABCControls$Symbol %in% controlGenes,]
caABCControls <- paste0(caABCControls$Chr,":",caABCControls$Start,"-",caABCControls$End)
caABCControls <- unique(caABCControls)
caABCControls <- caABCControls[!(caABCControls %in% caABCCases)]
write.table(x = caABCControls, file = paste0("../../data/memeEval/ame/",opt$caABCName,"_",opt$geneList,"/","regRegionsControl.txt"),quote = FALSE,append = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
