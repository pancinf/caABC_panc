#This script evaluates ABC models on eQTL data by the Sign test and also shows top differences for relevant enhancers

###
###
###Libraries and parameters
set.seed(1234)
options(scipen=999)
library(data.table)
library(reshape2)
library(plyranges)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(forcats)
library(optparse)
library(rstatix)
library(plyr) ######################### LOAD BEFORE DPLYR #########################
library(dplyr)

###
###
###Arguments
option_list = list(
  
  make_option("--inputList", type="character", default=NULL,
              help="Name of input list to use"),
  
  make_option("--minDist", type="numeric", default=NULL,
              help="Define minimum distance to TSS to define an enhancer. Set to -1 to consider all promoters"),
 
  make_option("--gtexCut", type="numeric", default=NULL,
              help="Minimum TPM value of GTEx expression for a gene to be included"),
  
  make_option("--caviarThres", type="numeric", default=NULL,
              help="eQTL cut-off"),

  make_option("--tissue", type="character", default=NULL,
              help="Tissue to filter TPM and caviar file on"),
			  
  make_option("--cap", type="numeric", default=NULL,
              help="Given cap for outliers in boxplot. Equally applied to upper/lower"),
  
  make_option("--compMod ", type="character", default=NULL,
              help="Tag. Which models are compared?"),

  make_option("--outDir", type="character", default=NULL,
              help="Name of output directory")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$inputList) | is.null(opt$minDist) |
    is.null(opt$gtexCut) | is.null(opt$caviarThres) | is.null(opt$tissue) |
    is.null(opt$cap) | is.null(opt$compMod) | is.null(opt$outDir)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Plot function
plotBox <- function(diffTable, plotMode,capVal){
	p <- ggplot(diffTable, aes(x = "",y =Diff)) +geom_boxplot(outlier.shape = NA, width = 0.4) + geom_jitter(cex = 0.5,width = 0.1, alpha = 0.5)+
		xlab(NULL) + ylab("Pairwise difference of ranks")+
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			legend.key=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=14))+
			scale_x_discrete()
	ggsave(filename = paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_plotMode", plotMode,capVal,"_BoxDiff.png"), plot = p,width = 3.5,height = 4.38,dpi = 600)
	ggsave(filename = paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_plotMode", plotMode,capVal,"_BoxDiff.pdf"), plot = p,width = 3.5,height = 4.38)
}

###
###
###Main

##
##Make output directory and read inputList
dir.create(paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir), showWarnings = FALSE)
inputList <- fread(paste0("../../data/GTEx/evalGTEx/inputLists/",opt$inputList), data.table = FALSE, header = TRUE)

##
##Read in gtex genes and filter on expression
gtex <- fread(input = "../../data/GTEx/raw/GTExTPM.gct",data.table = FALSE, header = TRUE)
gtex <- gtex[,colnames(gtex) %in% c("Name","Description",opt$tissue)]
colnames(gtex) <- c("ensemblID","symbol","tissueTPM")
gtex$ensemblID <- gsub("\\..*","",gtex$ensemblID)
gtex <- gtex[gtex$tissueTPM >= opt$gtexCut,]

##
##Process caviar data
caviar <- fread(input = "../../data/GTEx/raw/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz", data.table = FALSE, header = TRUE)
colnames(caviar)[2] <- "ensemblID_CAVIAR"
caviar$ensemblID_CAVIAR <- gsub("\\..*","",caviar$ensemblID_CAVIAR)
caviar <- caviar[caviar$Probability >= opt$caviarThres,]
caviar <- caviar[caviar$TISSUE == opt$tissue,]
caviar <- caviar[caviar$ensemblID_CAVIAR %in% gtex$ensemblID,]

##
##Process GTEx gtf data
gtexGtf <- fread(input = "../../data/GTEx/raw/gencodeGTEx.gtf", data.table = FALSE, header = FALSE)
gtexGtf <- gtexGtf[gtexGtf$V3 == "gene",]
gtexGtfColsplit <- colsplit(gtexGtf$V9, pattern = "; ",names = c("ensemblID","transcript_id","gene_type","symbol","remaining"))
gtexGtfColsplit <- gtexGtfColsplit[,c(1,4)]

gtexGtfColsplit$ensemblID <- gsub(pattern = "gene_id \"",replacement = "",x = gtexGtfColsplit$ensemblID)
gtexGtfColsplit$ensemblID <- gsub(pattern = "\"",replacement = "",x = gtexGtfColsplit$ensemblID)
gtexGtfColsplit$ensemblID <- gsub(pattern = "\\..*",replacement = "",x = gtexGtfColsplit$ensemblID)

gtexGtfColsplit$symbol <- gsub(pattern = "gene_name \"",replacement = "",x = gtexGtfColsplit$symbol)
gtexGtfColsplit$symbol <- gsub(pattern = "\"",replacement = "",x = gtexGtfColsplit$symbol)
gtexGtfColsplit$symbol <- gsub(pattern = "\\..*",replacement = "",x = gtexGtfColsplit$symbol)

gtexGtfColsplit <- paste0(gtexGtfColsplit$ensemblID,"_",gtexGtfColsplit$symbol)
gtexGtf <- data.frame(chrGtf = gtexGtf$V1, startGtf = gtexGtf$V4, endGtf = gtexGtf$V5, strand = gtexGtf$V7, geneSymbol = gtexGtfColsplit)
gtexGtfPlus <- gtexGtf[gtexGtf$strand == "+",]
gtexGtfPlus$endGtf <- NULL
gtexGtfPlus$strand <- NULL
colnames(gtexGtfPlus)[2] <- "TSS_PosGtf"
gtexGtfMinus <- gtexGtf[gtexGtf$strand == "-",]
gtexGtfMinus$startGtf <- NULL
gtexGtfMinus$strand <- NULL
colnames(gtexGtfMinus)[2] <- "TSS_PosGtf"
gtexGtf <- rbind(gtexGtfPlus,gtexGtfMinus)

##
##Merge ABC and CAVIAR data
pltRec <- c()
pltAbc <- c()
pltRank <- c()
pltPerc <- c()
pltName <- c()
plteQTL <- c()
plteQTLMatch <- c()
pltEnsemblID <- c()
pltpeakEnsemblID <- c()

##
##Make table of overlapping regions for both models
for(i in 1:nrow(inputList)){
  abc <- fread(inputList$samplePath[i],data.table = FALSE, header = TRUE)
  colnames(abc) <- c("chrom","start","end","ensemblID","symbol","peakID","signalValue","contact","adaptedActivity","scaledContact","intergenicScore","TSS_dist","ABC_Score")

  abc <- abc[abc$TSS_dist > opt$minDist,]
  abc$ensemblID <- gsub("\\..*","",abc$ensemblID)
  abc <- abc[abc$ensemblID %in% gtex$ensemblID,]
  abc$modelName <-inputList$modelName[i]
  abc$peakEnsemblID <- paste0(abc$ensemblID,"_",abc$peakID)
  abc$geneSymbol <- paste0(abc$ensemblID,"_",abc$symbol)

  #Remove regions outside of GTEx eQTL annotation borders
  abc <- abc %>%
    left_join(gtexGtf, by = "geneSymbol",relationship = "many-to-many") %>%
    filter(end > (TSS_PosGtf - 1000000) & end < (TSS_PosGtf + 1000000))
  if(exists("abcFull")){
    abcFull <- rbind(abcFull,abc)   
  }else{
    abcFull <- abc
  }
}
abcFullFreq <- data.frame(table(abcFull$peakEnsemblID))
abcFullFreq <- abcFullFreq[abcFullFreq$Freq == nrow(inputList),]
abcFull <- abcFull[abcFull$peakEnsemblID %in% as.character(abcFullFreq$Var1),]

##
##Read in single model
for(i in 1:nrow(inputList)){
  abc <- fread(inputList$samplePath[i],data.table = FALSE, header = TRUE)
  colnames(abc) <- c("chrom","start","end","ensemblID","symbol","peakID","signalValue","contact","adaptedActivity","scaledContact","intergenicScore","TSS_dist","ABC_Score")
  abc <- makeGRangesFromDataFrame(df = abc,seqnames.field = "chrom",start.field = "start",end.field = "end",starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  abc <- abc[abc$TSS_dist >  opt$minDist,]
  abc$ensemblID <- gsub("\\..*","",abc$ensemblID)
  abc$peakEnsemblID <- paste0(abc$ensemblID,"_",abc$peakID)
  abc$geneSymbol <- paste0(abc$ensemblID,"_",abc$symbol)

  #Remove non-shared enhancers
  abc <- abc[abc$peakEnsemblID %in% abcFull$peakEnsemblID,]
  abc <- abc[abc$ensemblID %in% gtex$ensemblID,]
  abc$modelName <- inputList$modelName[i]
  abc$rankABC <- rank(x = -abc$ABC_Score,ties.method = "max")
  abc$perc <- abc$rankABC / max(abc$rankABC)
    
  #Match to caviar eQTLs (removing duplicates)
  caviarShort <- caviar[caviar$ensemblID_CAVIAR %in% abc$ensemblID,]
  caviarShort <- makeGRangesFromDataFrame(df = caviarShort,seqnames.field = "CHROM",start.field = "POS",end.field = "POS",starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE)
  abcCaviarShort <- data.frame(join_overlap_inner(caviarShort,abc))
  abcCaviarShort <- abcCaviarShort[order(abcCaviarShort$rankABC),]
  abcCaviarShort$eQTLMatch <- FALSE
  abcCaviarShort$eQTLMatch[abcCaviarShort$ensemblID == abcCaviarShort$ensemblID_CAVIAR] <- TRUE
  abcCaviarShort <- abcCaviarShort[!duplicated(abcCaviarShort[,'peakEnsemblID']),]
  abcCaviarShort <- abcCaviarShort[abcCaviarShort$eQTLMatch == TRUE,]
  recall <- 0
  abcCaviarShort$recall <- recall
  for(j in 1:nrow(abcCaviarShort)){
    if(abcCaviarShort$eQTLMatch[j] == TRUE){
    recall <- recall + 1/(nrow(abcCaviarShort[abcCaviarShort$eQTLMatch == TRUE,]))
   }
    abcCaviarShort$recall[j] <- recall
  }

  pltRec <- c(pltRec,abcCaviarShort$recall)
  pltPerc <- c(pltPerc,abcCaviarShort$perc)
  pltAbc <- c(pltAbc,abcCaviarShort$ABC_Score)
  pltRank <- c(pltRank,abcCaviarShort$rankABC)
  pltName <- c(pltName,abcCaviarShort$modelName)
  plteQTL <- c(plteQTL,abcCaviarShort$eQTL)
  plteQTLMatch <- c(plteQTLMatch,abcCaviarShort$eQTLMatch)
  pltEnsemblID <- c(pltEnsemblID,abcCaviarShort$ensemblID)
  pltpeakEnsemblID <- c(pltpeakEnsemblID,abcCaviarShort$peakEnsemblID)
  }
finalDF <- data.frame(recall = pltRec, abc = pltAbc, rankABC = pltRank, perc = pltPerc, model = pltName, eQTL = plteQTL,eQTLMatch = plteQTLMatch, ensemblID = pltEnsemblID, peakEnsemblID = pltpeakEnsemblID)
finalDF <- finalDF[order(finalDF$eQTL, finalDF$ensemblID),]
write.table(finalDF, paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_finalDF.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Write and plot paired differences + their statistics + number of "+","0" and "-"
finalDF$ID <- paste0(finalDF$peakEnsemblID,"_",finalDF$eQTL)
differences <- finalDF %>%
  group_by(ID)  %>%
  reframe(Diff = combn(rankABC,2,diff))
write.table(differences, paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_differences.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(data.frame(quantile = names(quantile(differences$Diff)), value = quantile(differences$Diff)), paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_differencesQuantiles.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(data.frame(table(sign(differences$Diff))), paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_differencesFrac.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Sign test
stat.test <- finalDF %>%
  sign_test(rankABC~model, detailed= TRUE)
if(quantile(differences$Diff,.5) > 0 & stat.test$estimate < 0){
	stat.test$estimate <- abs(stat.test$estimate)
	conf.low <- -1 * stat.test$conf.high
	conf.high <- -1 * stat.test$conf.low
	stat.test$conf.low <- conf.low
	stat.test$conf.high <- conf.high
}else if(quantile(differences$Diff,.5) < 0 & stat.test$estimate > 0){
	stat.test$estimate <- -(stat.test$estimate)
	conf.low <- -1 * stat.test$conf.high
	conf.high <- -1 * stat.test$conf.low
	stat.test$conf.low <- conf.low
	stat.test$conf.high <- conf.high
}
write.table(stat.test, paste0("../../data/GTEx/evalGTEx/ev/",opt$outDir,"/",opt$compMod,"_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_signTable.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Make plots

#Full Boxplot
plotBox(diffTable = differences, plotMode = "Full", capVal = "")

#Boxplot - Outlier Capped
upperY <- quantile(differences$Diff, 1 - opt$cap)
lowerY <- quantile(differences$Diff,opt$cap)
differencesCapped <- differences
differencesCapped$Diff[differencesCapped$Diff > upperY] <- upperY
differencesCapped$Diff[differencesCapped$Diff < lowerY] <- lowerY
plotBox(diffTable = differencesCapped, plotMode = "Capped", capVal = opt$cap)
