#This script compares cage seq signals for longest transcript and MANE transcript based TSS

###
###
###Libraries
options(scipen = 999)
set.seed(1234)
library(data.table)
library(ggplot2)
library(optparse)
library(GenomicRanges)
library(dplyr)
library(rstatix)

###
###
###Arguments
option_list = list(
  
  make_option("--maneTranscripts", type="character", default=NULL,
              help="Name of mane genes file"),
  
  make_option("--longestTranscripts", type="character", default=NULL,
              help="Name of longest genes file"),
  
  make_option("--cagePlus", type="character", default=NULL,
              help="Name of cage seq plus strand file"),
  
  make_option("--cageMinus", type="character", default=NULL,
              help="Name of cage seq minus strand file"),
			  
  make_option("--tissue", type="character", default=NULL,
              help="Tissue to filter GTEx on"),
			  
  make_option("--gtexCut", type="numeric", default=NULL,
              help="minimum TPM threshold of expression"),
  
  make_option("--minDist", type="numeric", default=5000,
              help="How far should the Longest-transcript and Mane-transcript TSS be distant from another to be considered in the analysis?"),
  
  make_option("--maxDist", type="numeric", default=250,
              help="Maximum Distance to look for the maximum cage seq signal around a TSS"),
			  
  make_option("--cap", type="numeric", default=NULL,
              help="Given cap for outliers in boxplot. Equally applied to upper/lower")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$maneTranscripts) | is.null(opt$longestTranscripts) | is.null(opt$cagePlus) | is.null(opt$cageMinus) |
    is.null(opt$tissue) | is.null(opt$gtexCut) | is.null(opt$cap)){
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Functions
plotBox <- function(diffTable, plotMode, capVal){
	p <- ggplot(diffTable, aes(x = "",y =Diff)) +geom_boxplot(outlier.shape = NA, width = 0.4) + geom_jitter(cex = 0.5,width = 0.1, alpha = 0.5)+
		xlab(NULL) + ylab("Pairwise difference of CAGE-Seq signal")+
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			legend.key=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=14))+
			scale_x_discrete()
	ggsave(filename = paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_maxDist",opt$maxDist,"_",opt$tissue,"_","TPMCut",opt$gtexCut,"_plotMode",plotMode,capVal,"_BoxDiff.png"),p,width = 3.5,height = 4.38,dpi = 600)
	ggsave(filename = paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_maxDist",opt$maxDist,"_",opt$tissue,"_","TPMCut",opt$gtexCut,"_plotMode",plotMode,capVal,"_BoxDiff.pdf"),p,width = 3.5,height = 4.38)
}

###
###
###Main

##
##Get expressed genes in GTEx
gtex <- fread("../../data/GTEx/raw/GTExTPM.gct",data.table = FALSE)
gtex <- gtex[,colnames(gtex) %in% c("Name","Description",opt$tissue)]
colnames(gtex) <- c("ensemblID","symbol","tissueTPM")
gtex$ensemblID <- gsub("\\..*","",gtex$ensemblID)
gtex <- gtex[gtex$tissueTPM >= opt$gtexCut,]
gtex$symGene <- paste0(gtex$symbol,"_",gtex$ensemblID)

##
##Subset to TSS (considering strand) and to expressed genes
MLT <- fread(paste0("../../data/ref/",opt$maneTranscripts), select = c(1:4,6), data.table = FALSE)
colnames(MLT) <- c("chr_MLT","startMLT","endMLT","symGene","strand")
MLT <- MLT[MLT$symGene %in% gtex$symGene,]
MLT <- MLT[order(MLT$symGene),]
lT <- fread(paste0("../../data/ref/",opt$longestTranscripts), select = c(1:4,6), data.table = FALSE)
colnames(lT) <- c("chr_lT","startLT","endLT","symGene","strand")
lT <- lT[lT$symGene %in% gtex$symGene,]
lT <- lT[order(lT$symGene),]

mergedT <- merge(MLT, lT, by = c("symGene","strand"))
mergedT <- mergedT[((abs(mergedT$startMLT - mergedT$startLT) > opt$minDist) & mergedT$strand == "+") | ((abs(mergedT$endMLT - mergedT$endLT) > opt$minDist) & mergedT$strand == "-"),]

mergedT$endMLT[mergedT$strand == "+"] <- mergedT$startMLT[mergedT$strand == "+"] + opt$maxDist
mergedT$endLT[mergedT$strand == "+"] <- mergedT$startLT[mergedT$strand == "+"] + opt$maxDist
mergedT$startMLT[mergedT$strand == "+"] <- mergedT$startMLT[mergedT$strand == "+"] - opt$maxDist
mergedT$startLT[mergedT$strand == "+"] <- mergedT$startLT[mergedT$strand == "+"] - opt$maxDist

mergedT$startMLT[mergedT$strand == "-"] <- mergedT$endMLT[mergedT$strand == "-"] - opt$maxDist
mergedT$startLT[mergedT$strand == "-"] <- mergedT$endLT[mergedT$strand == "-"] - opt$maxDist
mergedT$endMLT[mergedT$strand == "-"] <- mergedT$endMLT[mergedT$strand == "-"] + opt$maxDist
mergedT$endLT[mergedT$strand == "-"] <- mergedT$endLT[mergedT$strand == "-"] + opt$maxDist

##
##Annotate maximum CAGE-Seq signal. If no cage seq signal present --> set to 0
cagePlus <- fread(paste0("../../data/cage/signals/",opt$cagePlus), data.table = FALSE)
cagePlus <- makeGRangesFromDataFrame(cagePlus,seqnames.field = "V1",start.field = "V2", end.field = "V3",keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)
cageMinus <- fread(paste0("../../data/cage/signals/",opt$cageMinus), data.table = FALSE)
cageMinus <- makeGRangesFromDataFrame(cageMinus,seqnames.field = "V1",start.field = "V2", end.field = "V3",keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)

maxMane <- vector("list",nrow(mergedT))
maxLongest <- vector("list",nrow(mergedT))
for(i in 1:nrow(mergedT)){
  maneTss=GRanges(seqnames=mergedT$chr_MLT[i],
                  ranges=IRanges(start = mergedT$startMLT[i], mergedT$endMLT[i]))
  longestTss=GRanges(seqnames=mergedT$chr_lT[i],
                  ranges=IRanges(start = mergedT$startLT[i], mergedT$endLT[i]))  
  if(mergedT$strand[i] == "+"){
  
    #Plus MANE
    cagePlusMane <- subsetByOverlaps(cagePlus, maneTss)
    cagePlusMane <- max(cagePlusMane$V4)
    if((is.infinite(cagePlusMane)) == TRUE){
      maxMane[i] <- 0
    }else if((!(is.infinite(cagePlusMane))) == TRUE){
      maxMane[i] <- cagePlusMane
    }
	
    #Plus longest
    cagePlusLongest <- subsetByOverlaps(cagePlus, longestTss)
    cagePlusLongest <- max(cagePlusLongest$V4)
    if((is.infinite(cagePlusLongest)) == TRUE){
      maxLongest[i] <- 0
    }else if((!(is.infinite(cagePlusLongest))) == TRUE){
      maxLongest[i] <- cagePlusLongest
    }
  }else if(mergedT$strand[i] == "-"){
  
    #Minus MANE
    cageMinusMane <- subsetByOverlaps(cageMinus, maneTss)
    cageMinusMane <- max(cageMinusMane$V4)
    if(is.infinite(cageMinusMane)){
      maxMane[i] <- 0
    }else if(!(is.infinite(cageMinusMane))){
      maxMane[i] <- cageMinusMane
    }
	
    #Minus longest
    cageMinusLongest <- subsetByOverlaps(cageMinus, longestTss)
    cageMinusLongest <- max(cageMinusLongest$V4)
    if(is.infinite(cageMinusLongest)){
      maxLongest[i] <- 0
    }else if(!(is.infinite(cageMinusLongest))){
      maxLongest[i] <- cageMinusLongest
    }
  }
}

mergedT$maxMane <- unlist(maxMane)
mergedT$maxLongest <- unlist(maxLongest)
mergedT$Diff <- mergedT$maxMane - mergedT$maxLongest

#0-index
mergedT$startMLT <- mergedT$startMLT -1
mergedT$startLT <- mergedT$startLT -1
mergedT$endMLT <- mergedT$endMLT -1
mergedT$endLT <- mergedT$endLT -1
write.table(mergedT,paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_maxDist",opt$maxDist,"_",opt$tissue,"_","TPMCut",opt$gtexCut,"_cageTableFull.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Write number of "+","0" and "-" + median of differences
write.table(data.frame(table(sign(mergedT$Diff))),paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_recallCut",opt$recallCut,"_differencesFrac.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)
write.table(median(mergedT$Diff),paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_caviarCut",opt$caviarThres,"_recallCut",opt$recallCut,"_differencesMedian.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

##
##Make Sign test
differences <- data.frame(symGene = c(mergedT$symGene,mergedT$symGene),model = c(rep("MANE",nrow(mergedT)),rep("Longest",nrow(mergedT))), cageSignal = c(mergedT$maxMane, mergedT$maxLongest))
stat.test <- differences %>%
  sign_test(cageSignal~model, detailed= TRUE)
write.table(stat.test,paste0("../../data/cage/cageResults/cage_minDist",opt$minDist,"_maxDist",opt$maxDist,"_",opt$tissue,"_","TPMCut",opt$gtexCut,"_signTable.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Make boxplots

#Full
plotBox(diffTable = mergedT, plotMode = "Full", capVal = "")

#Outlier capped
upperY <- quantile(mergedT$Diff, 1 - opt$cap)
lowerY <- quantile(mergedT$Diff,opt$cap)
mergedT$Diff[mergedT$Diff > upperY] <- upperY
mergedT$Diff[mergedT$Diff < lowerY] <- lowerY
plotBox(diffTable = mergedT, plotMode = "Capped", capVal = opt$cap)
