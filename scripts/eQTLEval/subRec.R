#This script generates rank-recall plots for ABC results + bed and IGV files of promoters and thresholded enhancers.

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
library(optparse)
library(tidyverse)
library(plyr)
library(dplyr)

###
###
###Arguments
option_list = list(
  
  make_option("--abcPath", type="character", default=NULL,
              help="Path to ABC file"),
  
  make_option("--minDist", type="numeric", default=NULL,
              help="Define minimum distance to TSS to define an enhancer. Set to -1 to consider all promoters"),
  
  make_option("--gtexCut", type="numeric", default=NULL,
              help="Minimum TPM value of GTEx expression for a gene to be included"),
			  
  make_option("--tissue", type="character", default=NULL,
              help="Tissue to filter TPM and caviar file on"),

  make_option("--promoterEnhancerFile", type="character", default=NULL,
              help="Path to EnhancerList.txt. Required for promoter definition"),

  make_option("--promQuant", type="numeric", default=NULL,
              help="Top quantile of promoter DNAse-seq to keep e.g. 0.6"),
  
  make_option("--caviarThres", type="numeric", default=NULL,
              help="eQTL cut-off"),
  
  make_option("--recCut", type="numeric", default=NULL,
              help="Desired eQTL recall cut-off"),
  
  make_option("--promoterFile", type="character", default=NULL,
              help="File with promoter positions (as in ABC). Provide without path to folder"),
   
  make_option("--outDir", type="character", default=NULL,
              help="Path to output directory"),
			  
  make_option("--addTag", type="character", default=NULL,
              help="Path to output directory")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$abcPath) | is.null(opt$minDist) | is.null(opt$gtexCut) |
    is.null(opt$tissue) | is.null(opt$promoterEnhancerFile) | is.null(opt$promQuant) | is.null(opt$caviarThres) |
    is.null(opt$recCut) | is.null(opt$outDir) | is.null(opt$addTag)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Main

##
##Make output directory
dir.create(paste0("../../data/GTEx/recSub/sub/",opt$outDir), showWarnings = FALSE)

##
##Read in gtex genes and filter on expression
gtex <- fread(input = "../../data/GTEx/raw/GTExTPM.gct",data.table = FALSE)
gtex <- gtex[,colnames(gtex) %in% c("Name","Description",opt$tissue)]
colnames(gtex) <- c("ensemblID","symbol","tissueTPM")
gtex$ensemblID <- gsub("\\..*","",gtex$ensemblID)
gtex <- gtex[gtex$tissueTPM >= opt$gtexCut,]

##
##Read in promoter open chrom information and filter on
promoters <- fread(opt$promoterEnhancerFile, header = TRUE, data.table = FALSE)
promoters <- promoters[,colnames(promoters) %in% c("chr","start","end","normalized_dhs","isPromoterElement","promoterSymbol")]
colnames(promoters)[colnames(promoters) == "chr"] <- "seqnames"
colnames(promoters)[colnames(promoters) == "promoterSymbol"] <- "ensemblID"
promoters <- promoters[promoters$ensemblID %in% gtex$ensemblID,]
promoters <- promoters[!(promoters$seqnames %in% c("chrX","chrY","chrM")),]
promoters <- promoters[promoters$isPromoterElement == TRUE,]
promoters <- promoters[promoters$normalized_dhs >= quantile(promoters$normalized_dhs,(1-opt$promQuant)),]
promoters$seqnames <- substr(promoters$seqnames,4,5)
promoters <- promoters %>%
	mutate(ensemblID = strsplit(as.character(ensemblID), ",")) %>%
	unnest(ensemblID)
promoters <- promoters[promoters$ensemblID != "",]
promoters$ensemblID <- gsub(pattern = ".*\\_",replacement = "", x = promoters$ensemblID)

#Form overlap of expressed and open chrom-promoter genes
gtexOpen <- gtex[gtex$ensemblID %in% promoters$ensemblID,]

##
##Process caviar data
caviar <- fread(input = "../../data/GTEx/raw/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz", data.table = FALSE)
colnames(caviar)[2] <- "ensemblID_CAVIAR"
caviar$ensemblID_CAVIAR <- gsub("\\..*","",caviar$ensemblID_CAVIAR)
caviar <- caviar[caviar$Probability >= opt$caviarThres,]
caviar <- caviar[caviar$TISSUE == opt$tissue,]
caviar <- caviar[caviar$ensemblID_CAVIAR %in% gtexOpen$ensemblID,]

##
##Process gtex gtf data
gtexGtf <- fread(input = "../../data/GTEx/raw/gencodeGTEx.gtf", data.table = FALSE, header = FALSE)
gtexGtf <- gtexGtf[gtexGtf$V3 == "gene",]
gtexGtfColsplit <- colsplit(gtexGtf$V9,pattern = "; ",names = c("ensemblID","transcriptID","geneType","symbol","remaining"))
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
##Read in maABC probe-model
abc <- fread(opt$abcPath,data.table = FALSE, select = c(1:6,12:13))
colnames(abc) <- c("chrom","start","end","ensemblID","symbol","peakID","TSS_dist","ABC_Score")

#Remove variants outside of GTEx eQTL detection borders
abc$ensemblID <- gsub("\\..*","",abc$ensemblID)
abc <- abc[abc$ensemblID %in% gtexOpen$ensemblID,]
abc$geneSymbol <- paste0(abc$ensemblID,"_",abc$symbol)
abc <- abc %>%
	left_join(gtexGtf, by = "geneSymbol") %>%
	filter(end > (TSS_PosGtf - 1000000) & end < (TSS_PosGtf + 1000000))

#Match maABC to caviar
abc <- abc[abc$TSS_dist >  opt$minDist,]
abc$rankABC <- rank(x = -abc$ABC_Score,ties.method = "max")
abc <- makeGRangesFromDataFrame(df = abc,seqnames.field = "chrom",start.field = "start",end.field = "end",starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
caviarShort <- caviar[caviar$ensemblID_CAVIAR %in% abc$ensemblID,]
caviarShort <- makeGRangesFromDataFrame(df = caviarShort,seqnames.field = "CHROM",start.field = "POS",end.field = "POS",starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE)
abcCaviarShort <- data.frame(join_overlap_inner(caviarShort,abc))
abcCaviarShort <- abcCaviarShort[order(abcCaviarShort$rankABC),]
abcCaviarShort$eQTLMatch <- FALSE
abcCaviarShort$eQTLMatch[abcCaviarShort$ensemblID == abcCaviarShort$ensemblID_CAVIAR] <- TRUE
abcCaviarShort <- abcCaviarShort[abcCaviarShort$eQTLMatch == TRUE,]
recall <- 0
abcCaviarShort$recall <- recall
for(j in 1:nrow(abcCaviarShort)){
  if(abcCaviarShort$eQTLMatch[j] == TRUE){
    recall <- recall + 1/(nrow(abcCaviarShort[abcCaviarShort$eQTLMatch == TRUE,]))
  }
  abcCaviarShort$recall[j] <- recall
}

##
##Plot curve
plotDF <- abcCaviarShort[abcCaviarShort$eQTLMatch == TRUE,]
p <- ggplot(plotDF, aes(rankABC, recall)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  geom_vline(xintercept = max(abcCaviarShort$rankABC[abcCaviarShort$recall <= opt$recCut]), linetype="dashed", 
             color = "green", size=0.5)+
  geom_hline(yintercept = opt$recCut, linetype="dashed", 
             color = "red", size=0.5)+			 
  xlim(0,round_any(max(abcCaviarShort$rankABC[abcCaviarShort$recall <= opt$recCut]) + 10000, 10000 + 1, f = ceiling) )+
  xlab("rank ABC") + ylab("HC-eQTL recall")
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"plotCurves_minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,".png"),p,width = 6,height = 4.38,dpi = 600)
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"plotCurves_minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,".pdf"),p,width = 6,height = 4.38)

##
##Write table and igv files (at desired cut-off)

#Promoters
write.table(promoters[,colnames(promoters) %in% c("seqnames","start","end","ensemblID")],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviar",opt$recCut,"_promoters.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Reduce promoters to target genes
targetGenes <- fread(input = "../../data/GTEx/recSub/targetGenes/targets.txt", data.table = FALSE, header = FALSE)
write.table(promoters[promoters$ensemblID %in% targetGenes$V1,colnames(promoters) %in% c("seqnames","start","end","ensemblID")],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviar",opt$recCut,"_promotersTargets.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Enhancer at defined ABC-eQTL recall cut-off
abc <- data.frame(abc)
abc$strand <- NULL
abc$width <- NULL
minABC <- min(abcCaviarShort$ABC_Score[abcCaviarShort$recall <= opt$recCut])
abc <- abc[abc$ABC_Score >= minABC,]

#Make link table
abcLink <- merge(abc,promoters, by = c("seqnames","ensemblID"))

#Make and write final enhancer table
abc <- data.frame(seqnames = abc[,colnames(abc) == "seqnames"],
                      start = abc[,colnames(abc) == "start"],
                      end = abc[,colnames(abc) == "end"],
                      ensemblID = abc[,colnames(abc) == "ensemblID"])
write.table(abc,paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_enhancers.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Reduce enhancers to target genes
write.table(abc[abc$ensemblID %in% targetGenes$V1,],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_enhancersTargets.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

##
##Go on with the link table
abcLink <- data.frame(chrProm = abcLink[,colnames(abcLink) == "seqnames"],
                  startProm = abcLink[,colnames(abcLink) == "startProm"],
                  endProm = abcLink[,colnames(abcLink) == "endProm"],
                  chrEnh = abcLink[,colnames(abcLink) == "seqnames"],
                  startEnh = abcLink[,colnames(abcLink) == "start"],
                  endEnh = abcLink[,colnames(abcLink) == "end"],
                  ensemblID = abcLink[,colnames(abcLink) == "ensemblID"],
                  score = abcLink[,colnames(abcLink) == "ABC_Score"],
                  strandProm = abcLink[,colnames(abcLink) == "strand"],
                  strandEnh = "*")
write.table(abcLink,paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_links.bedpe"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Reduce link table to defined target genes
write.table(abcLink[abcLink$ensemblID %in% targetGenes$V1,],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$addTag,"minDist",opt$minDist,"_TPMCut",opt$gtexCut,"_promQuant",opt$promQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_linksTargets.bedpe"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)
