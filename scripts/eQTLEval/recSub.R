#This script generates rank-recall plots for caABC results + bed and IGV files of promoters and thresholded enhancers.

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
library(plyr)
library(tidyverse)
library(dplyr)

###
###
###Arguments
option_list = list(
  
  make_option("--caABCpath", type="character", default=NULL,
              help="Path to caABC file"),

  make_option("--sample", type="character", default=NULL,
              help="Sample name"),
 
  make_option("--minDist", type="numeric", default=NULL,
              help="Define minimum distance to TSS to define an enhancer. Set to -1 to consider all promoters"),
  
  make_option("--gtexCut", type="numeric", default=NULL,
              help="Minimum TPM value of GTEx expression for a gene to be included"),
			  
  make_option("--tissue", type="character", default=NULL,
              help="Tissue to filter TPM and caviar file on"),

  make_option("--caviarThres", type="numeric", default=NULL,
              help="eQTL cut-off"),

  make_option("--recCut", type="numeric", default=NULL,
              help="Desired eQTL recall cut-off"),

  make_option("--promoterEnhancerFile", type="character", default=NULL,
              help="Path to EnhancerList.txt. Required for promoter definition"),

  make_option("--topPromQuant", type="numeric", default=NULL,
              help="Top quantile of promoter DNaseSeq to keep e.g. 0.85"),

  make_option("--cGtf", type="character", default=NULL,
              help="Path to canonical GTF file"),
   
  make_option("--outDir", type="character", default=NULL,
              help="Path to output directory")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$caABCpath) | is.null(opt$sample) | is.null(opt$minDist) | is.null(opt$gtexCut) |
    is.null(opt$tissue) | is.null(opt$caviarThres) | is.null(opt$recCut) | is.null(opt$promoterEnhancerFile) | is.null(opt$topPromQuant) |
    is.null(opt$outDir) | is.null(opt$cGtf)) {
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
##Read in and filter for gtex expression
gtex <- fread(input = "../../data/GTEx/raw/GTExTPM.gct",data.table = FALSE)
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
caviar$CHROM <- paste0("chr", caviar$CHROM)

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
##Read in caABC sample-model
caABC <- fread(opt$caABCpath,data.table = FALSE, select = c(1:6,12:13), header = TRUE)
colnames(caABC) <- c("chrom","start","end","ensemblID","symbol","peakID","TSS_dist","ABC_Score")
caABC$chrom <- paste0("chr", caABC$chrom)

#Remove variants outside of GTEx eQTL detection borders
caABC$ensemblID <- gsub("\\..*","",caABC$ensemblID)

caABC$geneSymbol <- paste0(caABC$ensemblID,"_",caABC$symbol)
caABC <- caABC %>%
	left_join(gtexGtf, by = "geneSymbol",relationship = "many-to-many") %>%
	filter(end > (TSS_PosGtf - 1000000) & end < (TSS_PosGtf + 1000000))

##
##Read in open chrom information and filter on
promotersDist <- caABC[caABC$TSS_dist == 0,c(1,2,3,4,6)]
promoters <- fread(opt$promoterEnhancerFile, header = TRUE, data.table = FALSE)
promoters <- promoters[,colnames(promoters) %in% c("chr","start","end","normalized_dhs","isPromoterElement","promoterSymbol")]
colnames(promoters)[colnames(promoters) == "chr"] <- "chrom"
colnames(promoters)[colnames(promoters) == "promoterSymbol"] <- "ensemblID"
colnames(promoters)[colnames(promoters) == "normalized_dhs"] <- "DNaseSeq_signal"
promoters <- promoters[!(promoters$chrom %in% c("chrX","chrY","chrM")),]
promoters <- promoters[promoters$isPromoterElement == TRUE,]
promoters <- promoters %>%
	mutate(ensemblID = strsplit(as.character(ensemblID), ",")) %>%
	unnest(ensemblID)
promoters <- promoters[promoters$ensemblID != "",]
promoters$ensemblID <- gsub(pattern = ".*\\_",replacement = "", x = promoters$ensemblID)
promoters <- promoters[promoters$ensemblID %in% gtex$ensemblID,]
promoters <- merge(promoters,promotersDist, by = c("chrom","start","end","ensemblID"))

#Make hist of log2 +1 dnaseseq signal and show distribution + threshold
p <- ggplot(promoters, aes(x=log2(DNaseSeq_signal +1))) + geom_histogram(bins=50) +
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
		xlab("log2(DNase-seq signal +1)")+
  geom_vline(xintercept = log2(quantile(promoters$DNaseSeq_signal,(1-opt$topPromQuant)) +1), linetype="dashed", 
             color = "black", linewidth=0.5)
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","plotCurves_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuantLine",opt$topPromQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,"_promAct.png"),p,width = 6,height = 4.38,dpi = 600)
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","plotCurves_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuantLine",opt$topPromQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,"_promAct.pdf"),p,width = 6,height = 4.38)
promoters <- promoters[promoters$DNaseSeq_signal >= quantile(promoters$DNaseSeq_signal,(1-opt$topPromQuant)),]

#Form overlap of expressed and open chrom-promoter genes and subset data
gtexOpen <- gtex[gtex$ensemblID %in% promoters$ensemblID,]
caviar <- caviar[caviar$ensemblID_CAVIAR %in% gtexOpen$ensemblID,]
caABC <- caABC[caABC$ensemblID %in% gtexOpen$ensemblID,]

##
##Match caABC to caviar (removing duplicates)
caABC <- caABC[caABC$TSS_dist >  opt$minDist,]
caABC$rank_caABC <- rank(x = -caABC$ABC_Score,ties.method = "max")
caABC <- makeGRangesFromDataFrame(df = caABC,seqnames.field = "chrom",start.field = "start",end.field = "end",starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
caviarShort <- caviar[caviar$ensemblID_CAVIAR %in% caABC$ensemblID,]
caviarShort <- makeGRangesFromDataFrame(df = caviarShort,seqnames.field = "CHROM",start.field = "POS",end.field = "POS",starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE)
caABCcaviarShort <- data.frame(join_overlap_inner(caviarShort,caABC))
caABCcaviarShort <- caABCcaviarShort[order(caABCcaviarShort$rank_caABC),]
caABCcaviarShort$eQTLMatch <- FALSE
caABCcaviarShort$eQTLMatch[caABCcaviarShort$ensemblID == caABCcaviarShort$ensemblID_CAVIAR] <- TRUE
caABCcaviarShort <- caABCcaviarShort[caABCcaviarShort$eQTLMatch == TRUE,]
caABCcaviarShort$peakEnsemblID <- paste0(caABCcaviarShort$ensemblID,"_",caABCcaviarShort$peakID)
caABCcaviarShort <- caABCcaviarShort[!duplicated(caABCcaviarShort[,'peakEnsemblID']),]
recall <- 0
caABCcaviarShort$recall <- recall
for(j in 1:nrow(caABCcaviarShort)){
  if(caABCcaviarShort$eQTLMatch[j] == TRUE){
    recall <- recall + 1/(nrow(caABCcaviarShort[caABCcaviarShort$eQTLMatch == TRUE,]))
  }
  caABCcaviarShort$recall[j] <- recall
}

##
##Plot curve
plotDF <- caABCcaviarShort[caABCcaviarShort$eQTLMatch == TRUE,]
p <- ggplot(plotDF, aes(rank_caABC, recall)) +
  geom_line(linewidth = 0.5) +
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
  geom_vline(xintercept = max(caABCcaviarShort$rank_caABC[caABCcaviarShort$recall <= opt$recCut]), linetype="dashed", 
             color = "green", linewidth=0.5)+
  geom_hline(yintercept = opt$recCut, linetype="dashed", 
             color = "red", linewidth=0.5)+			 
  xlim(0,round_any(max(caABCcaviarShort$rank_caABC[caABCcaviarShort$recall <= opt$recCut]) + 10000, 10000 + 1, f = ceiling) )+
  xlab("rank caABC") + ylab("HC-eQTL recall")
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","plotCurves_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,".png"),p,width = 6,height = 4.38,dpi = 600)
ggsave(filename = paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","plotCurves_minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviar",opt$caviarThres,"_recCut",opt$recCut,".pdf"),p,width = 6,height = 4.38)

##
##Write table and igv files (at desired cut-off)

#promoters
write.table(promoters[,colnames(promoters) %in% c("chrom","start","end","ensemblID")],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_recCut",opt$recCut,"_promoters.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Reduce promoters to target genes
targetGenes <- fread(input = "../../data/GTEx/recSub/targetGenes/targets.txt", data.table = FALSE, header = FALSE)
write.table(promoters[promoters$ensemblID %in% targetGenes$V1,colnames(promoters) %in% c("chrom","start","end","ensemblID")],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_recCut",opt$recCut,"_promotersTargets.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Enhancer at defined ABC-eQTL recall cut-off
caABC <- data.frame(caABC)
caABC$strand <- NULL
caABC$width <- NULL
colnames(caABC)[colnames(caABC) == "seqnames"] <- "chrom"
minABC <- min(caABCcaviarShort$ABC_Score[caABCcaviarShort$recall <= opt$recCut])
caABC <- caABC[caABC$ABC_Score >= minABC,]
caABC$start <- caABC$start -1

#Make link table
caABClink <- merge(caABC,promoters, by = c("chrom","ensemblID"))

#Make and write final enhancer table
caABC <- data.frame(chrom = caABC[,colnames(caABC) == "chrom"],
                      start = caABC[,colnames(caABC) == "start"],
                      end = caABC[,colnames(caABC) == "end"],
                      ensemblID = caABC[,colnames(caABC) == "ensemblID"])
write.table(caABC,paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_enhancers.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

#Reduce enhancers to target genes
write.table(caABC[caABC$ensemblID %in% targetGenes$V1,],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_enhancersTargets.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

##
##Go on with the link table
caABClink <- data.frame(chrProm = caABClink[,colnames(caABClink) == "chrom"],
                  startProm = caABClink[,colnames(caABClink) == "start.x"],
                  endProm = caABClink[,colnames(caABClink) == "end.x"],
                  chrEnh = caABClink[,colnames(caABClink) == "chrom"],
                  startEnh = caABClink[,colnames(caABClink) == "start.y"],
                  endEnh = caABClink[,colnames(caABClink) == "end.y"],
                  ensemblID = caABClink[,colnames(caABClink) == "ensemblID"],
                  score = caABClink[,colnames(caABClink) == "rank_caABC"],
                  strandProm = caABClink[,colnames(caABClink) == "strand"],
                  strandEnh = "*")
write.table(caABClink,paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_links.bedpe"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

##
##Reduce link table to defined target genes
write.table(caABClink[caABClink$ensemblID %in% targetGenes$V1,],paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_linksTargets.bedpe"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

##
##Write list of considered genes together with their median TPM values
gtex <- gtex[gtex$ensemblID %in% caABC$ensemblID,]
write.table(gtex, paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_gtexExp.txt"),sep = "\t",append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE)

##
##Cut Canonical-GTF file down to expressed genes

#GTF
ctGtf <- fread(opt$cGtf,data.table = FALSE, header = FALSE)
ctGtfGene <- colsplit(ctGtf$V9,pattern = "; ",names = c("ensemblID","remaining"))
ctGtfGene$ensemblID <- gsub("gene_id \"","",ctGtfGene$ensemblID)
ctGtfGene$ensemblID <- gsub("\\..*","",ctGtfGene$ensemblID)
ctGtf <- ctGtf[ctGtfGene$ensemblID %in% caABC$ensemblID,]
ctGtfGene <- ctGtfGene[ctGtfGene$ensemblID %in% caABC$ensemblID,]
ctBed <- data.frame(ctGtf$V1,ctGtf$V4 -1,ctGtf$V5, ctGtfGene$ensemblID)

#Bed
write.table(ctBed, paste0("../../data/GTEx/recSub/sub/",opt$outDir,"/",opt$sample,"_","minDist",opt$minDist,"_",opt$tissue,"_TPMCut",opt$gtexCut,"_topPromQuant",opt$topPromQuant,"_caviarCut",opt$caviarThres,"_recCut",opt$recCut,"_ctExpGene.bed"),sep = "\t",append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE)

