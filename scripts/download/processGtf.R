#This scripts generates input Transcripts, TSS and a gtf-file for ABC, M-AD-ABC and gABC
###
###
###Libraries
set.seed(1234)
options(scipen=999)
library(reshape2)
library(data.table)
library(optparse)

###
###
###External arguments
option_list = list(
  make_option("--ref", type="character", default=NULL, 
              help="path to reference directory"),
 
  make_option("--gencodeName", type="character", default=NULL, 
              help="Nime of the gencode file without ending"))  
			  
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$ref) | is.null(opt$gencodeName)) {
  print("You need to specify all arguments.")
  print_help(opt_parser)
  quit()
}

###
###
###Functions

##
##Function writes at the end results in ABC-required gene and tss files.
writeGeneFiles <- function(gencodeData, transcriptType){
	abcGenes <- gencodeData[,c(1,4,5,9,7)]
	abcGenes <- colsplit(gencodeData$V9, "; ", c("gene_id","gene_type", "gene_name"))
	abcGenes <- data.frame(chr = gencodeData$V1,start = gencodeData$V4, end = gencodeData$V5, name = abcGenes[,3],
						score = 0, strand = gencodeData$V7,Ensembl_ID = abcGenes[,1], gene_type = abcGenes[,2])
	
	abcGenes$Ensembl_ID <- gsub("\\..*","",abcGenes$Ensembl_ID)
	abcGenes$Ensembl_ID <- gsub("gene_id \"","",abcGenes$Ensembl_ID)
	abcGenes$Ensembl_ID <- gsub("\"","",abcGenes$Ensembl_ID)
	
	abcGenes$name <- gsub("\\..*","",abcGenes$name)
	abcGenes$name <- gsub("gene_name \"","",abcGenes$name)
	abcGenes$name <- gsub("\"","",abcGenes$name)
	abcGenes$name <- gsub(";","", abcGenes$name)
	abcGenes$name <- paste0(abcGenes$name,"_",abcGenes$Ensembl_ID)
	abcGenes$gene_type <- gsub("\\..*","",abcGenes$gene_type)
	abcGenes$gene_type <- gsub("gene_type \"","",abcGenes$gene_type)
	abcGenes$gene_type <- gsub("\"","",abcGenes$gene_type)
	
	write.table(x = abcGenes,file = paste0(opt$ref,opt$gencodeName,"_",transcriptType,"_gene.bed"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)
	
	abcGenes$end[abcGenes$strand == "+"] <- abcGenes$start[abcGenes$strand == "+"] +250
	abcGenes$start[abcGenes$strand == "+"] <- abcGenes$start[abcGenes$strand == "+"] -250
	abcGenes$start[abcGenes$strand == "-"] <- abcGenes$end[abcGenes$strand == "-"] -250
	abcGenes$end[abcGenes$strand == "-"] <- abcGenes$end[abcGenes$strand == "-"] +250
	write.table(x = abcGenes,file = paste0(opt$ref,opt$gencodeName,"_",transcriptType,"_tss.bed"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

###
###
###Main

##
##Process MANE
mane_raw <- read.delim(paste0(opt$ref,"mane.bed"),check.names = FALSE)
mane_raw <- cbind(mane_raw[,c(1,2,3,19)],0,mane_raw[,c(6,18,13,20,25)])
colnames(mane_raw) <- c("chrom","start","end","symbol","score","strand","Ensembl_ID","Ensembl_transcript","gene_type","maneStat")
mane_raw <- mane_raw[mane_raw$maneStat == "MANE Select",]
mane_raw$maneStat <- NULL
mane_raw$Ensembl_ID <- gsub("\\..*","",mane_raw$Ensembl_ID)
mane_raw$Ensembl_transcript <- gsub("\\..*","",mane_raw$Ensembl_transcript)
mane_raw$chrom <- gsub("chr","",mane_raw$chrom)
mane_raw$chrom[mane_raw$chrom == "X"] <- "23"
mane_raw$chrom[mane_raw$chrom == "Y"] <- "24"
mane_raw$chrom <- as.numeric(mane_raw$chrom)
mane_raw <- mane_raw[order(mane_raw$chrom,mane_raw$start),]
mane_raw$chrom <- paste0("chr",mane_raw$chrom)
mane_raw <- mane_raw[mane_raw$chrom != "chrNA",]
mane_raw$chrom[mane_raw$chrom == "chr23"] <- "chrX"
mane_raw$chrom[mane_raw$chrom == "chr24"] <- "chrY"
annoTable <- mane_raw[,c("Ensembl_ID","Ensembl_transcript")]
annoTable$geneTransPair <- paste0(mane_raw$Ensembl_ID,"_",mane_raw$Ensembl_transcript)

##
###Subset gencode to MANE transcripts where possible

gencode <- fread(paste0(opt$ref,opt$gencodeName,".gtf"),data.table = FALSE)
colnames(gencode)[1] <- "V1"

#make 0-based
gencode$V4 <- gencode$V4 -1

##
##Transcripts
gencodeTranscript <- gencode[gencode$V3 == "transcript",]
gencode_transcripts <- colsplit(gencodeTranscript$V9, "; ", c("gene_id","transcript_id","gene_type","gene_id","remaining"))

geneTranscript <- data.frame(Ensembl_transcript = gencode_transcripts$transcript_id,Ensembl_ID = gencode_transcripts$gene_id)

geneTranscript$Ensembl_transcript <- gsub("\\..*","",geneTranscript$Ensembl_transcript)
geneTranscript$Ensembl_transcript <- gsub("transcript_id ", "", geneTranscript$Ensembl_transcript)
geneTranscript$Ensembl_transcript <- gsub("\"", "", geneTranscript$Ensembl_transcript)
geneTranscript$Ensembl_ID <- gsub("\\..*","",geneTranscript$Ensembl_ID)
geneTranscript$Ensembl_ID <- gsub("gene_id ", "", geneTranscript$Ensembl_ID)
geneTranscript$Ensembl_ID <- gsub("\"", "", geneTranscript$Ensembl_ID)

geneTranscript$geneTransPair <- paste0(geneTranscript$Ensembl_ID,"_",geneTranscript$Ensembl_transcript)
geneTranscript$index <- 1:nrow(geneTranscript)
annoTableMerge <- merge(annoTable,geneTranscript, by = c("Ensembl_ID","Ensembl_transcript","geneTransPair"))
annoGeneExclude <- unique(geneTranscript$Ensembl_ID[geneTranscript$Ensembl_ID %in% annoTableMerge$Ensembl_ID])
gencodeTranscript <- gencodeTranscript[annoTableMerge$index,]

gencode_transcripts <- colsplit(gencodeTranscript$V9, "; ", c("gene_id","transript_id","gene_type","gene_name","remaining"))
gencodeTranscript$V9 <- paste0(gencode_transcripts$gene_id,"; ",gencode_transcripts$gene_type,"; ",gencode_transcripts$gene_name,";")
gencodeTranscript$V3 <- gsub("transcript", "gene",gencodeTranscript$V3)

##
##Remaining genes
gencodeGene <- gencode[gencode$V3 == "gene",]
gencode_genes <- colsplit(gencodeGene$V9, "; ", c("gene_id","gene_type", "gene_name","remaining"))
gene_id <- gsub("\\..*","",gencode_genes$gene_id)
gene_id <- gsub("gene_id \"","",gene_id)
gene_id <- gsub("\"","",gene_id)

gencodeGene$V9 <- paste0(gencode_genes$gene_id,"; ",gencode_genes$gene_type,"; ",gencode_genes$gene_name,";")

#Save for LT 
gencodeGeneLT <- gencodeGene
gencodeGeneLT <- gencodeGeneLT[order(gencodeGeneLT$V1,gencodeGeneLT$V4,gencodeGeneLT$V5),]
write.table(x = gencodeGeneLT,file = paste0(opt$ref,opt$gencodeName,"_LT.gtf"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)

#Combine transcript and gene
gencodeGene <- gencodeGene[!(gene_id %in% annoGeneExclude),]
gencode <- rbind(gencodeTranscript,gencodeGene)
gencode <- gencode[order(gencode$V1,gencode$V4,gencode$V5),]

write.table(x = gencode,file = paste0(opt$ref,opt$gencodeName,"_MLT.gtf"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)

##
##Write gene files
writeGeneFiles(gencodeData = gencode, transcriptType = "MLT")
writeGeneFiles(gencodeData = gencodeGeneLT, transcriptType = "LT")
