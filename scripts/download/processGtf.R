#This scripts generates input Transcripts, TSS and a gtf-file for ABC, gABC and caABC
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
writeGeneFiles <- function(gencodeData, feature, transcriptType){
	abcGenes <- gencodeData[,c(1,4,5,9,7)]
	if(feature == "transcript"){
		abcGenes <- colsplit(gencodeData$V9, "; ", c("gene_id","transcript_id", "gene_type", "gene_name","remaining"))
		abcGenes <- data.frame(chr = gencodeData$V1,start = gencodeData$V4, end = gencodeData$V5, name = abcGenes[,4],
                                                score = 0, strand = gencodeData$V7,Ensembl_ID = abcGenes[,1], gene_type = abcGenes[,3])
	}else if (feature == "gene"){
		abcGenes <- colsplit(gencodeData$V9, "; ", c("gene_id", "gene_type", "gene_name","remaining"))
		abcGenes <- data.frame(chr = gencodeData$V1,start = gencodeData$V4, end = gencodeData$V5, name = abcGenes[,3],
                                                score = 0, strand = gencodeData$V7,Ensembl_ID = abcGenes[,1], gene_type = abcGenes[,2])
	}
	abcGenes$remaining <- NULL
	
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
##Read gtf

gencode <- fread(paste0(opt$ref,opt$gencodeName,".gtf"),data.table = FALSE, header = FALSE)
colnames(gencode)[1] <- "V1"

##
##Write gtf files

#Write longest transcript file
write.table(x = gencode[gencode$V3 == "gene",],file = paste0(opt$ref,opt$gencodeName,"_LT.gtf"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)

#Write canonical transcript file
write.table(x = gencode[gencode$V3 == "transcript" & grepl(pattern = "Ensembl_canonical", x = gencode$V9),],file = paste0(opt$ref,opt$gencodeName,"_CT.gtf"),sep = "\t",append = FALSE,quote = FALSE,row.names = FALSE,col.names = FALSE)

##
##Write TSS and gene files

#make 0-based
gencode$V4 <- gencode$V4 -1

#write
writeGeneFiles(gencodeData = gencode[gencode$V3 == "transcript" & grepl(pattern = "Ensembl_canonical", x = gencode$V9),], feature = "transcript", transcriptType = "CT")
writeGeneFiles(gencodeData = gencode[gencode$V3 == "gene",], feature = "gene", transcriptType = "LT")

