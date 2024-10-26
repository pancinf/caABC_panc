#!/bin/bash
#This script downloads and prepares reference files

###
###
###Download gencode V44 GTF file
wget -O ../../data/ref/gencode.v44.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -c ../../data/ref/gencode.v44.annotation.gtf.gz > ../../data/ref/gencode.v44.annotation.gtf

echo "---"
echo "---"
echo "---"
echo "GENCODE GTF downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Subset gencode GTF to longest transcript + canonical transcript
Rscript --vanilla processGtf.R --ref ../../data/ref/ --gencodeName gencode.v44.annotation

echo "---"
echo "---"
echo "---"
echo "Longest-transcript and Canonical-transcript based GTF created"
echo "---"
echo "---"
echo "---"

