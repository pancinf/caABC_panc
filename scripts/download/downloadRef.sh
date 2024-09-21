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
###Integrate MANE information

###
###
###Get bigBedToBed
wget -O ../../software/ucscTools/bigBedToBed https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x ../../software/ucscTools/bigBedToBed

echo "---"
echo "---"
echo "---"
echo "bigBedToBed downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Process MANE
wget -O ../../data/ref/mane.bb https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb
../../software/ucscTools/bigBedToBed ../../data/ref/mane.bb ../../data/ref/mane.bed

echo "---"
echo "---"
echo "---"
echo "MANE downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Subset gencode GTF to longest transcript + MANE
Rscript --vanilla processGtf.R --ref ../../data/ref/ --gencodeName gencode.v44.annotation

echo "---"
echo "---"
echo "---"
echo "Longest-Transcript and MANE based GTF created"
echo "---"
echo "---"
echo "---"

