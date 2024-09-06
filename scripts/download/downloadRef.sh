#!/bin/bash
#This script downloads and prepares reference files

###
###
###Download gencode V44 GTF file
wget -O ../../data/ref/gencode.v44.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -c ../../data/ref/gencode.v44.annotation.gtf.gz > ../../data/ref/gencode.v44.annotation.gtf

###
###
###Integrate MANE information

###
###
###Get bigBedToBed
wget -O ../../software/ucscTools/bigBedToBed https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x ../../software/ucscTools/bigBedToBed

###
###
###Process MANE
wget -O ../../data/ref/mane.bb https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb
../../software/ucscTools/bigBedToBed ../../data/ref/mane.bb ../../data/ref/mane.bed

###
###
###Subset gencode GTF to longest transcript + MANE version
Rscript --vanilla processGtf.R --ref ../../data/ref/ --gencodeName gencode.v44.annotation

