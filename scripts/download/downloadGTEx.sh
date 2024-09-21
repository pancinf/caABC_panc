#!/bin/bash
#This script downloads GTEx files (including CAVIAR) for V8.

###
###
###Download TPM data
wget -O ../../data/GTEx/raw/GTExTPM.gct.gz https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
gunzip -c ../../data/GTEx/raw/GTExTPM.gct.gz > ../../data/GTEx/raw/GTExTPM.gct

echo "---"
echo "---"
echo "---"
echo "TPMs downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Download GTEx GTF file
wget -O ../../data/GTEx/raw/gencodeGTEx.gtf https://storage.googleapis.com/adult-gtex/references/v8/reference-tables/gencode.v26.GRCh38.genes.gtf

echo "---"
echo "---"
echo "---"
echo "GTEx GTF downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Download GTEx CAVIAR file
wget -O ../../data/GTEx/raw/GTEx_v8_finemapping_CAVIAR.tar https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/fine-mapping-cis-eqtl/GTEx_v8_finemapping_CAVIAR.tar
tar -xvf ../../data/GTEx/raw/GTEx_v8_finemapping_CAVIAR.tar -C ../../data/GTEx/raw/

echo "---"
echo "---"
echo "---"
echo "GTEx CAVIAR downloaded"
echo "---"
echo "---"
echo "---"