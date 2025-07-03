#/!bin/bash
#This script conducts TFBS enrichment analysis

if [ -z $1 ]
then
echo "caABC file name (without .txt ending)"
echo "Exit now"
exit 1
elif [ -z $2 ]
then
echo "Name of group of genes to check for TFBS enrichment."
echo "Exit now"
exit 1
else
echo "All inputs given"
fi

mkdir -p ../../data/memeEval/ame/
mkdir -p ../../data/memeEval/ame/${1}_${2}/
mkdir -p ../../data/memeEval/tfbs/
wget -O ../../data/memeEval/tfbs/H12CORE_meme_format.meme https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12CORE/formatted_motifs/H12CORE_meme_format.meme -nc

##
##Get and index reference
wget -O ../../data/ref/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -nc
gunzip -c ../../data/ref/hg38.fa.gz > ../../data/ref/hg38.fa
samtools faidx ../../data/ref/hg38.fa

##
##Extract caABC ame
Rscript --vanilla extractGenesCaABCAme.R --caABCName ${1} --geneList ${2}
##
##Make fasta of regions
for i in Case Control
do
samtools faidx ../../data/ref/hg38.fa -r ../../data/memeEval/ame/${1}_${2}/regRegions${i}.txt -o ../../data/memeEval/ame/${1}_${2}/hg38RegRegions${i}.fa 
##
##Mask repeats with Dust
dust ../../data/memeEval/ame/${1}_${2}/hg38RegRegions${i}.fa > ../../data/memeEval/ame/${1}_${2}/hg38RegRegions${i}Dust.fa
done

##
##Apply AME
ame --seed 1234 --evalue-report-threshold 10 -oc ../../data/memeEval/ame/${1}_${2}/ameOutHocomoco --control ../../data/memeEval/ame/${1}_${2}/hg38RegRegionsControlDust.fa ../../data/memeEval/ame/${1}_${2}/hg38RegRegionsCaseDust.fa ../../data/memeEval/tfbs/H12CORE_meme_format.meme
