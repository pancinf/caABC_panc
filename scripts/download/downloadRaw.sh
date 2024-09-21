#!/bin/bash
#This script downloads and preprocesses pancreatic ENCODE data for ABC models..

###
###
###Paths + names
SamplesOldDnase="ENCFF942NPB ENCFF135KGR ENCFF979ZTM ENCFF456GDF"
SamplesOldHic="ENCFF586MQY ENCFF140ROU ENCFF653HCO ENCFF705MKK"
SamplesOldH3K27ac="ENCFF077ETR ENCFF217SAY ENCFF721BMD ENCFF589NAM ENCFF939BNJ"
SamplesNewDnaseHic="sample1 sample2 sample3 sample4"
SamplesNewH3K27ac="sample1A sample1B sample2 sample3 sample4"
arrSamplesOldDnase=($SamplesOldDnase)
arrSamplesOldHic=($SamplesOldHic)
arrSamplesOldH3K27ac=($SamplesOldH3K27ac)
arrSamplesNewDnaseHic=($SamplesNewDnaseHic)
arrSamplesNewH3K27ac=($SamplesNewH3K27ac)

##
##Main

##
##Dnase + Hi-C
for sample in 0 1 2 3
do
        wget -O ../../data/raw/"${arrSamplesNewDnaseHic[$sample]}"Dnase.bam https://www.encodeproject.org/files/"${arrSamplesOldDnase[$sample]}"/@@download/"${arrSamplesOldDnase[$sample]}".bam
        samtools index ../../data/raw/"${arrSamplesNewDnaseHic[$sample]}"Dnase.bam

	wget -O ../../data/raw/"${arrSamplesNewDnaseHic[$sample]}"Hic.hic https://www.encodeproject.org/files/"${arrSamplesOldHic[$sample]}"/@@download/"${arrSamplesOldHic[$sample]}".hic
done

echo "---"
echo "---"
echo "---"
echo "DNase-Seq + HI-C downloaded"
echo "---"
echo "---"
echo "---"

##
##H3K27ac
for sample in 0 1 2 3 4
do
	wget -O ../../data/raw/"${arrSamplesNewH3K27ac[$sample]}"H3K27ac.bam https://www.encodeproject.org/files/"${arrSamplesOldH3K27ac[$sample]}"/@@download/"${arrSamplesOldH3K27ac[$sample]}".bam
	samtools index ../../data/raw/"${arrSamplesNewH3K27ac[$sample]}"H3K27ac.bam
done

echo "---"
echo "---"
echo "---"
echo "H3K27ac downloaded"
echo "---"
echo "---"
echo "---"