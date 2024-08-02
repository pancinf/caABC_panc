#!/bin/bash
#This script downloads and preprocesses pancreatic ENCODE data for ABC models..

##
###Paths + names
ProbesOldDnase="ENCFF942NPB ENCFF135KGR ENCFF979ZTM ENCFF456GDF"
ProbesOldHic="ENCFF586MQY ENCFF140ROU ENCFF653HCO ENCFF705MKK"
ProbesOldH3k27ac="ENCFF077ETR ENCFF217SAY ENCFF721BMD ENCFF589NAM ENCFF939BNJ"
ProbesNewDnaseHic="probe1 probe2 probe3 probe4"
ProbesNewH3k27ac="probe1A probe1B probe2 probe3 probe4"
arrProbesOldDnase=($ProbesOldDnase)
arrProbesOldHic=($ProbesOldHic)
arrProbesOldH3k27ac=($ProbesOldH3k27ac)
arrProbesNewDnaseHic=($ProbesNewDnaseHic)
arrProbesNewH3k27ac=($ProbesNewH3k27ac)


##
##Dnase + Hi-C
for probe in 0 1 2 3
do
        wget -O ../../data/raw/"${arrProbesNewDnaseHic[$probe]}"Dnase.bam https://www.encodeproject.org/files/"${arrProbesOldDnase[$probe]}"/@@download/"${arrProbesOldDnase[$probe]}".bam
        samtools index ../../data/raw/"${arrProbesNewDnaseHic[$probe]}"Dnase.bam

	wget -O ../../data/raw/"${arrProbesNewDnaseHic[$probe]}"Hic.hic https://www.encodeproject.org/files/"${arrProbesOldHic[$probe]}"/@@download/"${arrProbesOldHic[$probe]}".hic
done

##
##H3k27ac
for probe in 0 1 2 3 4
do
	wget -O ../../data/raw/"${arrProbesNewH3k27ac[$probe]}"H3k27ac.bam https://www.encodeproject.org/files/"${arrProbesOldH3k27ac[$probe]}"/@@download/"${arrProbesOldH3k27ac[$probe]}".bam
	samtools index ../../data/raw/"${arrProbesNewH3k27ac[$probe]}"H3k27ac.bam
done
