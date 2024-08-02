#!/bin/bash

#This script generates ABC(-like) models from ENCODE DNAse, H3K27AC and Hi-C data
###
###
probeName=""
h3k27acCol=""
acCol=""
gencodeVers=""

while getopts "p:h:n:g:" o;
do
	case $o in
		p )	probeName=${OPTARG};;
		h )	h3k27acCol=${OPTARG};;
		n )	acCol=${OPTARG};;
		g )	gencodeVers=${OPTARG};;
		? )
			echo "Given optin is not valid: -${OPTARG}"
			exit 1;;
esac
done

if [ -z ${probeName} ]
then
echo "No probeName given. Please provide with the -p argument"
exit 1
elif [ -z ${h3k27acCol} ]
then
echo "No Path(s) to H3K27AC-file(s) given. Please provide with the -h argument"
exit 1
elif [ -z ${acCol} ]
then
echo "No activity column number without qNorm given. Please provide with the -n argument (35 for probes with 1 H3K27AC file and 41 for probes with 2)"
exit 1
elif [ -z ${gencodeVers} ]
then
echo "No gencode version given. Please provide with the -g argument"
exit 1
fi

##
##
echo "Processed probe: ${probeName}"

###
###
###MACS2

cd ../../data
mkdir -p macsOut/
mkdir -p macsOut/${probeName}/

macs2 callpeak \
-t raw/${probeName}Dnase.bam \
-n ${probeName}.macs2 \
-f BAM \
-g hs \
-p 0.1 \
--call-summits \
--outdir macsOut/${probeName} \
--shift -75 --extsize 150 --nomodel

#Sort narrowPeak file
bedtools sort -faidx ref/sizes.genome -i macsOut/${probeName}/${probeName}.macs2_peaks.narrowPeak > macsOut/${probeName}/${probeName}.macs2_peaks.narrowPeak.sorted

###
###
###ABC-Regions
mkdir -p abcReg/

for transcriptType in MLT LT
do
	mkdir -p abcReg/${probeName}_${transcriptType}

	python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/makeCandidateRegions.py \
	--narrowPeak macsOut/${probeName}/${probeName}.macs2_peaks.narrowPeak.sorted \
	--accessibility raw/${probeName}Dnase.bam \
	--outDir abcReg/${probeName}_${transcriptType} \
	--chrom_sizes ref/sizes.genome \
	--chrom_sizes_bed ref/sizes.genome.bed \
	--regions_blocklist ../software/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_unified_blacklist.bed \
	--peakExtendFromSummit 250 \
	--nStrongestPeaks 150000 \
	--regions_includelist ref/gencode.v${gencodeVers}.annotation_${transcriptType}_tss.bed
done

###
###
###ABC-Activities
mkdir -p abcAct/
for transcriptType in MLT LT
do
	mkdir -p abcAct/${probeName}_${transcriptType}

	python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py \
	--candidate_enhancer_regions abcReg/${probeName}_${transcriptType}/${probeName}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
	--genes ref/gencode.v${gencodeVers}.annotation_${transcriptType}_gene.bed \
	--DHS raw/${probeName}Dnase.bam \
	--H3K27ac ${h3k27acCol} \
	--default_accessibility_feature DHS \
	--chrom_sizes ref/sizes.genome \
	--chrom_sizes_bed ref/sizes.genome.bed  \
	--cellType Pancreas \
	--outdir abcAct/${probeName}_${transcriptType}
	Rscript --vanilla ../scripts/abc/extract_activities.R --act abcAct/${probeName}_${transcriptType}/EnhancerList.txt --actOut abcAct/${probeName}_${transcriptType}/EnhancerList_act.bed \
	--actCol ${acCol}
done

###
###
###Extract Hi-C data
mkdir -p abcHic/
mkdir -p abcHic/${probeName}

python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/juicebox_dump_SCALE.py \
--hic_file raw/${probeName}Hic.hic \
--juicebox "java -jar ../software/juicer/juicer_tools_2.13.06.jar" \
--resolution 5000 \
--outdir abcHic/${probeName} 

###
###
###Final models
mkdir -p finalModels/

##
##ABC
STARE_ABCpp -b abcAct/${probeName}_LT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation_LT.gtf -o finalModels/ABC_${probeName}_LT -f abcHic/${probeName} -k 5000 -t 0 -c 1 -q False -i 5_tss

##
##gABC
STARE_ABCpp -b abcAct/${probeName}_LT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation.gtf -o finalModels/gABC_${probeName}_LT -f abcHic/${probeName} -k 5000 -t 0 -c 1 -q True -i all_tss

##
##M-AD-ABC
STARE_ABCpp -b abcAct/${probeName}_MLT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation_MLT.gtf -o finalModels/maABC_${probeName}_MLT -f abcHic/${probeName} -k 5000 -t 0 -c 1 -q True -i 5_tss
