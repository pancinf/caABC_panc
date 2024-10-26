#!/bin/bash

#This script generates ABC(-like) models from ENCODE DNAse, H3K27AC and Hi-C data
###
###
sampleName=""
h3k27acPath=""
acCol=""
gencodeVers=""

while getopts "p:h:n:g:" o;
do
	case $o in
		p )	sampleName=${OPTARG};;
		h )	h3k27acPath=${OPTARG};;
		n )	acCol=${OPTARG};;
		g )	gencodeVers=${OPTARG};;
		? )
			echo "Given optin is not valid: -${OPTARG}"
			exit 1;;
esac
done

if [ -z ${sampleName} ]
then
echo "No sampleName given. Please provide with the -p argument"
exit 1
elif [ -z ${h3k27acPath} ]
then
echo "No Path(s) to H3K27AC-file(s) given. Please provide with the -h argument"
exit 1
elif [ -z ${acCol} ]
then
echo "No activity column number without qNorm given. Please provide with the -n argument (34 for samples with 1 H3K27AC file and 40 for samples with 2)"
exit 1
elif [ -z ${gencodeVers} ]
then
echo "No gencode version given. Please provide with the -g argument"
exit 1
fi

##
##
echo "Processed sample: ${sampleName}"

###
###
###MACS2

cd ../../data
mkdir -p macsOut/
mkdir -p macsOut/${sampleName}/

macs2 callpeak \
-t raw/${sampleName}Dnase.bam \
-n ${sampleName}.macs2 \
-f BAM \
-g hs \
-p 0.1 \
--call-summits \
--outdir macsOut/${sampleName} \
--shift -75 --extsize 150 --nomodel

#Sort narrowPeak file
bedtools sort -faidx ref/sizes.genome -i macsOut/${sampleName}/${sampleName}.macs2_peaks.narrowPeak > macsOut/${sampleName}/${sampleName}.macs2_peaks.narrowPeak.sorted

echo "---"
echo "---"
echo "---"
echo "Peaks called"
echo "---"
echo "---"
echo "---"

###
###
###ABC-Regions
mkdir -p abcReg/

for transcriptType in CT LT
do
	mkdir -p abcReg/${sampleName}_${transcriptType}

	python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/makeCandidateRegions.py \
	--narrowPeak macsOut/${sampleName}/${sampleName}.macs2_peaks.narrowPeak.sorted \
	--accessibility raw/${sampleName}Dnase.bam \
	--outDir abcReg/${sampleName}_${transcriptType} \
	--chrom_sizes ref/sizes.genome \
	--chrom_sizes_bed ref/sizes.genome.bed \
	--regions_blocklist ../software/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_unified_blacklist.bed \
	--peakExtendFromSummit 250 \
	--nStrongestPeaks 150000 \
	--regions_includelist ref/gencode.v${gencodeVers}.annotation_${transcriptType}_tss.bed
done

echo "---"
echo "---"
echo "---"
echo "Regions extracted"
echo "---"
echo "---"
echo "---"

###
###
###ABC-Activities
mkdir -p abcAct/
for transcriptType in CT LT
do
	mkdir -p abcAct/${sampleName}_${transcriptType}

	python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py \
	--candidate_enhancer_regions abcReg/${sampleName}_${transcriptType}/${sampleName}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
	--genes ref/gencode.v${gencodeVers}.annotation_${transcriptType}_gene.bed \
	--DHS raw/${sampleName}Dnase.bam \
	--H3K27ac ${h3k27acPath} \
	--default_accessibility_feature DHS \
	--chrom_sizes ref/sizes.genome \
	--chrom_sizes_bed ref/sizes.genome.bed  \
	--cellType Pancreas \
	--outdir abcAct/${sampleName}_${transcriptType}
	Rscript --vanilla ../scripts/abc/extract_activities.R --act abcAct/${sampleName}_${transcriptType}/EnhancerList.txt --actOut abcAct/${sampleName}_${transcriptType}/EnhancerList_act.bed \
	--actCol ${acCol}
done

echo "---"
echo "---"
echo "---"
echo "Activities calculated"
echo "---"
echo "---"
echo "---"

###
###
###Extract HI-C data
mkdir -p abcHic/
mkdir -p abcHic/${sampleName}

python ../software/ABC-Enhancer-Gene-Prediction/workflow/scripts/juicebox_dump_SCALE.py \
--hic_file raw/${sampleName}Hic.hic \
--juicebox "java -jar ../software/juicer/juicer_tools_2.13.06.jar" \
--resolution 5000 \
--outdir abcHic/${sampleName}

echo "---"
echo "---"
echo "---"
echo "HI-C extracted"
echo "---"
echo "---"
echo "---"

###
###
###Final models
mkdir -p finalModels/

##
##ABC
STARE_ABCpp -b abcAct/${sampleName}_LT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation_LT.gtf -o finalModels/ABC_${sampleName}_LT -f abcHic/${sampleName} -k 5000 -t 0 -c 1 -q False -i 5_tss

##
##gABC
STARE_ABCpp -b abcAct/${sampleName}_LT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation.gtf -o finalModels/gABC_${sampleName}_LT -f abcHic/${sampleName} -k 5000 -t 0 -c 1 -q True -i all_tss

##
##caABC
STARE_ABCpp -b abcAct/${sampleName}_CT/EnhancerList_act.bed -n 4 -a ref/gencode.v${gencodeVers}.annotation_CT.gtf -o finalModels/caABC_${sampleName}_CT -f abcHic/${sampleName} -k 5000 -t 0 -c 1 -q True -i 5_tss

echo "---"
echo "---"
echo "---"
echo "Final models generated"
echo "---"
echo "---"
echo "---"
