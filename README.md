# caABC - an Improved and Interpretable Activity-By-Contact based Model for Predicting Enhancer-Promoter Interactions

## Background

The Activity-By-Contact model (Fulco et al. 2019) is one of the most promising approaches for defining tissue-specific enhancer-promoter (EP) interactions.
This model combines enhancer signal (DNase- or ATAC-seq + H3K27ac) with Hi-C data. The higher the signal and/or the higher the contact frequency, the higher an enhancer-promoter prediction will be scored (normalized for other potential enhancers of the respective target gene).
However, this model definition leads to an "over-mapping" of highly active enhancers to many nearby genes.
The recently published gABC model (Hecker et al. 2023) targets this issue by adding a further adaption step (normalizing for Hi-C-interaction signals to all nearby loci).<br />
Also Hecker et al. compared model performance between using longest-transcript based promoters and averaging accross predictions considering all promoters for each gene.
Considering all promoters gave an improvement and was integrated as part of gABC.
Thus with the final prediction table of gABC for each regulatory region only the mean distance to all considered promoters and a merged ABC score is given.

Consequences:
* One can not extract each respective promoters influence on the gene-wise prediction score. 
* One can not seperate between enhancer and promoter by distance at all and has to consider both for analyses (as it has been done by Hecker et al.)
* In practical terms for e.g. GWAS studies when using this model one has to work with all cosidered transcripts for the respective promoters instead of a single representative transcript, thus making analysis complex.

We developed a new model (caABC), which includes the gABC adaption step. For improved interpretation caABC does not consider all promoters but ENSEMBL canonical-transcript based ones (https://www.ensembl.org/info/genome/genebuild/canonical.html).

In our model tissue, the pancreas, caABC shows **significant outperformance** compared to a longest-transcript based ABC model while being **more interpretable** than gABC.

## Workflow to generate our results in pancreas tissue

### Install and download
* Install environments
	* Install miniconda https://docs.anaconda.com/miniconda/miniconda-install/
	* Install environments from our provided **abc.yml** and **evalAbc.yml** by the command **conda env create -f environment.yml**
* Clone ABC
	* Clone ABC from https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction?tab=readme-ov-file to the software folder
	* We used version 25e7de9007a5880efc8e5340adb40a3d324f81c4
* Raw data
	* Required conda environment: abc
	* Go to the scripts/download folder and run ./downloadAll.sh
	* This will also download juicer_tools required for Hi-C extraction
	* Take some tea - this will take time :)

### Generate ABC-based models
* Description: Generates ABC, gABC and caABC predictions from raw epigenomic data.
* Required conda environment: abc
* script: abc.sh
* Input arguments
	* -p: name of sample (as provided when downloading with the download script)
	* -n: column number of activity column. Choose 34 for samples with one H3K27AC file and 40 for samples with two respective files.
	* -g: gencode version (if you want to reproduce our results use 44 as it is done in the download script).
	* -h: path to H3K27AC files. Structure: "raw/" + filename. Comma-separated for multiple files.
* Example command (for a sample with 2 H3K27AC-files)
	* ./abc.sh -p sample1 -n 40 -g 44 -h raw/sample1AH3K27ac.bam,raw/sample1BH3K27ac.bam

### CAGE-seq evaluation
* Description: Compare CAGE-Seq signal strength between canonical-transcript and longest-transcript based promoters.
* Required conda environment: evalAbc
* script: cageEval.R
* Input arguments
	* --canonicalTrans: Name of canonical-transcript based transcript file (should be in your "data/ref" folder after downloading).
	* --longestTrans: Name of longest-transcript based transcript file (should be in your "data/ref" folder after downloading).
	* --cagePlus: Name of plus strand CAGE-seq file (should be in your "data/cage/signals" folder after downloading).
	* --cageMinus: Name of plus strand CAGE-seq file (should be in your "data/cage/signals" folder after downloading).
	* --typeAnno: Path to "*_gene.bed" file in the "data/ref" folder. Required to assign gene-types.
	* --keepTypes: Gene-types to keep (Default = "protein_coding,lncRNA").
	* --minDist: Minimum distance (bp) between matched longest-transcript and canonical-transcript based promoters to consider the gene (Default = 5000).
	* --maxDist: Maximum distance (bp) of a CAGE-seq signal from the TSS to consider it as a promoter signal (Default = 250).
	* --cap: capping value for plotting e.g. 0.05 means that the values above the 95% quantile and below the 5% quantile are set to 95% and 5% for plotting respectively.
	* --tissue: Tissue name of interest (matched to GTEx naming).
	* --gtexCut: Minimum GTEx V8 median-TPM value to consider a gene expressed (Default = 1).
* **Command to reproduce our results**: Rscript --vanilla cageEval.R --canonicalTrans gencode.v44.annotation_CT_gene.bed --longestTrans gencode.v44.annotation_LT_gene.bed --cagePlus pancreas-Adult.plus.bedGraph --cageMinus pancreas-Adult.minus.bedGraph --typeAnno ../../data/ref/gencode.v44.annotation_CT_gene.bed --keepTypes "protein_coding,lncRNA" --minDist 5000 --maxDist 250 --cap 0.05 --tissue Pancreas --gtexCut 1

### Compare ABC, gABC and caABC models
* Description: Pairwise compare ABC, gABC and caABC models capability to predict fine-mapped eQTLs.
* Required conda environment: evalAbc
* script: evalEQTL.R
* Input arguments
	* --inputList: Name of file containing comparison information (should be located in "data/GTEx/evalGTEx/inputLists/"). The structure of the inputList can be seen in the data/ListExamples folder (inputListCompare.txt)
	* --typeAnno: Path to "*_gene.bed" file in the "data/ref" folder. Required to assign gene-types.
	* --keepTypes: Gene-types to keep (Default = "protein_coding,lncRNA").
	* --minDist: minimum distance (bp) to TSS to consider a regulatory region an enhancer (Default = 250 / note that this argument in practice always works as "larger than", so in the default mode it considers all elements as enhancers which are at least 251 bp away from the TSS).
	* --gtexCut: Minimum GTEx V8 median-TPM value to consider a gene expressed (Default = 1).
	* --caviarThres: Minimum fine-mapping probability to consider an eQTL causal (Default = 0.8).
	* --tissue: Tissue name of interest (matched to GTEx naming).
 	* --cap: capping value for plotting e.g. 0.05 means that the values above the 95% quantile and below the 5% quantile are set to 95% and 5% for plotting respectively.
 	* --compMod: Tag to indicate which models for which sample are compared.
 	* --outDir: Name of output-directory.
  * **Command to reproduce our results**: Rscript --vanilla evalEQTL.R --inputList abc_gabc_Sample1.txt --typeAnno ../../data/ref/gencode.v44.annotation_CT_gene.bed --keepTypes "protein_coding,lncRNA" --minDist -1 --gtexCut 1 --caviarThres 0.8 --tissue Pancreas --cap 0.10 --compMod abc_gabc_Sample1 --outDir sample1

### Generate enhancer-promoter predictions per sample
* Description: Generate files for genes, promoters and respective enhancers for the tissue of interest. Files can be easily loaded in IGV.
* Required conda environment: evalAbc
* Input arguments
	* --caABCPath: Path to caABC prediction file.
	* --sample: Sample name 
	* --typeAnno: Path to "*_gene.bed" file in the "data/ref" folder. Required to assign gene-types.
	* --keepTypes: Gene-types to keep (Default = "protein_coding,lncRNA").
	* --minDist: minimum distance (bp) to TSS to consider a regulatory region an enhancer (Default = 250).
	* --gtexCut: Minimum GTEx V8 median-TPM value to consider a gene expressed.
	* --tissue: Tissue name of interest (matched to GTEx naming).
 	* --caviarThres: Minimum fine-mapping probability to consider an eQTL causal.
 	* --recCut: eQTL recall cut-off.
 	* --cGtf: Path to CT-GTF file.
 	* --outDir: Name of output directory.
* **Command to reproduce our results**: Rscript --vanilla recSub.R --caABCpath ../../data/finalModels/caABC_sample1_CT_ABCpp_scoredInteractions_c4.txt.gz --sample sample1 --typeAnno ../../data/ref/gencode.v44.annotation_CT_gene.bed --keepTypes "protein_coding,lncRNA" --minDist 250 --gtexCut 1 --tissue Pancreas --caviarThres 0.8 --recCut 0.7 --cGtf ../../data/ref/gencode.v44.annotation_CT.gtf --outDir sample1

### Generate replication-based enhancer-promoter predictions
* Description: Make a unified enhancer, promoter and gene list. Enhancers are only kept if found in at least 2 samples. Files can be easily loaded in IGV.
* Required conda environment: evalAbc
* Input Arguments
	* --caABCRepliList: Name of file containing promoter, enhancer and gene paths to use for replication filtering located at "data/GTEx/caAbcRepli/caAbcRepliInList/". The structure of the inputList can be seen in the data/ListExamples folder (inputListRepli.txt)
	* --minOvBpEnh: Minimum overlap (bp) to consider an enhancer replicated (Default = 250).
	* --geneAnnoFile: Path to "*_gene.bed" file in the "data/ref" folder. Required to assign gene-types.
	* --outTag: Tag set in the beginning of output file names.
* **Command to reproduce our results**: Rscript --vanilla caAbcRepli.R --caABCRepliList pancreas_caAbcRepliInList.txt --geneAnnoFile gencode.v44.annotation_CT_gene.bed --minOvBpEnh 250 --outTag caABC_pancreas

## Novel predictions in other tissues
* With small adaptions one can use this code to generate novel caABC models for other tissues based on publicaly available ENCODE data. For comparison of models and final EP-List generation (at a certain cut-off) we use fine-mapped tissue-specific eQTLs from the GTEx V8 CAVIAR dataset. You cannot use our framework if your tissue of interest is not represented in the GTEx CAVIAR eQTL dataset.
* You have to
	* change the files to download in the download script(s).
		* DNase-seq - sample(s) of interest
		* H3K27AC - sample(s) of interest
		* HI-C - sample(s) of interest
		* CAGE-seq - UCSC tissue of interest.
	* keep in mind that there might be samples with more than 1 H3K27AC file from ENCODE and consider this in the raw-download script (downloadRaw.sh).
	* assign the correct tissue name from GTEx for correct eQTL and expression extraction.
* You can speed up and reduce complexity by ignoring the CAGE-seq evaluation and model comparison if you just want to get the final caABC predictions.
