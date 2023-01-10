#!/bin/bash

## A script to rerun the pb-cpg-tools command on every line in a .bed file and provide methylation visualization.

# BED file for genomic segments of interest.
		# Make argument to enter file name in bash command
HK=../genes.bed

# Flanking region around the segments of interest.
FLANK=10000

# Path where the directory for pb-cpg-tools documents are kept.
PT="/Users/shalvichirmade/Documents/pb-cpg-tools & Visualizations/pb-CpG-tools-1.1.0/"

# Path where the hg38.fa file is kept. Please make sure it is unzipped.
HG38="/Users/shalvichirmade/Documents/pb-cpg-tools & Visualizations/Script/hg38.fa"

# Create a folder in your current directory where all the extra files are kept.
mkdir samtools_pb_cpg_tools_files
EXFILE="./samtools_pb_cpg_tools_files/"


# Conduct analysis.
while read line;
do
	#echo $line
	CHR=$(echo $line | cut -f1 -d " ")
	START=$(echo $line | cut -f2 -d " ")
	END=$(echo $line | cut -f3 -d " ")

	START_FLANK="$((${START}-${FLANK}))"
	END_FLANK="$((${END}+${FLANK}))"

	# Make sure start value is not negative.
	if [ $START_FLANK -lt 1 ]
	then
		START_FLANK=1
	else
		:
	fi

	SEG=${CHR}:${START_FLANK}-${END_FLANK}
	FILE=${CHR}_${START_FLANK}_${END_FLANK}
	# echo ${CHR} and ${START} ${END}
	# echo ${SEG} and ${FILE}

	# Check to see if the output file already exists for each segment. If yes, skips the rest of the steps.
	if [ -f "${EXFILE}cpg_${FILE}.combined.reference.bed" ]; 
	then
    	echo "${EXFILE}cpg_${FILE}.combined.reference.bed already exists and process will be skipped."
	else

		# Download methylation data for each genome segment of interest.
		samtools view -b -h https://downloads.pacbcloud.com/public/dataset/HG002-CpG-methylation-202202/HG002.GRCh38.haplotagged.bam ${SEG} > ${EXFILE}${FILE}.bam

		echo "File created for" ${FILE}

		# Create an indexed file for segment.
		samtools index -b ${EXFILE}${FILE}.bam ${EXFILE}${FILE}.bam.bai

		echo "Index created for" ${FILE}

		# Using pb-cpg-tools to conduct an analysis of the segment: CpG/5mC data from PacBio HiFi reads aligned to a reference genome.
		python "${PT}aligned_bam_to_cpg_scores.py" -b ${EXFILE}${FILE}.bam -f "${HG38}" -o ${EXFILE}cpg_${FILE} -p model -m reference -d "${PT}pileup_calling_model/"

		echo "Analysis conducted for" ${FILE} 
		echo ""


	fi

	Rscript bashscript_GC_Methylation_Visuzalization.R ${EXFILE} cpg_${FILE}.combined.reference.bed ${CHR} 

done <$HK

# Combine all .combined.reference.bed files together. The resulting file will be in your current directory while all the extra files are in samtools_pb_cpg_tools_files.
# cat ${EXFILE}*.combined.reference.bed > all_combined_reference.bed

# echo "Combined all .combined.reference.bed files together."
