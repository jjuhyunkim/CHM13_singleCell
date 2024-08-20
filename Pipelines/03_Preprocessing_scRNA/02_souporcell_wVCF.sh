#! /bin/bash

cellRangerOut=$1 # the directory of Cellragner output
outDir=$2 # output directory
k=$3 # number of multiplexed samples.
vcf=$4 # VCF files including multiplexed samples
knownSample=$5 # name of sample in VCF
fasta=$6 # single cell reference
souporcellSif=${dir_to_souporcell_env}/souporcell2.5.sif

bam="possorted_genome_bam.bam"
barcode=`basename $(ls ${cellRangerOut}/*_cell_barcodes.csv)`

thread=15

echo -e "${dir_to_souporcell_env}/souporcell_pipeline.py \
        -i ${cellRangerOut}/${bam} \
        -b ${cellRangerOut}/${barcode} \
	--known_genotypes ${vcf} \
	--known_genotypes_sample_names ${knownSample} \
        -f ${fasta} -t ${thread} -o ${outDir} -k ${k} ${RESET}"

${dir_to_souporcell_env}/souporcell_pipeline.py \
	-i ${cellRangerOut}/${bam} \
	-b ${cellRangerOut}/${barcode} \
	--known_genotypes ${vcf} \
        --known_genotypes_sample_names ${knownSample} \
	-f ${fasta} -t ${thread} -o ${outDir} -k ${k}
