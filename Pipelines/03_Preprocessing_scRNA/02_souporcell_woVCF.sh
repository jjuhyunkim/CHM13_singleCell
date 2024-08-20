#! /bin/bash

cellRangerOut=$1
outDir=$2
k=$3

bam="possorted_genome_bam.bam"
barcode=`basename $(ls ${cellRangerOut}/*_cell_barcodes.csv)`

fasta=$6 # single cell reference
souporcellSif=${dir_to_souporcell_env}/souporcell2.5.sifthread=15

echo -e "${CYAN}CMD:singularity exec ${souporcellSif} \
        souporcell_pipeline.py \
        -i ${cellRangerOut}/${inputBam} \
        -b ${cellRangerOut}/${barcode} \
        -f ${fasta} -t ${thread} -o ${outDir} -k ${k} ${RESET}"

singularity exec --fakeroot --privileged ${souporcellSif} \
	souporcell_pipeline.py \
	-i ${cellRangerOut}/${bam} \
	-b ${cellRangerOut}/${barcode} \
	-f ${fasta} -t ${thread} -o ${outDir} -k ${k}
