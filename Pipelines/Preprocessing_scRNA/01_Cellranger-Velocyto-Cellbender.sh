#!/bin/bash

# v1 20.07.01
# v2 20.08.03
# v3 21.01.06
# v4 22.11.28

## Arguments ##
lib=$1
sample=$2
pool=$3
wd=$4
expectedCell=$5

thread=10

## Software ##
cellranger=${path_to_cellragner}/cellragner
souporcell=${path_to_souporcell}/souporcell_pipeline.py
samtools=${path_to_samtools}/samtools

## DBs ##
RNA_ref=${path_to_10X_reference}

echo $RNA_ref
############################################### Running #######################################################################


source ~/anaconda3/etc/profile.d/conda.sh


# CellRanger Count #
echo "Start Cellranger"
cd ${wd}
${cellranger} count --id=${sample} --fastqs=${wd}/fastq/ --sample=${lib} --transcriptome=${RNA_ref} --localcores=8 --localmem=64

## suporcell
if [ ${pool} != 1 ] ; then
echo "souporcell start!"
conda activate GMIsouporcell
mkdir -p ${wd}/souporcell
gzip -dcf ${wd}/${lib}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ${wd}/${lib}/outs/filtered_feature_bc_matrix/barcodes.tsv
${souporcell} -i ${wd}/${lib}/outs/possorted_genome_bam.bam -b ${wd}/${lib}/outs/filtered_feature_bc_matrix/barcodes.tsv -f ${rna_ref}/fasta/genome.fa -t ${thread} -o ${wd}/souporcell -k ${pool}
conda deactivate
echo "souporcell done!"
else
echo "Souporcell not Run bc of Pooling X"
fi

## velocyto
echo "Start Velocyto"
conda activate scvelo
mkdir ${wd}/velocyto
cd ${wd}/velocyto
velocyto run10x -m ~/Files/repeatmasker/GRCh38_rmsk.gtf ${wd}/${sample} ${RNA_ref}/genes/genes.gtf
conda deactivate


# cell bender
echo "Start CellBender"
conda activate cellbender
mkdir ${wd}/cellbender
cellbender remove-background --input ${wd}/${sample}/outs/raw_feature_bc_matrix --output ${wd}/cellbender/${sample}.cellbender.h5 --expected-cells ${expectedCell} 
conda deactivate

## v1 20.07.01
## v2 20.08.03
## v3 21.01.06 
# GMI server route


