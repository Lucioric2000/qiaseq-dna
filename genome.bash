#!/bin/bash
#Index the genome fasta file, using samtools and bwa, only if does not exists a file with a md5 hash identical to a hash annotated in a file generated after
#a successful bwa run below
if [[ -n $1 ]];
then
	qiagen_parent_folder=$1
    conda_home=$2
	conda_env=$3
	GENOME_BUILD_ALT_NAME=$4
else
	qiagen_parent_folder=/srv/qgen
    conda_home=/root/conda
    conda_env=python37
	GENOME_BUILD_ALT_NAME=hg19
fi

mkdir -p ${qiagen_parent_folder}/data/genome/${GENOME_BUILD_ALT_NAME}
#md5sum -c ${qiagen_parent_folder}/data/genome/ucsc.${GENOME_BUILD_ALT_NAME}.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
ls ${qiagen_parent_folder}/data/genome/ucsc.${GENOME_BUILD_ALT_NAME}.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
## Download genome files
    cd ${qiagen_parent_folder}/data/genome/${GENOME_BUILD_ALT_NAME} && \
        wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.${GENOME_BUILD_ALT_NAME}.dict \
            https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.${GENOME_BUILD_ALT_NAME}.fa.gz && \
        gunzip ucsc.${GENOME_BUILD_ALT_NAME}.fa.gz && \
        samtools faidx ucsc.${GENOME_BUILD_ALT_NAME}.fa && \
        bwa index ucsc.${GENOME_BUILD_ALT_NAME}.fa && \
        md5sum -b ucsc.${GENOME_BUILD_ALT_NAME}.fa.pac > ucsc.${GENOME_BUILD_ALT_NAME}.fa.pac.md5 )
cd ${qiagen_parent_folder}/qiaseq-dna
