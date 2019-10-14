#Index the genome fasta file, using samtools and bwa, only if does not exists a file with a md5 hash identical to a hash annotated in a file generated after
#a successful bwa run below
if [[ -n $1 ]];
then
	qiagen_parent_folder=$1
	conda_home=$2
	conda_preffix=$3
	conda_env=$4
else
	qiagen_parent_folder=/srv/qgen
	conda_home=/srv/conda
	conda_preffix=/srv/conda
	conda_env=base
fi

#md5sum -c ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
ls ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
## Download genome files
    wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P ${qiagen_parent_folder}/data/genome/
    cd ${qiagen_parent_folder}/data/genome && \
        gunzip ucsc.hg19.fa.gz  && \
        ${conda_home}/bin/samtools faidx ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa && \
        ${conda_home}/bin/bwa index ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa
        md5sum -b ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac > ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 )
cd ${qiagen_parent_folder}/qiaseq-dna
