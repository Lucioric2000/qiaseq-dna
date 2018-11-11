#Index the genome fasta file, using samtools and bwa, only if does not exists a file with a md5 hash identical to a hash annotated in a file generated after
#a successful bwa run below
srv_qiagen=/srv/qgen
#md5sum -c ${srv_qiagen}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
ls ${srv_qiagen}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
## Download genome files
    wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P ${srv_qiagen}/data/genome/
    cd ${srv_qiagen}/data/genome && \
        gunzip ucsc.hg19.fa.gz  && \
        ${CONDA_PREFIX}/bin/samtools faidx ${srv_qiagen}/data/genome/ucsc.hg19.fa && \
        ${CONDA_PREFIX}/bin/bwa index ${srv_qiagen}/data/genome/ucsc.hg19.fa
        md5sum -b ${srv_qiagen}/data/genome/ucsc.hg19.fa.pac > ${srv_qiagen}/data/genome/ucsc.hg19.fa.pac.md5 )
cd ${srv_qiagen}/qiaseq-dna
