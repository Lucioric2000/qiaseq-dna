#!/bin/bash
#You should run this example using the command:
#source NEB_S2_and_cellline_example.txt
#Go to the project dir
cd /srv/qgen/qiaseq-dna
#Activate the conda environment
source /root/conda/bin/activate base
#Make hard links from the sample files NEB_S2_L001_R?_001.fastq.gz to files named Cellline_S10_R?_001.fastq.gz, where ? stands for 1 and for 2
ln /srv/qgen/example/NEB_S2_L001_R1_001.fastq.gz /srv/qgen/example/Cellline_S10_L001_R1_001.fastq.gz
ln /srv/qgen/example/NEB_S2_L001_R2_001.fastq.gz /srv/qgen/example/Cellline_S10_L001_R2_001.fastq.gz
time python run_qiaseq_dna.py forcellline.txt v2 single out2v6_{0} NEB_S2 Cellline_S10 &> run_smcounterv2_cellline.log &