#!/bin/bash
## Download annotations for SnpEff (4.3)
conda_home=/srv/conda
# /databases/v4_3/snpEff_v4_3_GRCh38.86.zip
#sudo mkdir -p ${conda_home}/share/snpeff-4.3.1t-1/data/GRCh37.75
sudo mkdir -p ${conda_home}/share/snpeff-4.3.1t-1/data/GRCh38.86
sudo wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh38.86.zip -P ${conda_home}/share/snpeff-4.3.1t-1/
sudo rm -rf ${conda_home}/share/snpeff-4.3.1t-1/data/
cd ${conda_home}/share/snpeff-4.3.1t-1/
sudo unzip snpEff_v4_3_GRCh38.86.zip
${conda_home}/jre/bin/java -jar ${conda_home}/share/snpeff-4.3.1t-1/snpEff.jar download GRCh38.86