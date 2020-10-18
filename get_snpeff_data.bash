#!/bin/bash
## Download annotations for SnpEff
#conda_home=/root/conda
#snpeff_version=4.3.1t-3
#snpeff_version_short=4_3
#genome_version=GRCh37.75
conda_home=$1
snpeff_version=$2
snpeff_version_short=$3
genome_version=$4

snpeff_folder=${conda_home}/share/snpeff-${snpeff_version}

sudo mkdir -p ${conda_home}/share/snpeff-${snpeff_version}/data/${genome_version}
sudo wget http://downloads.sourceforge.net/project/snpeff/databases/v${snpeff_version_short}/snpEff_v${snpeff_version_short}_${genome_version}.zip -P ${snpeff_folder}/
sudo rm -rf ${snpeff_folder}/data/
sudo chmod 777 -R ${snpeff_folder}
cd ${snpeff_folder}
sudo unzip snpEff_v${snpeff_version_short}_${genome_version}.zip
sudo ${conda_home}/bin/java -jar ${snpeff_folder}/snpEff.jar download ${genome_version}
#example: ./get_snpeff_data.bash /root/conda 4.3.1t-1 4_3 GRCh37.75