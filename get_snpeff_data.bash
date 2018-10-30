#!/bin/bash
## ${conda_home} 4.3.1t-1 4_3 GRCh37.75
## Download annotations for SnpEff (4.3)
#conda_home=/srv/conda
#snpeff_version=4.3.1t-1
#snpeff_version_short=4_3
#genome_version=GRCh37.75
conda_home=$1
snpeff_version=$2
snpeff_version_short=$3
genome_version=$4

sudo mkdir -p ${conda_home}/share/snpeff-${snpeff_version}/data/${genome_version}
sudo wget http://downloads.sourceforge.net/project/snpeff/databases/v${snpeff_version_short}/snpEff_v${snpeff_version_short}_${genome_version}.zip -P ${conda_home}/share/snpeff-${snpeff_version}/
sudo rm -rf ${conda_home}/share/snpeff-${snpeff_version}/data/
cd ${conda_home}/share/snpeff-${snpeff_version}/
sudo unzip snpEff_v${snpeff_version_short}_${genome_version}.zip
#sudo chmod 777 ${conda_home}/share/snpeff-${snpeff_version}
sudo ${conda_home}/jre/bin/java -jar ${conda_home}/share/snpeff-${snpeff_version}/snpEff.jar download ${genome_version}