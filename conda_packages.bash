#!/bin/bash
################ Install python modules ################
if [[ -n $1 ]];
then
	conda_home=$1
	conda_preffix=$2
	conda_env=$3
	snpeff_version=$4
	snpsift_version=$5
else
	conda_home=/root/conda
	conda_preffix=/root/conda
	conda_env=base
	snpeff_version=4.3.1t
	snpsift_version=4.3.1t
fi
## Install some modules with conda
#This includes R (rstudio) and biopython
echo Installing conda packages in environment ${conda_env}, which reside in the folder ${conda_home}.
source ${conda_home}/bin/activate ${conda_env}
#${conda_home}/bin/conda config --add channels defaults
#${conda_home}/bin/conda config --add channels bioconda
#${conda_home}/bin/conda config --add channels conda-forge

# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version. Thus picard and jdk version should be pinned, at least by picard version is hardcoded in a configuration file
sudo ${conda_home}/bin/conda install -y -n ${conda_env} python=3.7
sudo ${conda_home}/bin/pip install cython edlib python-Levenshtein
sudo ${conda_home}/bin/pip install --upgrade pip
sudo ${conda_home}/bin/pip install statistics msgpack-python python_http_client==1.2.3 smtpapi==0.3.1 PyHamcrest==1.9.0 openpyxl edlib
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda perl bwa pysam samtools biopython scipy bedtools htslib cutadapt picard snpeff=${snpeff_version} snpsift #=${snpsift_version} MySQL-python
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c conda-forge r-ggplot2 r-essentials r-mass r-scales r-extrafont r-plyr r-gridextra r-naturalsort r-tidyverse r-stringi pyicu
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c cyclus java-jre
#sudo ${conda_home}/bin/conda create -y -n python37 python=3.7 cython samtools
#sudo ${conda_home}/bin/conda install -y -n python37 -c bioconda samtools
#sudo ${conda_home}/envs/python37/bin/pip install edlib
