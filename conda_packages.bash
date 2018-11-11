#!/bin/bash
################ Install python modules ################
if [[ -n $1 ]];
then
	conda_home=$1
	conda_preffix=$2
	conda_env=$3
else
	conda_home=/srv/conda
	conda_preffix=/srv/conda
	conda_env=base
fi
## Install some modules with conda
#This includes R (rstudio) and biopython
echo Installing conda packages in environment ${conda_env}, which reside in the folder ${conda_preffix}.
source ${conda_home}/bin/activate ${conda_env}
${conda_home}/bin/conda config --add channels defaults
${conda_home}/bin/conda config --add channels bioconda
${conda_home}/bin/conda config --add channels conda-forge

# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version. Thus picard and jdk version should be pinned, at least by picard version is hardcoded in a configuration file
sudo ${conda_preffix}/bin/pip install --upgrade pip
sudo ${conda_preffix}/bin/pip install statistics msgpack-python python_http_client==1.2.3 smtpapi==0.3.1 PyHamcrest==1.9.0
#sudo ${conda_home}/bin/conda install -n ${conda_env} -c tsnyder gcc_linux-cos6-x86_64 tsnyder gxx_linux-cos6-x86_64 tsnyder binutils_linux-cos6-x86_64
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda perl bwa pysam samtools biopython scipy MySQL-python bedtools htslib cutadapt picard=2.18.15 snpeff snpsift
#sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda perl-cpan-shell perl-app-cpanminus
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c cyclus java-jdk=8.45.14
sudo ${conda_home}/bin/conda install -y -n ${conda_env} openpyxl gxx_linux-64 gcc_linux-64 gfortran_linux-64
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda r-ggplot2 rstudio r-essentials r-mass r-scales r-extrafont r-plyr
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c conda-forge r-gridextra r-naturalsort
