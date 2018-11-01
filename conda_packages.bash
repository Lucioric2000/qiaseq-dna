#!/bin/bash
################ Install python modules ################
conda_home=$1
conda_preffix=$2
conda_env=$3
#conda_home=/srv/conda
#conda_env=base
## Install some modules with conda
#This includes R (rstudio) and biopython
source ${conda_home}/bin/activate base
${conda_home}/bin/conda config --add channels defaults
${conda_home}/bin/conda config --add channels bioconda
${conda_home}/bin/conda config --add channels conda-forge

# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version. Thus picard and jdk version should be pinned, at least by picard version is hardcoded in a configuration file
sudo ${conda_preffix}/bin/pip install --upgrade pip
sudo ${conda_preffix}/bin/pip install statistics msgpack-python python_http_client==1.2.3 smtpapi==0.3.1 PyHamcrest==1.9.0
#sudo ${conda_home}/bin/conda install -n ${conda_env} -c tsnyder gcc_linux-cos6-x86_64 tsnyder gxx_linux-cos6-x86_64 tsnyder binutils_linux-cos6-x86_64
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda perl perl-cpan-shell perl-app-cpanminus
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c cyclus java-jdk=8.45.14
sudo ${conda_home}/bin/conda install -y -n ${conda_env} openpyxl gxx_linux-64 gcc_linux-64 gfortran_linux-64
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c bioconda r-ggplot2 bedtools htslib cutadapt picard=2.18.15 snpeff snpsift bwa pysam samtools biopython rstudio r-essentials r-mass r-scales r-extrafont r-plyr samtools scipy MySQL-python
sudo ${conda_home}/bin/conda install -y -n ${conda_env} -c conda-forge r-gridextra r-naturalsort
