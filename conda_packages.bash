#!/bin/bash
################ Install python modules ################
conda_home=/srv/conda
## Install some modules with conda
#This includes R (rstudio) and biopython
source ${conda_home}/bin/activate base
${conda_home}/bin/conda config --add channels defaults
${conda_home}/bin/conda config --add channels bioconda
${conda_home}/bin/conda config --add channels conda-forge

# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version
sudo ${conda_home}/bin/pip install --upgrade pip
sudo ${conda_home}/bin/pip install statistics msgpack-python python_http_client==1.2.3 smtpapi==0.3.1 PyHamcrest==1.9.0
#sudo ${conda_home}/bin/conda install -c tsnyder gcc_linux-cos6-x86_64 tsnyder gxx_linux-cos6-x86_64 tsnyder binutils_linux-cos6-x86_64
sudo ${conda_home}/bin/conda install -c bioconda perl-app-cpanminus
sudo ${conda_home}/bin/conda install -c cyclus java-jdk=8.45.14
sudo ${conda_home}/bin/conda install openpyxl gxx_linux-64 gcc_linux-64 gfortran_linux-64
sudo ${conda_home}/bin/conda install -c bioconda bedtools htslib cutadapt picard snpeff snpsift bwa pysam samtools biopython rstudio r-essentials r-mass r-scales r-extrafont r-plyr samtools scipy MySQL-python
