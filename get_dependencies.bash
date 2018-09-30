#!/bin/bash
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
sudo yum install git unzip cpan wget gcc bzip2 python-devel nano expat-devel openssl-devel
srv_qiagen=/srv/qgen
sudo mkdir ${srv_qiagen}
sudo chmod -R 777 ${srv_qiagen}
cd ${srv_qiagen}

qseqdnamatch=`expr match "$(pwd)" '.*\(qiaseq-dna\)'`
if [[ $qseqdnamatch -eq "qiaseq-dna" ]]
then
    echo "Already in qiaseq-dna folder."
else
    echo "Not in qiaseq-dna folder."
    git clone --recursive https://github.com/Lucioric2000/qiaseq-dna
    cd qiaseq-dna
fi
#Sets up a script with the environment variables needed
#To uninstall:
#sudo rm -rf /srv/qgen; sudo rm -rf /opt/conda; sudo rm -rf /srv/conda
#Declare the location of the conda installation
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
condabin=`which conda`
if [ -z $condabin ]
then
    conda_home=/opt/conda
    conda_home=/srv/conda
    #Install the Miniconda Python pachages manager
    echo "Next, the Miniconda package will be downloaded and installed"
    echo "You should install it at the default location shown"
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    chmod +x Miniconda2-latest-Linux-x86_64.sh
    sudo sh Miniconda2-latest-Linux-x86_64.sh -p $conda_home -u
    rm Miniconda2-latest-Linux-x86_64.sh
    #Make the updated shell path available in this session:
    source ~/.bashrc
else
    conda_home=${condabin%/bin/conda}
    echo "Conda installation found at $conda_home. Script will use tht installation."
fi

################ Install python modules ################
## Install some modules with conda
#This includes R (rstudio) and biopython
source ${conda_home}/bin/activate base
# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version
sudo ${conda_home}/bin/conda install -c cyclus java-jdk=8.45.14
sudo ${conda_home}/bin/conda install openpyxl gxx_linux-64
sudo ${conda_home}/bin/pip install --upgrade pip
sudo ${conda_home}/bin/pip install statistics msgpack-python python_http_client==1.2.3 smtpapi==0.3.1 PyHamcrest==1.9.0
#conda install bedtools=2.25.0 htslib=1.3.1 cutadapt=1.10 picard=1.97 snpeff=4.2 bwa=0.7.15 pysam=0.9.0 java-jdk=8.45.14 samtools 1.5
sudo ${conda_home}/bin/conda install -c bioconda bedtools htslib cutadapt picard snpeff snpsift bwa pysam samtools biopython rstudio r-essentials r-mass r-scales r-extrafont r-plyr samtools scipy MySQL-python

mkdir -p ${srv_qiagen}/bin/downloads && mkdir -p ${srv_qiagen}/data/genome && mkdir -p ${srv_qiagen}/data/annotation && mkdir -p ${srv_qiagen}/example/

################ Install various version specific 3rd party tools ################

wget https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar -P ${srv_qiagen}/bin/
cd ${srv_qiagen}/bin/ && tar -xvf ssw.tar.gz
sudo bash -c "echo '/srv/qgen/bin/ssw/src'>/etc/ld.so.conf.d/ssw.conf"
#Import the configuration files into the system
sudo ldconfig

## Download and install 3rd party libraries
wget https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz -P ${srv_qiagen}/bin/downloads/
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf py-editdist-0.3.tar.gz && \
    cd py-editdist-0.3 && sudo ${conda_home}/bin/python setup.py install && \
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf sendgrid-v2.2.1.tar.gz && \
    cd sendgrid-python-2.2.1 && sudo ${conda_home}/bin/python setup.py install

################ R packages ################
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" | sudo tee /root/.Rprofile >/dev/null
#sudo ${conda_home}/bin/Rscript -e "install.packages('plyr')"
#sudo ${conda_home}/bin/Rscript -e "install.packages('MASS')"
sudo ${conda_home}/bin/Rscript -e "install.packages('ggplot2')"
sudo ${conda_home}/bin/Rscript -e "install.packages('gridExtra')"
sudo ${conda_home}/bin/Rscript -e "install.packages('naturalsort')"
#sudo ${conda_home}/bin/Rscript -e "install.packages('scales')"
#sudo ${conda_home}/bin/Rscript -e "install.packages('extrafont')"

## Perl
sudo cpan install CPAN
sudo cpan reload cpan
sudo cpan install CPAN::Meta CPAN::Meta::YAML ExtUtils::CBuilder Module::Metadata Parse::CPAN::Meta Perl::OSType TAP::Harness JSON::PP

sudo cpan install Module::Runtime HTTP::Date Test::Pod XML::Twig
sudo cpan install IO::Socket::SSL DateTime DBI DBD::SQLite Env::Path File::chdir Getopt::Long::Descriptive Sort:Naturally Config::IniFiles Data::Dump::Color Data::Table::Excel Hash::Merge File::Slurp

################ TVC binaries ################
mkdir -p ${srv_qiagen}/bin/TorrentSuite/
wget https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tmap \
         https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tvc \
     -P ${srv_qiagen}/bin/TorrentSuite/
chmod 775 ${srv_qiagen}/bin/TorrentSuite/tmap ${srv_qiagen}/bin/TorrentSuite/tvc
