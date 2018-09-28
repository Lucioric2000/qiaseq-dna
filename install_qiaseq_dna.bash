#!/bin/bash
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
sudo yum install git unzip cpan wget gcc bzip2 python-devel
srv_qiagen=/srv/qgen
sudo mkdir ${srv_qiagen}
sudo chmod 777 ${srv_qiagen}
cd ${srv_qiagen}

qseqdnamatch=`expr match "$(pwd)" '.*\(qiaseq-dna\)'`
if [[ $qseqdnamatch -eq "qiaseq-dna" ]]; then
    echo "Already in qiaseq-dna folder."
else;
    echo "Not in qiaseq-dna folder."
    git clone --recursive https://github.com/Lucioric2000/qiaseq-dna
    cd qiaseq-dna
fi
#Sets up a script with the environment variables needed
#To uninstall:
#sudo rm -rf /srv/qgen; sudo rm -rf /opt/conda
#Declare the location of the conda installation
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
condabin=`which conda`
if [ -z $condabin ]
then
    conda_home=/opt/conda
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
source activate base
# Picard 1.97 was not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version
sudo ${conda_home}/bin/conda install -c cyclus java-jdk=8.45.14
sudo ${conda_home}/bin/conda install openpyxl
sudo ${conda_home}/bin/pip install --upgrade pip
sudo ${conda_home}/bin/pip install statistics
#conda install bedtools=2.25.0 htslib=1.3.1 cutadapt=1.10 picard=1.97 snpeff=4.2 bwa=0.7.15 pysam=0.9.0 java-jdk=8.45.14 samtools 1.5
sudo ${conda_home}/bin/conda install -c bioconda bedtools htslib cutadapt picard snpeff snpsift bwa pysam samtools biopython rstudio samtools scipy MySQL-python

#Output the contents of ~/.bashrc plus the content enclosed in qoutes (which is in a string representation that handles newline characters) to the file ~/.bashrc.new.qiagen
#echo -e "\n#Shell environment for qiagen\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${srv_qiagen}/bin/ssw/src/"|cat ~/.bashrc ->~/.bashrc.new.qiagen
#Move the file ~/.bashrc.new.qiagen to ~/.bashrc (overwriting the existent ~/.bashrc withouk asking for confirmation)
#mv -f ~/.bashrc.new.qiagen ~/.bashrc
##Calls this script
#bash -c /etc/profile.d/qiagen.sh
#Make the directories if don't exist
mkdir -p ${srv_qiagen}/bin/downloads && \
    mkdir -p ${srv_qiagen}/data/genome && \
    mkdir -p ${srv_qiagen}/data/annotation && \
    mkdir -p ${srv_qiagen}/example/

#Update the updatable software in Debian repositories
#sudo apt-get -y update

################ Install various version specific 3rd party tools ################

wget https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar -P ${srv_qiagen}/bin/
cd ${srv_qiagen}/bin/ && tar -xvf ssw.tar.gz
sudo bash -c "echo '/srv/qgen/bin/ssw/src'>/etc/ld.so.conf.d/ssw.conf"
#Import the configuration files into the system
sudo ldconfig

## Download and install 3rd party libraries
wget https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz -P ${srv_qiagen}/bin/downloads/
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf py-editdist-0.3.tar.gz && \
    cd py-editdist-0.3 && python setup.py install && \
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf sendgrid-v2.2.1.tar.gz && \
    cd sendgrid-python-2.2.1 && python setup.py install

################ R packages ################
sudo echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > /root/.Rprofile
sudo ${conda_home}/bin/Rscript -e "install.packages('plyr')"
sudo ${conda_home}/bin/Rscript -e "install.packages('MASS')"
sudo ${conda_home}/bin/Rscript -e "install.packages('ggplot2')"
sudo ${conda_home}/bin/Rscript -e "install.packages('gridExtra')"
sudo ${conda_home}/bin/Rscript -e "install.packages('naturalsort')"
sudo ${conda_home}/bin/Rscript -e "install.packages('scales')"
sudo ${conda_home}/bin/Rscript -e "install.packages('extrafont')"

## Perl
sudo cpan Module::Runtime
sodo cpan XML::Twig DateTime DBI DBD::SQLite Env::Path File::chdir Getopt::Long::Descriptive Sort:Naturally Config::IniFiles Data::Dump::Color Data::Table::Excel Hash::Merge File::Slurp

################ TVC binaries ################
mkdir -p ${srv_qiagen}/bin/TorrentSuite/
wget https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tmap \
         https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tvc \
     -P ${srv_qiagen}/bin/TorrentSuite/
chmod 775 ${srv_qiagen}/bin/TorrentSuite/tmap ${srv_qiagen}/bin/TorrentSuite/tvc

################ Add data directory ################

## Download Annotation files
wget https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz \
         https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz.tbi \
         https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz \
         https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz.tbi \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz.tbi \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat_TRF.bed \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL_RepeatMasker.bed \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/bkg.error.v2.RData \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL.full.bed \
     https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat.full.bed \
      -P ${srv_qiagen}/data/annotation/

## Download annotations for SnpEff (4.3)
sudo mkdir -p ${conda_home}/share/snpeff-4.3.1t-1/data/GRCh37.75
sudo wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip -P ${conda_home}/share/snpeff-4.3.1t-1/
sudo rm -rf ${conda_home}/share/snpeff-4.3.1t-1/data/
cd ${conda_home}/share/snpeff-4.3.1t-1/
sudo unzip snpEff_v4_3_GRCh37.75.zip
/opt/conda/jre/bin/java -jar /opt/conda/share/snpeff-4.3.1t-1/snpEff.jar download GRCh37.75

## Annotation file
wget https://storage.googleapis.com/qiaseq-dna/data/annotation/refGene.txt \
         -P ${srv_qiagen}/data/annotation/


## Add example fastqs and files
wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz \
         https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz \
     https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt  \
     https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed \
         -P ${srv_qiagen}/example/

## Add test files for smCounterv2
mkdir -p ${srv_qiagen}/test_smcounter-v2/
wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed \
         https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam \
     https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
     https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai \
     -P ${srv_qiagen}/test_smcounter-v2/

## Download genome files
wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P ${srv_qiagen}/data/genome/
cd ${srv_qiagen}/data/genome && \
    gunzip ucsc.hg19.fa.gz  && \
    ## Index the fasta using samtools
    samtools faidx ${srv_qiagen}/data/genome/ucsc.hg19.fa && \ 
    ${conda_home}/bin/bwa index ${srv_qiagen}/data/genome/ucsc.hg19.fa


## Annotation file
wget https://storage.googleapis.com/qiaseq-dna/data/annotation/refGene.txt \
         -P ${srv_qiagen}/data/annotation/



#time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 single out1 NEB_S2 &> run_v1.log &
#time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2 NEB_S2 &> run_v2.log &
#python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 tumor-normal tumor_readset normal_readset > run_v1_tn.log 2>&1 &

