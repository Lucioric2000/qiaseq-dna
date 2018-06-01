#!/bin/bash
#Install the Miniconda Python pachages manager
sudo yum install git unzip r-base
echo "Next, the Miniconda package will be downloaded and installed"
echo "You should install it as the miniconda3 subdirectory of your home directory"
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
sh Miniconda2-latest-Linux-x86_64.sh
#rm Miniconda2-latest-Linux-x86_64.sh
#Make the updated shell path available in this session
#source ~/.bashrc

#Code for installing the qiagen-dna software and the example in BASH
#Sets up a script with the environment variables needed
srv_qiagen=/srv/qiagen
sudo echo -e "#Shell environment for qiagen\nexport PYTHONPATH=$PYTHONPATH:${srv_qiagen}/code/qiaseq-dna/\nexport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${srv_qiagen}/bin/ssw/src/">/etc/profile.d/qiagen.sh
#Calls thiss script
bash -c /etc/profile.d/qiagen.sh
#Declare the location of the conda installaction
conda_home=~/miniconda2
#Make the directories if don't exist
sudo mkdir -p ${srv_qiagen}/code && \
    sudo mkdir -p ${srv_qiagen}/bin/downloads && \
    sudo mkdir -p ${srv_qiagen}/data/genome && \
    sudo mkdir -p ${srv_qiagen}/data/annotation && \
    sudo mkdir -p ${srv_qiagen}/example/

#Update the updatable software in Debian repositories
sudo apt-get -y update
#Install R
#sudo apt-get -y install r-base

################ Install various version specific 3rd party tools ################
#This includes R (rstudio) and biopython
#Pysam in the dockerfile is said to be 0.9.0, but the version that colud be installed is 0.6
conda install -c bioconda bedtools=2.25.0 htslib=1.3.1 cutadapt=1.10 snpeff=4.2 bwa=0.7.15 rstudio biopython pysam
#Picard 1.97 wa not found in the default conda ditribution
################ Update openjdk ################
## note : picard gets updated to match jdk version
conda install -c cyclus java-jdk=8.45.14
################ Install python modules ################
## Install some modules with conda
conda install scipy MySQL-python openpyxl pysam=0.9.0
pip install statistics

wget https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar -P ${srv_qiagen}/bin/
cd ${srv_qiagen}/bin/ && tar -xvf ssw.tar.gz    

## Download and install 3rd party libraries
wget https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz -P ${srv_qiagen}/bin/downloads/
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf py-editdist-0.3.tar.gz && \
    cd py-editdist-0.3 && python setup.py install && \
    cd ${srv_qiagen}/bin/downloads/ && tar -xvf sendgrid-v2.2.1.tar.gz && \
    cd sendgrid-python-2.2.1 && python setup.py install

################ R packages ################
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
sudo Rscript -e "install.packages('plyr')"


################ Add latest samtools version for sort by Tag feature ################
wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -O ${srv_qiagen}/bin/downloads/samtools-1.5.tar.bz2 && \
    cd ${srv_qiagen}/bin/downloads/ tar -xvf samtools-1.5.tar.bz2 && \
    cd samtools-1.5  && mkdir -p ${srv_qiagen}/bin/samtools-1.5 && ./configure --prefix ${srv_qiagen}/bin/samtools-1.5 && make && make install

################ Add data directory ################
## Download genome files
sudo wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P ${srv_qiagen}/data/genome/
cd ${srv_qiagen}/data/genome && \
    sudo gunzip ucsc.hg19.fa.gz  && \
    ## Index the fasta using samtools
    ${srv_qiagen}/bin/samtools-1.5/bin/samtools faidx ${srv_qiagen}/data/genome/ucsc.hg19.fa && \ 
    ## Run bwa to generate index files 
    ${conda_home}/bin/bwa index ${srv_qiagen}/data/genome/ucsc.hg19.fa

## Download Annotation files
sudo wget https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz \
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

## Download annotation using SnpEff command
sudo wget http://downloads.sourceforge.net/project/snpeff/databases/v4_2/snpEff_v4_2_GRCh37.75.zip -P ${conda_home}/share/snpeff-4.2-0/
rm -rf ${conda_home}/share/snpeff-4.2-0/data/
cd ${conda_home}/share/snpeff-4.2-0/
sudo unzip snpEff_v4_2_GRCh37.75.zip
#Frees up some space
rm snpEff_v4_2_GRCh37.75.zip

## The command below is not working anymore because of some certificate issue (debug later)
#RUN /opt/conda/jre/bin/java -jar ${conda_home}/share/snpeff-4.2-0/snpEff.jar download GRCh37.75
################ Modules for CNV Analysis ################
## Perl
sudo cpan DateTime
sudo cpan DBI
sudo cpan DBD::SQLite
sudo cpan Env::Path
sudo cpan File::chdir
sudo cpan Getopt::Long::Descriptive
sudo cpan Sort:Naturally
sudo cpan Config::IniFiles
sudo cpan Data::Dump::Color
sudo cpan Data::Table::Excel
sudo cpan Hash::Merge
sudo cpan File::Slurp

## R
sudo Rscript -e "install.packages('MASS')"
sudo Rscript -e "install.packages('ggplot2')"
sudo Rscript -e "install.packages('gridExtra')"
sudo Rscript -e "install.packages('naturalsort')"
sudo Rscript -e "install.packages('scales')"
sudo Rscript -e "install.packages('ggplot2')"
sudo Rscript -e "install.packages('extrafont')"
## Annotation file
sudo wget https://storage.googleapis.com/qiaseq-dna/data/annotation/refGene.txt \
         -P ${srv_qiagen}/data/annotation/

################ TVC binaries ################
sudo mkdir -p ${srv_qiagen}/bin/TorrentSuite/
sudo wget https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tmap \
    	 https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tvc \
	 -P ${srv_qiagen}/bin/TorrentSuite/
sudo chmod 775 ${srv_qiagen}/bin/TorrentSuite/tmap ${srv_qiagen}/bin/TorrentSuite/tvc


## Add example fastqs and files
sudo wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz \
    	 https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz \
	 https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt  \
	 https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed \
    	 -P ${srv_qiagen}/example/

## Add test files for smCounterv2
sudo mkdir -p ${srv_qiagen}/test_smcounter-v2/
sudo wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed \
         https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai \
	 -P ${srv_qiagen}/test_smcounter-v2/
