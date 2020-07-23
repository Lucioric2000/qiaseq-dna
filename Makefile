VERSION := $(shell cat version.dat)
SOURCE=qiaseq-dna
qiagen_parent_folder=/srv/qgen
qiagen_code_parent_folder=/srv/qgen/code
conda_home=/root/conda
conda_env=base
ROOT=/root
DISTRO := $(shell bash -c "yum --help&>/dev/null && echo centos || echo ubuntu")
GENOME_BUILD=GRCh37.75
GENOME_BUILD_ALT_NAME=hg19
#4.3.1t-3 4_3
SNPEFF_VERSION=4.3.1t-0
SNPEFF_VERSION_SHORT_STRING=4_3 
SNPSIFT_VERSION=4.3.1t-2
SNPEFF_VERSION_WO_BUILD=4.3.1t
SNPSIFT_VERSION_WO_BUILD=4.3.1t
#SNPEFF_VERSION=4.5covid19-1
#SNPEFF_VERSION_SHORT_STRING=4_5_covid19_1 
#SNPEFF_VERSION_SHORT_STRING=4_5_covid19_1 
archive:
	sudo rm -f $$p/$(SOURCE)-*.tar.gz $$p/Miniconda*
	p=`pwd` && rm -f $$p/$(SOURCE)-$(VERSION).tar.gz && tar --transform="s@^@$(SOURCE)-$(VERSION)/@" -cvzf $$p/$(SOURCE)-$(VERSION).tar.gz *
version:
	@echo $(VERSION)
libraries_centos:
	sudo yum -y install git unzip cpan wget gcc gcc-c++ bzip2 python2-devel nano expat-devel openssl-devel perl perl-CPAN perl-devel curl perl-App-cpanminus python3 python3-pip python3-libs python3-tools python3-devel
	#sudo pip3 install edlib
libraries_ubuntu:
	sudo apt-get install -y git unzip wget gcc build-essential bzip2 python2-dev nano perl curl python3 cpanminus python3-pip python3-dev #expat-dev openssl-dev perl-dev perl-App-cpanminus python3-libs python3-tools
	#sudo pip3 install edlib
toqgz: archive
	cp ./$(SOURCE)-$(VERSION).tar.gz ./install_$(SOURCE)-v$(VERSION).bash $(qiagen_parent_folder)
toqiaseq: archive
	cp ./$(SOURCE)-$(VERSION).tar.gz ./install_$(SOURCE)-v$(VERSION).bash $(qiagen_parent_folder)
	cd $(qiagen_parent_folder) && tar -xvzf $(qiagen_parent_folder)/$(SOURCE)-$(VERSION).tar.gz && sudo rm -rf $(SOURCE)-old
	if [ -e $(qiagen_parent_folder)/$(SOURCE) ]; then sudo mv -f $(qiagen_parent_folder)/$(SOURCE) $(qiagen_parent_folder)/$(SOURCE)-old; fi
	cd $(qiagen_parent_folder) && sudo mv $(SOURCE)-$(VERSION) $(SOURCE)
toroot: archive
	sudo cp ./$(SOURCE)-$(VERSION).tar.gz ./install_$(SOURCE)-v$(VERSION).bash $(ROOT)/
	sudo chmod +x $(ROOT)/install_$(SOURCE)-v$(VERSION).bash
takeback:
	sudo cp -R ${qiagen_parent_folder}/${SOURCE} ${qiagen_code_parent_folder}
update:
	mkdir -p $(qiagen_parent_folder)/qiaseq-dna
	cp -rf $(qiagen_parent_folder)/code/qiaseq-dna/* $(qiagen_parent_folder)/qiaseq-dna/

install:
	make libraries_$(DISTRO)
	make conda_install
	#make install_python27_env_if_needed
	make modules_and_snpeff
	make thirdparty_tools
	make nirvana
	make data_files
	make help
nirvana:
	bash ./install_Nirvana.bash
conda_install:# clean
	#Install the Miniconda Python pachages manager
	echo "Next, the Miniconda package will be downloaded and installed at the folder ${conda_home}."
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	chmod +x Miniconda3-latest-Linux-x86_64.sh
	sudo bash Miniconda3-latest-Linux-x86_64.sh -p ${conda_home} -u -b
	rm Miniconda3-latest-Linux-x86_64.sh
	${conda_home}/bin/conda init bash
	#Make the updated shell path available in this session:
	#source ~/.bashrc
	#source ${conda_home}/bin/activate ${conda_env}
	echo "Conda was installed in the ${conda_home} folder. The environment that will be used is ${conda_env}."
install_python27_env_if_needed:
	source ${conda_home}/bin/activate ${conda_env} && (
	#echo Activated conda environment ${conda_env}
	apyversion=$(python -c "import sys;print(sys.version)")
	echo "Python version: ${apyversion}"
	if [[ $apyversion =~ ^2.7 ]]
	then
	    #echo Python version is 2.7 in the base conda environment
	    echo base
	else
	    conda_env=python2.7
	    #Conda environment not found: creating it
	    sudo ${conda_home}/bin/conda create -n ${conda_env} python=2.7;
	    source ${conda_home}/bin/activate ${conda_env};
	    #echo Created and activated the conda environment ${conda_env}
	    echo python2.7
	fi
	 ) || ( 
	#Conda environment not found: creating it
	sudo ${conda_home}/bin/conda create -n ${conda_env} python=2.7;
	source ${conda_home}/bin/activate ${conda_env};
	#echo Created and activated the conda environment ${conda_env}
	echo base
	)

modules_and_snpeff:
	################ Install python modules ################
	## Install some modules with conda
	#This includes R (rstudio) and biopython
	bash -c "source ${conda_home}/bin/activate" && ./conda_packages.bash ${conda_home} ${conda_home} ${conda_env} ${SNPEFF_VERSION_WO_BUILD} ${SNPSIFT_VERSION_WO_BUILD}
	./install_perl_modules.bash ${conda_home} ${conda_home} ${conda_env}
	#Install Snpeff, Snpsift and their data
	./get_snpeff_data.bash ${conda_home} ${SNPEFF_VERSION} ${SNPEFF_VERSION_SHORT_STRING} ${GENOME_BUILD} #If you wanted to use the GRCh38, you should replace GRCh37.75 to GRCh38.86 in this line

thirdparty_tools:
	################ Install various version specific 3rd party tools ################
	wget https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar -P ${qiagen_parent_folder}/bin/
	cd ${qiagen_parent_folder}/bin/ && tar -xvf ssw.tar.gz
	sudo bash -c "echo '${qiagen_parent_folder}/bin/ssw/src'>/etc/ld.so.conf.d/ssw.conf"
	#Import the configuration files into the system
	sudo ldconfig
	mkdir -p ${qiagen_parent_folder}/bin/downloads
	wget https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz -P ${qiagen_parent_folder}/bin/downloads/
	cd ${qiagen_parent_folder}/bin/downloads/ && tar -xvf py-editdist-0.3.tar.gz && \
	cd py-editdist-0.3 && sudo ${conda_home}/bin/python setup.py install && \
	cd ${qiagen_parent_folder}/bin/downloads/ && tar -xvf sendgrid-v2.2.1.tar.gz && \
	cd sendgrid-python-2.2.1 && sudo ${conda_home}/bin/python setup.py install
	################ TVC binaries ################
	mkdir -p ${qiagen_parent_folder}/bin/TorrentSuite/
	wget https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tmap \
	https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tvc \
	-P ${qiagen_parent_folder}/bin/TorrentSuite/
	chmod 775 ${qiagen_parent_folder}/bin/TorrentSuite/tmap ${qiagen_parent_folder}/bin/TorrentSuite/tvc

#Data files:
data_files: annotations examples testfiles genomes
annotations:
	## Download Annotation files
	mkdir -p ${qiagen_parent_folder}/data/annotation
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
	https://storage.googleapis.com/qiaseq-dna/data/annotation/refGene.txt \
	-P ${qiagen_parent_folder}/data/annotation/

examples:
	## Add example fastqs and files
	mkdir -p ${qiagen_parent_folder}/example/
	wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz -P ${qiagen_parent_folder}/example/
	wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz -P ${qiagen_parent_folder}/example/
	wget https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt -P ${qiagen_parent_folder}/example/
	wget https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed -P ${qiagen_parent_folder}/example/
	ln -f ${qiagen_parent_folder}/example/NEB_S2_L001_R1_001.fastq.gz ${qiagen_parent_folder}/example/Cellline_S10_L001_R1_001.fastq.gz
	ln -f ${qiagen_parent_folder}/example/NEB_S2_L001_R2_001.fastq.gz ${qiagen_parent_folder}/example/Cellline_S10_L001_R2_001.fastq.gz

testfiles:
	## Add test files for smCounterv2
	mkdir -p ${qiagen_parent_folder}/test_smcounter-v2/
	wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed -P ${qiagen_parent_folder}/test_smcounter-v2/
	wget https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam -P ${qiagen_parent_folder}/test_smcounter-v2/
	wget https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt -P ${qiagen_parent_folder}/test_smcounter-v2/
	wget https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai -P ${qiagen_parent_folder}/test_smcounter-v2/

genomes:
	#Index the genome fasta file, using samtools and bwa, only if does not exists a file with a md5 hash identical to a hash annotated in a file generated after
	#a successful bwa run below
	bash -c "./genome.bash ${qiagen_parent_folder} ${conda_home} base ${GENOME_BUILD_ALT_NAME}"

#help and debug:
versions:
	python --version
help:
	#echo "Smcounter2 is now installed."
	@echo "To run smcounter (version 1 or 2), you should first:"
	@echo "1. activate the Anaconda enviroment ${conda_env}, using the command source ${conda_home}/bin/activate ${conda_env}"
	@echo "2. Go to the directory ${qiagen_parent_folder}/qiaseq-dna"
	@echo "After that, you can execute a smcounter using a command of this form:"
	@echo "time python run_qiaseq_dna.py <config file> <v1/v2> <single/tumor-normal> <output file path> <sample1>[ <sample2>[ ...<samplen>]] &> <log_file> &"
	@echo "To see a description of the command line options (except time and &> <log_file> &) you may execute the command:"
	@echo "python run_qiaseq_dna.py --help (from the directory ${qiagen_parent_folder}/qiaseq-dna and with the conda environment ${conda_env} activated)."
	@echo "Now such help message will be displayed:"
	@echo cd ${qiagen_parent_folder}/qiaseq-dna && bash -c "source ${conda_home}/bin/activate && python run_qiaseq_dna.py --help"
	@echo you may run smcounter with commands like that:
	@echo
	@echo 1. Smcounterv1:
	@echo "time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 single out1 NEB_S2 &> run_smcounterv1.log &"
	@echo
	@echo 2. Smcounterv2:
	@echo "time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2 NEB_S2 &> run_smcounterv2.log &"
	@echo
	@echo 3. Tumor-normal analysis:
	@echo "time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 tumor-normal tumor_readset normal_readset > run_smcounterv1_tn.log 2>&1 &"
	@echo
	@echo 4. Multiple samples:
	@echo "time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2_multisamples_{0} sample1 sample2 sample3 (...) samplen &> run_v6.log &"
	@echo For instance:
	@echo "time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2_multisamples_{0}_{1} NEB_S2 NEB_S2 &> run_smcounterv2_multiplesamples.log &"
test1:
	cd $(qiagen_parent_folder)/$(SOURCE) && bash -c "source ${conda_home}/bin/activate && time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 single out_smcounterv1 NEB_S2 &> run_smcounterv1.log"
test2single:
	cd $(qiagen_parent_folder)/$(SOURCE) && bash -c "source ${conda_home}/bin/activate; time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out_smcounterv2 NEB_S2 &> run_smcounterv2.log"
testtn:
	cd $(qiagen_parent_folder)/$(SOURCE) && bash -c "source ${conda_home}/bin/activate; time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 tumor-normal tumor_readset normal_readset > run_smcounterv1_tn.log 2>&1"
testmulti:
	cd $(qiagen_parent_folder)/$(SOURCE) && bash -c "source ${conda_home}/bin/activate base  && time python run_qiaseq_dna.py forcellline.txt v2 single out_smcounterv2_{0}_{1}_cellline NEB_S2 Cellline_S10 &> run_smcounterv2_multisamples_cellline.log"
testx:
	cd $(qiagen_parent_folder)/$(SOURCE) && bash -c "source ${conda_home}/bin/activate; pwd; which python"


#Cleaning
procclean:
	pkill python
	pkill conda
clean:
	sudo rm -rf ${conda_home} /srv/qgen/conda ~/.conda ~/.condarc ~/.continuum /srv/qgen/qiaseq-dna
condaclean:
	sudo rm -rf ${conda_home} /srv/qgen/conda ~/.conda ~/.condarc ~/.continuum

