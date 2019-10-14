VERSION=14.0
SOURCE=qiaseq-dna
qiagen_parent_folder=/srv/qgen
conda_home=/root/conda
conda_env=base

archive:
	sudo rm -f $$p/$(SOURCE)-*.tar.gz $$p/Miniconda*
	p=`pwd` && rm -f $$p/$(SOURCE)-$(VERSION).tar.gz && tar --transform="s@^@$(SOURCE)-$(VERSION)/@" -cvzf $$p/$(SOURCE)-$(VERSION).tar.gz *
libraries:
	sudo yum -y install git unzip cpan wget gcc gcc-c++ bzip2 python-devel nano expat-devel openssl-devel perl perl-CPAN perl-devel curl gcc perl-App-cpanminus
toroot: archive
	cp ./$(SOURCE)-$(VERSION).tar.gz ./install_${SOURCE}-v${VERSION}.bash ${qiagen_parent_folder}
	cp ./Makefile ${qiagen_parent_folder}/${SOURCE}

install:
	make conda_install
	#make install_python27_env_if_needed
	make pythonmodules
	make thirdparty_tools
	make data_files
	make help
conda_install: clean
	#Install the Miniconda Python pachages manager
	echo "Next, the Miniconda package will be downloaded and installed"
	echo "You should install it at the default location shown"
	wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
	chmod +x Miniconda2-latest-Linux-x86_64.sh
	sudo bash Miniconda2-latest-Linux-x86_64.sh -p ${conda_home} -u -b
	rm Miniconda2-latest-Linux-x86_64.sh
	${conda_home}/bin/conda init --bash
	#Make the updated shell path available in this session:
	#source ~/.bashrc
	#source ${conda_home}/bin/activate ${conda_env}
	echo Conda was installed in the ${conda_home} folder. The environment that will be used is ${conda_env}.
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

pythonmodules:
	################ Install python modules ################
	## Install some modules with conda
	#This includes R (rstudio) and biopython
	./conda_packages.bash ${conda_home} ${CONDA_PREFIX} ${conda_env}
	./install_perl_modules.bash ${conda_home} ${CONDA_PREFIX} ${conda_env}
	./get_snpeff_data.bash ${conda_home} 4.3.1t-3 4_3 GRCh37.75 #If you wanted to use the GRCh38, you should replace GRCh37.75 to GRCh38.86 in this line

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
	cd py-editdist-0.3 && sudo ${CONDA_PREFIX}/bin/python setup.py install && \
	cd ${qiagen_parent_folder}/bin/downloads/ && tar -xvf sendgrid-v2.2.1.tar.gz && \
	cd sendgrid-python-2.2.1 && sudo ${CONDA_PREFIX}/bin/python setup.py install
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
	wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz \
	https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz \
	https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt  \
	https://github.com/qiaseq/qiaseq-dna/files/2405163/CDHS-13593Z-900.primer3.txt  \
	https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed \
	https://github.com/qiaseq/qiaseq-dna/files/2405163/CDHS-13593Z-900.roi.txt \
	-P ${qiagen_parent_folder}/example/
	mv ${qiagen_parent_folder}/example/CDHS-13593Z-900.roi.txt ${qiagen_parent_folder}/example/CDHS-13593Z-900.roi.bed
	ln /srv/qgen/example/NEB_S2_L001_R1_001.fastq.gz /srv/qgen/example/Cellline_S10_L001_R1_001.fastq.gz
	ln /srv/qgen/example/NEB_S2_L001_R2_001.fastq.gz /srv/qgen/example/Cellline_S10_L001_R2_001.fastq.gz

testfiles:
	## Add test files for smCounterv2
	mkdir -p ${qiagen_parent_folder}/test_smcounter-v2/
	wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed \
	https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam \
	https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
	https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai \
	-P ${qiagen_parent_folder}/test_smcounter-v2/

genomes:
	#Index the genome fasta file, using samtools and bwa, only if does not exists a file with a md5 hash identical to a hash annotated in a file generated after
	#a successful bwa run below
	## Download genome files
	mkdir -p ${qiagen_parent_folder}/data/genome
	#ls ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || (
	md5sum -c ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 &>/dev/null && echo found bwa results file with the expected hash || ( \
	wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
	https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P ${qiagen_parent_folder}/data/genome/ && \
	cd ${qiagen_parent_folder}/data/genome && gunzip ucsc.hg19.fa.gz && \
	${CONDA_PREFIX}/bin/samtools faidx ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa && \
	${CONDA_PREFIX}/bin/bwa index ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa && \
	md5sum -b ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac > ${qiagen_parent_folder}/data/genome/ucsc.hg19.fa.pac.md5 )

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
	@cd ${qiagen_parent_folder}/qiaseq-dna && python run_qiaseq_dna.py --help
	@echo you may run smcounter with commands like that:
	@echo
	@echo 1. Smcounterv1:
	@echo time python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 single out1 NEB_S2 &> run_smcounterv1.log &
	@echo
	@echo 2. Smcounterv2:
	@echo time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2 NEB_S2 &> run_smcounterv2.log &
	@echo
	@echo 3. Tumor-normal analysis:
	@echo python run_qiaseq_dna.py run_sm_counter_v1.params.txt v1 tumor-normal tumor_readset normal_readset > run_smcounterv1_tn.log 2>&1 &
	@echo
	@echo 4. Multiple samples:
	@echo time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2_multisamples_{0} sample1 sample2 sample3 (...) samplen &> run_v6.log &
	@echo For instance:
	@echo time python run_qiaseq_dna.py run_sm_counter_v2.params.txt v2 single out2_multisamples_{0}_{1} NEB_S2 NEB_S2 &> run_smcounterv2_multiplesamples.log &


#Cleaning
procclean:
	pkill python
	pkill conda
clean:
	sudo rm -rf ${conda_home} /srv/qgen/conda ~/.conda ~/.condarc ~/.continuum


