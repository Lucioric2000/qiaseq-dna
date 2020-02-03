#!/bin/bash
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
sudo yum -y install git unzip cpan wget gcc gcc-c++ bzip2 python-devel nano expat-devel openssl-devel perl perl-CPAN perl-devel curl gcc perl-App-cpanminus
#srv_qiagen=/srv/qgen
qseq_folder_title=qiaseq-dna
parent_folder=/srv/qgen
qseq_folder=${parent_folder}/${qseq_folder_title}
version=15.4
parent_folder_permissions=755

function install(){
    /usr/bin/make clean
    /usr/bin/make install
    source ~/.bashrc
}
qseqdnamatch=`expr match "$(pwd)" '.*\(qiaseq-dna\)'`
echo "Installing ${qseq_folder_title} version ${version}"
if [[ $qseqdnamatch = "qiaseq-dna" ]]
then
    echo "Already in qiaseq-dna folder ($(pwd))."
    #sudo chmod -R ${parent_folder_permissions} ${parent_folder}
    install
else
    p=$(pwd);
    echo "Not in qiaseq_dna folder, but in $p."
    #sudo mkdir -p ${parent_folder}
    #sudo chmod -R ${parent_folder_permissions} ${parent_folder}
    #cd "${parent_folder}" && tar -xvzf ${qseq_folder}-${version}.tar.gz
    tar -xvzf ${qseq_folder_title}-${version}.tar.gz
    if [ -e ${qseq_folder}-old ]; then sudo rm -rf ${qseq_folder}-old; fi
    if [ -e ${qseq_folder} ]; then sudo mv -Tf ${qseq_folder} ${qseq_folder}-old; fi
    if [ -e ${qseq_folder_title} ]; then sudo rm -rf ${qseq_folder_title}; fi
    #sudo mv ${qseq_folder_title}-${version} ${qseq_folder_title}
    #sudo mv -T ${qseq_folder_title} ${qseq_folder}
    sudo mv -T ${qseq_folder_title}-${version} ${qseq_folder}
    cd "${qseq_folder}" && install $@
fi
#Sets up a script with the environment variables needed
#To uninstall:
#sudo rm -rf /srv/qgen; sudo rm -rf /opt/conda; sudo rm -rf /root/conda

################# R packages ################
#echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" | sudo tee /root/.Rprofile >/dev/null

