#!/bin/bash
#Code for installing the qiagen-dna software and the example in BASH
#Install the packages needed to start (Note that to get his file you should have installed git earlier, buy the word git stays here for
#   informative purposes: no hurt for re-trying to install it)
sudo yum -y install git unzip cpan wget gcc gcc-c++ bzip2 python-devel nano expat-devel openssl-devel perl perl-CPAN perl-devel curl gcc perl-App-cpanminus
#srv_qiagen=/srv/qgen
qseq_folder_title=qiaseq-dna
parent_folder=/srv/qgen
qseq_folder=${parent_folder}/${qseq_folder_title}
versvar=${0##.*install_qiaseq-dna-v}
version=${versvar%%.bash} #to have this variable running, one have to execute only the versioned verison of this script

parent_folder_permissions=755

function install(){
    /usr/bin/make condaclean
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
    sudo tar -xvzf ${qseq_folder_title}-${version}.tar.gz
    if [ -e ${qseq_folder}-old ]; then sudo rm -rf ${qseq_folder}-old; fi
    if [ -e ${qseq_folder} ]; then mv -Tf ${qseq_folder} ${qseq_folder}-old; fi
    mkdir -p ${parent_folder}
    sudo mv -T ${qseq_folder_title}-${version} ${qseq_folder}
    sudo chown ${USER}:${USER} -R ${qseq_folder}
    cd "${qseq_folder}" && install $@
fi
#Sets up a script with the environment variables needed
#To uninstall:
#sudo rm -rf /srv/qgen; sudo rm -rf /opt/conda; sudo rm -rf /root/conda

################# R packages ################
#echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" | sudo tee /root/.Rprofile >/dev/null

