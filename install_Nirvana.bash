#!/bin/bash

# adjust these paths to reflect where you have downloaded the Nirvana data files
# In this example, we assume that the Cache, References, and SupplementaryDatabase
# folders have been downloaded into the NIRVANA_ROOT folder.

# In addition to downloading the Nirvana data files, make sure you have .NET Core 2.0
# installed on your computer:
# https://www.microsoft.com/net/download/core

#. NET Core install:

sudo rpm -Uvh https://packages.microsoft.com/config/rhel/7/packages-microsoft-prod.rpm
sudo yum update
sudo yum install aspnetcore-runtime-2.1 dotnet-sdk-2.1

# Nirvana install

TOP_DIR=~
NIRVANA_ROOT=$TOP_DIR/Nirvana
NIRVANA_BIN=$NIRVANA_ROOT/bin/Release/netcoreapp2.1/Nirvana.dll
DATA_DIR=$NIRVANA_ROOT/Data
DOWNLOADER_BIN=$NIRVANA_ROOT/bin/Release/netcoreapp2.1/Downloader.dll
NIRVANA_TAG=v3.5.0

# just change this to GRCh38 if you want to set everything up for hg38
GENOME_ASSEMBLY=GRCh37
#SA_VERSION=44
#CACHE_VERSION=26
#REF_VERSION=5

#CACHE_DIR=$DATA_DIR/Cache/$CACHE_VERSION/$GENOME_ASSEMBLY
#SA_DIR=$DATA_DIR/SupplementaryDatabase/$SA_VERSION
#REF_DIR=$DATA_DIR/References/$REF_VERSION

#CACHE_TGZ=$DATA_DIR/v${CACHE_VERSION}.tar.gz
#SA_TGZ=$DATA_DIR/v${SA_VERSION}_${GENOME_ASSEMBLY}.tar.gz
#REF_TGZ=$DATA_DIR/v${REF_VERSION}.tar.gz

#CACHE_TEST=$CACHE_DIR/Ensembl.transcripts.ndb
#SA_TEST=$SA_DIR/$GENOME_ASSEMBLY/chr1.nsa
#REF_TEST=$REF_DIR/Homo_sapiens.${GENOME_ASSEMBLY}.Nirvana.dat
#SOURCE_TEST=$NIRVANA_ROOT/Nirvana.sln

SA_DIR=$DATA_DIR/SupplementaryAnnotation/$GENOME_ASSEMBLY
REF_DIR=$DATA_DIR/References
CACHE_DIR=$DATA_DIR/Cache/$GENOME_ASSEMBLY
REF_TEST=$REF_DIR/Homo_sapiens.${GENOME_ASSEMBLY}.Nirvana.dat


# =====================================================================

YELLOW='\033[1;33m'
RESET='\033[0m'

echo -ne $YELLOW
echo " _   _ _                             "
echo "| \ | (_)                            "
echo "|  \| |_ _ ____   ____ _ _ __   __ _ "
echo "| . \` | | '__\ \ / / _\` | '_ \ / _\` |"
echo "| |\  | | |   \ V / (_| | | | | (_| |"
echo "|_| \_|_|_|    \_/ \__,_|_| |_|\__,_|"
echo -e $RESET

# create the data directories
create_dir() {
    if [ ! -d $1 ]
    then
        mkdir -p $1
    fi
}
sudo_create_dir() {
    if [ ! -d $1 ]
    then
        sudo mkdir -p $1
    fi
}
download() {
    FILENAME=$1
    TOP_URL=$2
    LOCAL_PATH=$SA_DIR/$FILENAME

    if [ ! -f  $LOCAL_PATH ]
    then
    echo "- downloading $FILENAME"
    pushd $SA_DIR
    curl -s -O http://d24njugtyo9gfs.cloudfront.net/$TOP_URL/$GENOME_ASSEMBLY/$FILENAME
    popd
    fi
}


# =====================
# unpack the data files
# =====================

unpack_file() {
    if [ ! -f $4 ]
    then
	pushd $2 > /dev/null
	echo -n "- unpacking $1 files... ($2) ($3) ($4)"
	sudo tar -xfz $3
	echo "finished."
	popd > /dev/null
	rm $3
    fi
}

# ==================================
# clone the Nirvana repo (if needed)
# ==================================

pushd $TOP_DIR

if [ ! -d $NIRVANA_ROOT ]
then
    git clone https://github.com/Illumina/Nirvana.git
    git checkout $NIRVANA_TAG
fi


# =============
# build Nirvana
# =============

pushd $NIRVANA_ROOT
if [ ! -f $NIRVANA_BIN ]
then
    #pushd Nirvana > /dev/null
    dotnet clean
    dotnet build -c Release -v diagnostic
    #popd > /dev/null
fi

# download the Nirvana source
if [ ! -f $REF_TEST ]
then
    create_dir $DATA_DIR
    echo dotnet $DOWNLOADER_BIN --ga $GENOME_ASSEMBLY --out $DATA_DIR
    dotnet $DOWNLOADER_BIN --ga $GENOME_ASSEMBLY --out $DATA_DIR
fi


# ==============================
# run Nirvana on a test VCF file
# ==============================

# download a test vcf file
if [ ! -f HiSeq.10000.vcf ]
then
    wget https://github.com/samtools/htsjdk/raw/master/src/test/resources/htsjdk/variant/HiSeq.10000.vcf
fi

# analyze it with Nirvana
COMMAND="dotnet $NIRVANA_BIN -c $CACHE_DIR/Both --sd $SA_DIR -r $REF_TEST -i HiSeq.10000.vcf -o HiSeq.10000.annotated"
echo Running $COMMAND
$COMMAND

popd
popd