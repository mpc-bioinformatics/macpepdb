#!/bin/bash

if [ "$1" == "--help" ]; then
    cat << EndOfMessage
This script check if all dependencies for a database update are available, downloads the last UniProt release and starts an update.
Please make sure you created the database once before, so the digestion parameter are saved in the database.
Parameter:
    -t Number of threads
    -d Database URL (postgresql://user:password@ip_or_name:port/database)
EndOfMessage
    exit 0
fi

while getopts ":t:d:" option
do
    case "$option"
    in
        t)
            THREADS=${OPTARG}
            ;;
        d)
            DATABASE_URL=${OPTARG}
            ;;
        ?)
            echo "${option} is missing."
            exit
            ;;
    esac
done

python -m macpepdb digestion -i foo.txt -f text -l ./digestion_test.log -u ./digestion_test.statistics.csv -s ./digestion_test.statistics.csv -e trypsin -c 2 -t 1 -d $DATABASE_URL --minimum-peptide-length 10 --maximum-peptide-length 20

DEPENDENCY_IS_MISSING=0

check_if_software_is_available() {
    which $1 > /dev/null 2>&1
    exit_code=$?
    if [ "$exit_code" != "0" ]
    then
        echo "Can't find '${1}' in your path."
        DEPENDENCY_IS_MISSING=1
    fi
}

# Check dependencies
## curl
check_if_software_is_available "curl"
## gunzip
check_if_software_is_available "gunzip"
## gunzip
check_if_software_is_available "grep"
## python
check_if_software_is_available "python"

# Check if MaCPep Python module is installed
python -m macpepdb --help > /dev/null 2>&1
exit_code=$?
if [ "$exit_code" != "0" ]
then
    echo "MaCPepDB is not installed."
    DEPENDENCY_IS_MISSING=1
fi

if [ "$DEPENDENCY_IS_MISSING" != "0" ]
then
    echo "There are missing dependencies. Make sure all necessary dependencies are installed. Skipping UniProt download and database update."
    exit 1
fi

echo "Removing old data..."
rm -rf ./logs
rm -f ./uniprot_sprot.dat
rm -f ./uniprot_sprot.dat

echo "Download current release of Swiss-Prot and TrEMBL from ftp.expasy.org ..."
curl -O ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
curl -O ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz 
echo "Uncompress both files. This may take a while for TrEMBL ..."
gunzip uniprot_sprot.dat.gz 
gunzip uniprot_trembl.dat.gz

mkdir ./logs

# The digestion paramter are saved in the database from the first manually digestion, they will override the ones set by the CLI so no need to set the correct parameter here. 
python -m macpepdb digestion -i ./uniprot_sprot.dat -i ./uniprot_trembl.dat -f text -l ./logs/uniprot_0.log -u ./logs/uniprot_unprocessible_0.txt -s ./logs/uniprot_statistics_0.csv -e trypsin -c 1  -t $THREADS -d $DATABASE_URL
echo "Check if there were unprocessible proteins ..."
NUMBER_OF_UNPROCESSIBLE_PROTEINS=$(grep -c ">" ./logs/uniprot_unprocessible_0.txt)
if [ "$NUMBER_OF_UNPROCESSIBLE_PROTEINS" == "0" ]
then
    echo "No unprocessible proteins found."
else
    echo "Found unprocessible proteins, digest them a single thread only."
    python -m macpepdb digestion -i ./logs/uniprot_unprocessible_0.txt -f fasta -l ./logs/uniprot_1.log -u ./logs/uniprot_unprocessible_1.txt -s ./logs/uniprot_statistics_1.csv -e trypsin -c 1 -t 1 -d $DATABASE_URL
    NUMBER_OF_UNPROCESSIBLE_PROTEINS=$(grep -c ">" ./logs/uniprot_unprocessible_1.txt)
    if [ "$NUMBER_OF_UNPROCESSIBLE_PROTEINS" == "0" ]
    then
        echo "All former unprocessible proteins were processed."
    else
        echo "Found unprocessible protein again. Please check logs in: ./logs/uniprot_1.log"
    fi
fi
