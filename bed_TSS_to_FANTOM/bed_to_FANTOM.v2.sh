#!/bin/bash

######### BED coordinates to FANTOM peaks ###############
#							#
# Uses input bed file with TSS coordinates to check for	#
# closest and highest FANTOM peak. 			#
# If no peak is found then original coordinates in	#
# input bed are kept.					#
# Path to FANTOM faile is required.			#
#							#
# # Dependencies					#
#    -bedtools 2.26.0					#
#    -python						#
#							#
# # How to use:						#
#    ./bed_TSS_to_FANTOM.sh -i test_files/input.bed \	#
#    -f test_files/mm9_to_mm10_liftover.bed		#
#							#
#########################################################

## v2
# RANGE modified.
#   Range is the +/- distance from TSS to look for FANTOM peak
#
# tmp/increased.bed file not used anymore.
#   distance is directly checked in the closestBed step 


# Add bedtools module
module add UHTS/Analysis/BEDTools/2.26.0

# Capture arguments from command line
DIST_UP=1000
DIST_DOWN=1000

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -i|--input-bed)
    INFILE="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--fantom)
    FANTOM="$2"
    shift # past argument
    shift # past value
    ;;
    -du|--dist_up_stream)
    RANGE="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

INFILE="${INFILE}"
FANTOM="${FANTOM}"
RANGE="${RANGE}"


# Check mandatory arguments
if [ -z $INFILE ]; then
	echo "-i and -f parameters must be specified"; exit 1;
fi
if [ -z $FANTOM ]; then
	echo "-i and -f parameters must be specified"; exit 1;
fi

#create tmp dir
mkdir -p tmp

# Change to FANTOM TSS when distance < 1000nt
# In cases of tie in clossenes to FANTOM peaks, highest scored peak is selected.
range=$RANGE
sortBed -i $INFILE > tmp/sorted.bed
closestBed -s -d -k 100 -a tmp/sorted.bed -b $FANTOM | awk -v x=$range '{ if ($NF<=x) {OFS="\t"; print $1,$13,$14,$4,$11,$6} else {OFS="\t"; print $1,$2,$3,$4,-1,$6}}' > tmp/intersect.bed
python /data/projects/p283_rna_and_disease/scripts/bed_TSS_to_FANTOM/selectFANTOMpeak.py tmp/intersect.bed

rm -R tmp
