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
#    -bedtools						#
#    -python						#
# 							#
# # How to use:						#
#    ./bed_TSS_to_FANTOM.sh -i test_files/input.bed \	#
#    -f test_files/mm9_to_mm10_liftover.bed		#
#							#
#########################################################

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
    DIST_UP="$2"
    shift # past argument
    shift # past value
    ;;
    -dd|--dist_down_stream)
    DIST_DOWN="$2"
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
DIST_UP="${DIST_UP}"
DIST_DOWN="${DIST_DOWN}"

# Check mandatory arguments
if [ -z $INFILE ]; then
	echo "-i and -f parameters must be specified"; exit 1;
fi
if [ -z $FANTOM ]; then
	echo "-i and -f parameters must be specified"; exit 1;
fi

#create tmp dir
mkdir -p tmp

# Increase TSS's
plus=$DIST_UP
minus=$DIST_DOWN
awk -v x=$plus -v y=$minus '{OFS="\t"; print $1,$2-x,$3+y,$4,$5,$6}' $INFILE > tmp/increased.bed
sortBed -i tmp/increased.bed > tmp/tmp.bed
mv tmp/tmp.bed tmp/increased.bed

# Change to FANTOM TSS when distance < 1000nt
# In cases of tie in clossenes to FANTOM peaks choose manually the one with higer score (column 6, third digit)
closestBed -s -d -a tmp/increased.bed -b $FANTOM | awk '{ if ($NF<=0) {OFS="\t"; print $1,$13,$14,$4,$11,$6} else {OFS="\t"; print $1,$2,$3,$4,$5,$6}}' > tmp/intersect.bed
python selectFANTOMpeak.py tmp/intersect.bed

rm -R tmp


# Extract unique candidate gene ids
#cut -f7 Data/lincRNA_candidates.tsv | sort | uniq > Data/lincRNA_gene_ids.txt
#cut -f6 Data/epifactor.tsv | sort | uniq > Data/epifactor_gene_ids.txt
# Extract annotation for selected gene IDs
#grep -f Data/lincRNA_gene_ids.txt Data/mm10_gencode.vM16.basic.annotation.gtf > input/lincRNA.gtf
#grep -f Data/epifactor_gene_ids.txt Data/mm10_gencode.vM16.basic.annotation.gtf > input/epifactor.gtf
# Create BED with TSS's from gtf's
#cat input/lincRNA.gtf | awk '$3=="gene"' | awk 'BEGIN {OFS = "\t"} {if($7=="+"){print $1,$4-1, $4, $10, 1, "+"} else {print $1,$5-1, $5, $10, 1, "-"}}' > input/lincRNA.bed
#cat input/epifactor.gtf | awk '$3=="gene"' | awk 'BEGIN {OFS = "\t"} {if($7=="+"){print $1,$4-1, $4, $10, 1, "+"} else {print $1,$5-1, $5, $10, 1, "-"}}' > input/epifactor.bed



# Extract genomic sequences upstream of TSS
#./extract_candidates.sh

# Extract sgRNAs for each genomic region
#python extract_seqs.py

# Obtain off-target information for each sgRNA
# Run next command in my laptop; MySQL off-targets database is already prepared there
#python off_targets.py out/epifactor_to_score.out out/epifactor_w_offtargets.out
#python off_targets.py out/lincRNA_to_score.out out/lincRNA_w_offtargets.out

# Score each sgRNA with Rule_Set_2
#python score.py out/epifactor_w_offtargets.out out/epifactor_avaluated.out
#python score.py out/lincRNA_w_offtargets.out out/lincRNA_avaluated.out

# Sort results
#python sortDoenchCRISPRa.py out/epifactor_avaluated.out input/epifactor_FANTOM_increased.bed out/epifactor_top10_selected.out
#python sortDoenchCRISPRa.py out/lincRNA_avaluated.out input/lincRNA_FANTOM_increased.bed out/lincRNA_top10_selected.out

# Result to BED format
#awk 'NR>1{print "chr"$2, $3, $4, $1, $12, $6}' out/lincRNA_top10_selected.out
#awk 'NR>1{print "chr"$2, $3, $4, $1, $12, $6}' out/epifactor_top10_selected.out



