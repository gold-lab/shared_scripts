
#######################################################################################################################
#!/bin/bash          
#######################################################################################################################


# lncRNA Annotator
# A script to perform genomic positional annotation of lncRNAs.

# Version 5
# April 2014

# By Rory Johnson, University of Bern: rory.johnson@dbmr.unibe.ch or roryjohnson1@gmail.com or www.gold-lab.org

# Usage: sh lncrna.annotator.5.sh project_name lncRNA_annotation_gtf entire_annotation_gtf

####
#VERSION 2: change intergenic lncRNA-protcod distance estimation so now the reported distance is that given by BEDTOOLS, ie the nearest to nearest distance, NOT the promoter-promoter distance as previously.
####
#VERSION 3: change the commands that create the lncRNA and pc BED files (genes, transcripts, exons) to be based on perl pattern matching, so that they are not sensitive to the order that things are defined in the description column of the GTF, which can sometimes change.
####
#VERSION 5: follows from V3 (except fragment from V4, where indicated). Remove the extra V4 functionality but here fixing the problem with unstranded transcripts.
# Add new classification for unstranded trxs.
###


# Important points:

# Input files: Two separate GTF format files - (1) lncRNA annotation and (2) full annotation, both available from gencodegenes.org.
# For example: gencode.v19.long_noncoding_RNAs.gtf and gencode.v19.annotation.gtf.

# LncRNAs are treated at transcript level, protein-coding genes are treated at gene level.

# Output files are saved in a folder called: <<project_name>>.output
# Principal output file: <<project_name>>.output/<<project_name>>.classification.output
# Format: LncRNA_Transcript_ID, Nearest_protein_coding_gene_ID, Positional_annotation_class, Distance (bp) 

# How distances are reported:
# Lncrna-protcod distance is reported differently for intergenic or genic lncRNAs:
# Intergenic-distance of nearest points
# Genic-distance of promoters: negative=lncRNA upstream of protcod, positive= lncRNA downstream of protcod

# LncRNA types: 
#(0) unstranded, intergenic
#(1) samestrand, lincRNA upstream  
#(2) divergent  
#(3) samestrand, protcod upstream  
#(4) convergent  
#(5) intronic_AS  
#(6) intronic_SS  
#(7) exonic_AS  
#(8) exonic_SS
#(9) unstranded, genic
# In cases where lncRNA can have two or more types, the higher-numbered type gets chosen.

# There is no minimum window size around protein coding genes.

# Note: lncRNA annotation must contain transcripts and exons, entire annotation must contain genes. Other features will be discarded.

# Note: if you need to create the lncRNA subannotation gtf file from the total file, use something like this command:
# cat gencode20_havana.gtf | awk '$14!~/protein|pseudo|TR_gene|IG_gene/' > gencode20_havana.lncRNA.gtf
#######################################################################################################################
# Inputs and usage

PROJECT_NAME=$1  # used to name output folder and files
LNCRNA_GTF=$2  # the lncRNA GTF file
TOTAL_GTF=$3   # the full annotation GTF file

[ $# -eq 0 ] && { echo "Usage: $0 Project_name LncRNA_GTF_file Entire_GTF_file"; exit 1; }

#########################################################################################################################################################
# Outputs
#########################################################################################################################################################

# Make output folder
OUTPUT_FOLDER="$PROJECT_NAME.output"
mkdir $OUTPUT_FOLDER

# this file holds an interim, redundant classification
DEFINITION_TEMP="$PROJECT_NAME.classification.temp"

#########################################################################################################################################################
# Create BED files 
#########################################################################################################################################################
# This section taken from V4 initially
# Remove unstranded entries for lncRNA

# Define a BED of lncRNA TSS-TTS boundaries, transcript level
LNCRNA_BOUNDARIES="$PROJECT_NAME.lncrna_boundaries.bed"
cat $LNCRNA_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/transcript_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "transcript"){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' > $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES
#cat $LNCRNA_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/transcript_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "transcript" && $line[6]=~/[+-]/){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' > $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES

#Define a BED of lncRNA exons
LNCRNA_EXONS="$PROJECT_NAME.lncrna_exons.bed"
#cat $LNCRNA_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/transcript_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "exon"){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' | sort -k4 > $OUTPUT_FOLDER/$LNCRNA_EXONS
cat $LNCRNA_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/transcript_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "exon" && $line[6]=~/[+-]/){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' | sort -k4 > $OUTPUT_FOLDER/$LNCRNA_EXONS

# Define a BED of protcod TSS-TTS boundaries, gene level
PROTCOD_BOUNDARIES="$PROJECT_NAME.protcod_boundaries.bed"
cat $TOTAL_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/gene_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "gene" && $line[8]=~/gene_type \"protein_coding\"/){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' > $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES

# Define a BED of protcod exons
PROTCOD_EXONS="$PROJECT_NAME.protcod_exons.bed"
cat $TOTAL_GTF | perl -ne '@line=split "\t"; ($geneid)=$line[8]=~/gene_id "(.*?)"/; if($line[0]!~/#/ && $line[2] eq "exon" && $line[8]=~/gene_type \"protein_coding\"/){$start=$line[3]-1; print "$line[0]\t$start\t$line[4]\t$geneid\t.\t$line[6]\n"}' > $OUTPUT_FOLDER/$PROTCOD_EXONS

#########################################################################################################################################################
# Intersections to create lncRNA categories
#########################################################################################################################################################
# Intergenic LincRNAs - these are lncRNAs that across their full length, do not overlap protcod boundaries.

#create the set of intergenic RNAs by negatively intersecting with the protein coding gene locations
CLASS_LINC_BED_TEMP="$PROJECT_NAME.class_lincRNA.bed.temp"
intersectBed -v -a $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES > $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP

#(Version2)distance is from Bedtools, $13 in the output

# Type1: same strand, lncRNA upstream
CLASS_LINC_BED="$PROJECT_NAME.class_lincRNA.bed"
closestBed -a $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -d | awk '($6=="+" && $12=="+" && $2<$8){dist=$13;print $4"\t"$10"\t1\t"dist} ($6=="-" && $12=="-" && $3>$9){dist=$13;print $4"\t"$10"\t1\t"dist}' > $OUTPUT_FOLDER/$DEFINITION_TEMP

# Type2: divergent
closestBed -a $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -d | awk '($6=="-" && $12=="+" && $2<$8){dist=$13;print $4"\t"$10"\t2\t"dist} ($6=="+" && $12=="-" && $3>$9){dist=$13;print $4"\t"$10"\t2\t"dist}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP

# Type3: same strand, lncRNA downstream
closestBed -a $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -d | awk '($6=="+" && $12=="+" && $2>$8){dist=$13;print $4"\t"$10"\t3\t"dist} ($6=="-" && $12=="-" && $3<$9){dist=$13;print $4"\t"$10"\t3\t"dist}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP

# Type4: convergent
closestBed -a $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -d | awk '($6=="+" && $12=="-" && $2<$8){dist=$13;print $4"\t"$10"\t4\t"dist} ($6=="-" && $12=="+" && $3>$9){dist=$13;print $4"\t"$10"\t4\t"dist}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP

# New in V5: Type0: unstranded, genic
closestBed -a $OUTPUT_FOLDER/$CLASS_LINC_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -d | awk '($6=="."){dist=$13;print $4"\t"$10"\t0\t"dist}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP


#########################################################################################################################################################
# Genic lncRNAs - overlap protcod in some way

# define GenicRNAs as those that overlap protein coding on either strand
CLASS_GENICRNA_BED_TEMP="$PROJECT_NAME.class_genicRNA.bed.temp"
intersectBed -a $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES | sort -k4 > $OUTPUT_FOLDER/$CLASS_GENICRNA_BED_TEMP

# this is the Genic RNA - protein coding overlap file, with the promoter-promoter distance in the last column ** NEGATIVE value indicates lncRNA promoter UPSTREAM of protcod promoter
# New in V5: for unstranded transcripts, distance is measured from centrepoint.
CLASS_GENICRNA_TEMP="$PROJECT_NAME.class_genicRNA.temp"
#intersectBed -a $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -wa -wb | sort -k4 | awk '$6=="+" && $12=="+"{print $0"\t"$2-$8} $6=="-" && $12=="-"{print $0"\t"$9-$3} $6=="+" && $12=="-"{print $0"\t"$9-$2} $6=="-" && $12=="+"{print $0"\t"$3-$8}'  > $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP
intersectBed -a $OUTPUT_FOLDER/$LNCRNA_BOUNDARIES -b $OUTPUT_FOLDER/$PROTCOD_BOUNDARIES -wa -wb | sort -k4 | awk '$6=="+" && $12=="+"{print $0"\t"$2-$8} $6=="-" && $12=="-"{print $0"\t"$9-$3} $6=="+" && $12=="-"{print $0"\t"$9-$2} $6=="-" && $12=="+"{print $0"\t"$3-$8} $6=="." && $12=="+"{print $0"\t"(($2+$3)/2-$8)} $6=="." && $12=="-"{print $0"\t"($9-($2+$3)/2)}' > $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP

# make another copy that will not subsequently be altered, for troubleshooting purposes
CLASS_GENICRNA_UNFILTERED_TEMP="$PROJECT_NAME.class_genicRNA.unfiltered.temp"
cp $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP $OUTPUT_FOLDER/$CLASS_GENICRNA_UNFILTERED_TEMP

# make a BED of Genic RNA exons
CLASS_GENICRNA_EXONS_BED_TEMP="$PROJECT_NAME.class_genicRNA_exons.bed.temp"
join -1 4 -2 4 $OUTPUT_FOLDER/$LNCRNA_EXONS $OUTPUT_FOLDER/$CLASS_GENICRNA_BED_TEMP | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6}' > $OUTPUT_FOLDER/$CLASS_GENICRNA_EXONS_BED_TEMP

# list the exonic_antisense lncRNAs, and their overlapping protcod: TYPE7
CLASS_EXAS_TEMP="$PROJECT_NAME.class_EXAS.temp"
intersectBed -wa -wb -a $OUTPUT_FOLDER/$CLASS_GENICRNA_EXONS_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_EXONS | awk '$6!=$12' | awk '$6=="+"{print $4"\t"$10"\t7"} $6=="-"{print $4"\t"$10"\t7"}' | sort -u > $OUTPUT_FOLDER/$CLASS_EXAS_TEMP

# list the lncRNAs overlapping a protcod gene exonically on same strand: TYPE8
CLASS_EXSS_TEMP="$PROJECT_NAME.class_EXSS.temp"
intersectBed -wa -wb -a $OUTPUT_FOLDER/$CLASS_GENICRNA_EXONS_BED_TEMP -b $OUTPUT_FOLDER/$PROTCOD_EXONS | awk '$6==$12' | awk '$6=="+"{print $4"\t"$10"\t8"} $6=="-"{print $4"\t"$10"\t8"}' | sort -u > $OUTPUT_FOLDER/$CLASS_EXSS_TEMP

#########################################################################################################################################################
# Construct the final definition file, with format: LNCRNA / nearest_PROTCOD / TYPE / PROM-PROM-DIST
#
# Be careful because many overlaps are nonredundant, ie a lncRNA may overlap several protein coding genes.
# Priority of type in this order: ExSS / ExAS / Intronic, and break ties within type by prom-prom distance
#############################
# Extract the EXSS, and the awk serves to select the correct lncRNA-protcod pairs in cases where the lncRNA overlaps multiple protcods. There may be redundant lncRNAs, ie that overlap multiple genes. Resolve these in the second awk command by selecting the closest promoter.
join -1 1 -2 4 $OUTPUT_FOLDER/$CLASS_EXSS_TEMP $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP |  awk '$2==$12{print $1"\t"$2"\t"$3"\t"$15}' | awk '{ if((data[$1])){ if(dist[$1]<0 && $4> dist[$1]){data[$1]=$0; dist[$1]=$4}; if(dist[$1]>=0 && $4<dist[$1]){data[$1]=$0; dist[$1]=$4}     } else {data[$1]=$0; dist[$1]=$4  };     }  END { for (a in data){ print data[a]}}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP

# Remove these from the genic lncRNA file:
JUNK_TEMP="$PROJECT_NAME.junk.temp"
join -1 4 -2 1 -v1 $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP $OUTPUT_FOLDER/$CLASS_EXSS_TEMP | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > $OUTPUT_FOLDER/$JUNK_TEMP
mv $OUTPUT_FOLDER/$JUNK_TEMP $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP

# Next extract the EXAS, and the awk serves to select the correct lncRNA-protcod pairs in cases where the lncRNA overlaps multiple protcods. There may be redundant lncRNAs, ie that overlap multiple genes. Resolve these in the second awk command by selecting the closest promoter.
join -1 1 -2 4 $OUTPUT_FOLDER/$CLASS_EXAS_TEMP $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP |  awk '$2==$12{print $1"\t"$2"\t"$3"\t"$15}' | awk '{ if((data[$1])){ if(dist[$1]<0 && $4> dist[$1]){data[$1]=$0; dist[$1]=$4}; if(dist[$1]>=0 && $4<dist[$1]){data[$1]=$0; dist[$1]=$4}     } else {data[$1]=$0; dist[$1]=$4  };     }  END { for (a in data){ print data[a]}}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP

# Remove these from the genic lncRNA file:
join -1 4 -2 1 -v1 $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP $OUTPUT_FOLDER/$CLASS_EXAS_TEMP | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > $OUTPUT_FOLDER/$JUNK_TEMP
mv $OUTPUT_FOLDER/$JUNK_TEMP $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP

# Now work out the intronic sense and antisense together. To resolve ambiguities, select the (promoter) nearest lncRNA-protcod pairs, then categorise them by whether theyre on the same strand (TYPE6) or opposite strand (TYPE5)
cat $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP | awk '{ if( line[$4]) {if(dist[$4]<0 && $13>dist[$4] && $13<-(dist[$4])){line[$4]=$0; dist[$4]=$13};  if(dist[$4]>0 && $13<dist[$4] && $13>-(dist[$4])){line[$4]=$0; dist[$4]=$13}      } else {  line[$4]=$0; dist[$4]=$13}       } END{for (a in line){print line[a]}}' | awk '$6==$12{print $4"\t"$10"\t6\t"$13} $6!=$12{print $4"\t"$10"\t5\t"$13}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP 

# New for this V5: collect genic unstranded transcripts and classify as 9 
cat $OUTPUT_FOLDER/$CLASS_GENICRNA_TEMP | awk '$6=="." {print $4"\t"$10"\t9\t"$13}' >> $OUTPUT_FOLDER/$DEFINITION_TEMP 

###############################################################################################################################
# Write results to final output file, resolving any redundancies across classes
#
# Finally, filter the temporary definition file, to resolve cases where a lncRNA belongs to two separate classifications. Resolve by prioritising between classifications, with higher numbered types having priority.

DEFINITION_FILE="$PROJECT_NAME.classification.output"

cat  $OUTPUT_FOLDER/$DEFINITION_TEMP | awk '{  if(data[$1]){if($3>class[$1]){ data[$1]=$0; class[$1]=$3;   }   } else { data[$1]=$0; class[$1]=$3;}     } END{for (a in data){print data[a]}} ' | sort > $OUTPUT_FOLDER/$DEFINITION_FILE
















