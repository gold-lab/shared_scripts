
# lncRNA Annotator
## A script to perform genomic positional annotation of lncRNAs

##### Version 5;  April 2014
##### By Rory Johnson, University of Bern: rory.johnson@dbmr.unibe.ch or roryjohnson1@gmail.com or www.gold-lab.org
#
### Usage: sh lncrna.annotator.5.sh project_name lncRNA_annotation_gtf entire_annotation_gtf
###
#### Input files: Two separate GTF format files - (1) lncRNA annotation and (2) full annotation, both available from gencodegenes.org.
##### For example: gencode.v19.long_noncoding_RNAs.gtf and gencode.v19.annotation.gtf.
##### GTF input files need to be previouly sorted
###
#### Output files are saved in a folder called: <<project_name>>.output
##### Principal output file: <<project_name>>.output/<<project_name>>.classification.output
##### Format: LncRNA_Transcript_ID, Nearest_protein_coding_gene_ID, Positional_annotation_class, Distance (bp) 
#
###
### Important points:

* Installation of BEDTOOLS is required for executing the script
* LncRNAs are treated at transcript level, protein-coding genes are treated at gene level.
* Lncrna-protcod distance is reported differently for intergenic or genic lncRNAs:
#####    - Intergenic-distance of nearest points
#####    - Genic-distance of promoters: negative=lncRNA upstream of protcod, positive= lncRNA downstream of protcod
* LncRNA types: 

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
#### In cases where lncRNA can have two or more types, the higher-numbered type gets chosen.
#### There is no minimum window size around protein coding genes.
#
#### Note: lncRNA annotation must contain transcripts and exons, entire annotation must contain genes. Other features will be discarded.
#### Note: if you need to create the lncRNA subannotation gtf file from the total file, use something like this command:
##### cat gencode20_havana.gtf | awk '$14!~/protein|pseudo|TR_gene|IG_gene/' > gencode20_havana.lncRNA.gtf
