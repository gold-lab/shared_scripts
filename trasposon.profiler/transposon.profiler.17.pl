
# Transposon Profiler

# Version 17

# A program to annotate Transposable and Repeat Elements in LncRNA

# Inputs are a GTF annotation and a RepeatMasker file

# Rory Johnson 




# Typical command: 
# perl transposon.profiler.16.pl -e input.files/gencode.v18.long_noncoding_RNAs.gtf -r input.files/repeat.masker.hg19.test -o test17.out

# Explanation of overlap classes:

#+1	+2	+3	+4	+5	+6
#"tss.plus", "splice_acc.plus","encompass.plus","splice_don.plus","tts.plus", "inside.plus"

##############################################################
#!/usr/bin/perl -w

use strict;

##############################
# Timing, see end
BEGIN { our $start_run = time(); }
###############################

my $output_folder;
my $repeat_masker_file;
my $lncrna_gtf_file;

# pass the program: the lncRNA GTF file, the repeatmasker annotation file, prefix for output folder.

my $USAGE = "\nUsage: $0 -e lncRNA_gtf_file -r repeatmasker_file -o output_prefix \n\n Important notes:\n(1) The gtf file must have genes, transcripts and exons \n\n";

unless (@ARGV) {
	print $USAGE;
	exit;	
}

for (my $i = 0; $i < scalar(@ARGV); ++$i) {
	if ( $ARGV[$i] eq '-e' ){
	$lncrna_gtf_file = $ARGV[$i + 1]; 
	}
	elsif ( $ARGV[$i] eq '-r' ) {
	$repeat_masker_file = $ARGV[$i + 1]; 
	}
	elsif ( $ARGV[$i] eq '-o' ) {
	$output_folder = $ARGV[$i + 1]; 
	}
}

#################################################################################################################################################
# Create output folder

$output_folder=$output_folder.".transposon_profiler.out";

system("mkdir $output_folder");

##############################################################################################################################################
# Important variables
my $min_bed_length=10;             #the minimum length in nt that a single overlap instance must have to be considered
my $min_overlap_count=20;              #minimum number of total overlap instances that a particular repeat must have to be considered
my $min_overlap_nt=500;               #this is the minimum total exonic overlap nt that a particular repeat type must have to be considered

##############################################################################################################################################
# Output files

# Merged exons of the lncRNA transcripts
my $exon_merged_bed=$output_folder.'/exon_merged.bed';
# Unmerged exons
my $exon_unmerged_bed=$output_folder.'/exon_unmerged.bed';
# Entire span of genes 
my $gene_bed_file=$output_folder.'/gene.bed';
#
my $intersect_exon_repeat_antisense_bed=$output_folder.'/intersect_bed_repeat_lncrna_exon.antisense.bed';
#
my $intersect_exon_repeat_bed=$output_folder.'/intersect_bed_repeat_lncrna_exon.bed';
#
my $intersect_exon_repeat_samesense_bed=$output_folder.'/intersect_bed_repeat_lncrna_exon.samesense.bed';
#
my $intersect_intron_repeat_antisense_bed=$output_folder.'/intersect_bed_repeat_lncrna_intron.antisense.bed';
#
my $intersect_intron_repeat_bed=$output_folder.'/intersect_bed_repeat_lncrna_intron.bed';
#
my $intersect_intron_repeat_samesense_bed=$output_folder.'/intersect_bed_repeat_lncrna_intron.samesense.bed';
# a BED of the entire transcript span
my $transcript_bed_file=$output_folder.'/transcript.bed';
# The entire RepBase as a BED file
my $adjusted_repeat_filename=$output_folder.'/repeat_masker_adjusted.bed';
# A BED of unambiguously stranded (NOT all) lncRNA introns
my $lncrna_intron_file=$output_folder.'/intron_stranded.bed';

#filename for the exon GTF
my $gtf_out=$output_folder.'/repeat_lncrna_exon_intersect.gtf';    
#filename for the intron GTF
my $gtf_out_intron=$output_folder.'/repeat_lncrna_intron_intersect.gtf';     
#these are the parts of repeats that intersect exons, but the actual parts that do not overlap
my $trimmings_bed_file=$output_folder.'/overlap_trimmings.bed';   
#overlap statistics for each repeat type
my $stats_filename=$output_folder.'/overlap_statistics.txt';
#summary statistics
my $summary_stats_filename=$output_folder.'/summary_statistics.txt';
#
my $te_overlap_type_table=$output_folder.'/te_overlap_type_table.txt';

my $profile_exon_file=$output_folder.'/repeat_exon.profile';
my $profile_intron_file=$output_folder.'/repeat_intron.profile';
my $profile_exon_ss_file=$output_folder.'/repeat_exon_ss.profile';
my $profile_intron_ss_file=$output_folder.'/repeat_intron_ss.profile';
my $profile_exon_as_file=$output_folder.'/repeat_exon_as.profile';
my $profile_intron_as_file=$output_folder.'/repeat_intron_as.profile';
my $candidate_list=$output_folder.'/repeat_candidate.list';

my $profile_temp=$output_folder.'/profile_temp';

my $repeat_class_file=$output_folder.'/repeat_class.tmp';


###################################################################################################

# Do separate analysis specifically for TEs that do not overlap the boundary of any annotated exon
# 
# This is to prevent confounding results eg insertion profiles, that arise from TEs containing promoters / TTS / splice sites.
# These no junction files are designated by "no_jn"

my $exon_junction_bed=$output_folder.'/exon_junction.bed';
my $adjusted_repeat_filename_no_jn=$output_folder.'/repeat_masker_adjusted.no_jn.bed';



my $intersect_exon_repeat_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_exon.no_jn.bed';
my $intersect_intron_repeat_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_intron.no_jn.bed';
my $gtf_out_exon_no_jn=$output_folder.'/repeat_lncrna_exon_intersect.no_jn.gtf';    #filename for the exon GTF
my $gtf_out_intron_no_jn=$output_folder.'/repeat_lncrna_intron_intersect.no_jn.gtf';     #filename for the intron GTF
my $intersect_exon_repeat_samesense_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_exon.samesense.no_jn.bed';
my $intersect_exon_repeat_antisense_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_exon.antisense.no_jn.bed';
my $intersect_intron_repeat_samesense_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_intron.samesense.no_jn.bed';
my $intersect_intron_repeat_antisense_bed_no_jn=$output_folder.'/intersect_bed_repeat_lncrna_intron.antisense.no_jn.bed';

##############################################################################################################################################


#Create a BED file holding all the RepeatMasker repeats

# *** check this

print "\n Converting RepeatMasker file\n";

my $rep_coord_start;
my $rep_coord_end;
my $repstrand;
my $repclass;
my $repstart;
my $repend;
my $repleft;

my %rep_class;   #this hash holds key-value repeat_name-repeat_class info
my %rep_family;  #similarly for repeat family

open (REPOUT, ">$adjusted_repeat_filename");
open (REPFILE, $repeat_masker_file);
while (<REPFILE>) {
 	chomp;
	my @line=split" ";
	$rep_class{$line[10]}=$line[12];
	$rep_family{$line[10]}=$line[11];

		if($line[9] eq '+'){
		my $hypothetical_start=$line[6]-$line[13]+1;
		my $hypothetical_end=$line[7]-$line[15];
		print REPOUT "$line[5]\t$line[6]\t$line[7]\t$line[5],$line[6],$line[7],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[15],$hypothetical_start,$hypothetical_end\t1\t$line[9]\n";
		}

		if($line[9] eq '-'){
		my $hypothetical_start=$line[6]+$line[13];
		my $hypothetical_end=$line[7]+$line[15]-1;
		print REPOUT "$line[5]\t$line[6]\t$line[7]\t$line[5],$line[6],$line[7],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[15],$hypothetical_start,$hypothetical_end\t1\t$line[9]\n";
		}
}
close (REPFILE); 
close (REPOUT);

# sort the Bed file
system("cat $adjusted_repeat_filename | sortBed > temp.bed");
system("mv temp.bed $adjusted_repeat_filename");

#############################################################################################################################################
# Create a merged exon BED out of the lncRNA GTF file.

print "\n Creating exon BED file\n";

exonGTFtoBED($lncrna_gtf_file,$exon_unmerged_bed);
my @results=`mergeBed -s -nms -i $exon_unmerged_bed`;
open (OUT, ">$exon_merged_bed"); # this code is simply because mergeBed seems to have a bug and outputs a BED lacking the score column
foreach my $line(@results){
	chomp $line;	
	my @line=split("\t",$line);
	print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t1\t$line[4]\n";
}
close(OUT);

transcriptGTFtoBED($lncrna_gtf_file,$transcript_bed_file);   # make a transcript bed file  ***need to check
geneGTFtoBED($lncrna_gtf_file,$gene_bed_file); #and a gene bed file  ***need to check 

#########################################################################################################################################
#create a filtered repeatmasker file that contains no repeats that overlap exon boundaries

#first, make BED file of all exon boundaries, using UNmerged exon bed

print "\n Creating non-junction RepeatMasker file\n";

system("cat $exon_unmerged_bed | perl -ne '\@line= split \" \"; \$a=\$line[1]-1; \$b=\$line[2]-1; print\"\$line[0]\t\$a\t\$line[1]\t\$line[3]a\t\$line[4]\t\$line[5]\n\$line[0]\t\$b\t\$line[2]\t\$line[3]b\t\$line[4]\t\$line[5]\n\";' > $exon_junction_bed");

#now filter
system("intersectBed -v -a $adjusted_repeat_filename -b $exon_junction_bed > $adjusted_repeat_filename_no_jn");



##############################################################################################################################################
# Create a stranded intron file, being the subtraction of the merged exons from the gene.bed file.

print "\n Creating stranded intron BED file\n";

my $intron_temp_bed=$output_folder.'/intron.temp.bed';  # all the introns
my $intron_exclude_temp_bed=$output_folder.'/intron.exclude.temp.bed'; # introns that overlap another intron on the opposite strand, and therefore do not have a unamibiguous strandedness

system("subtractBed -a $gene_bed_file -b $exon_merged_bed > $intron_temp_bed");
system("intersectBed -a $intron_temp_bed -b $intron_temp_bed -wb | awk '\$4!=\$10 {print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $intron_exclude_temp_bed"); 
system("intersectBed -v -a $intron_temp_bed -b $intron_exclude_temp_bed | sortBed > $lncrna_intron_file");  # remove ambiguous introns
##############################################################################################################################################
# intersections, imposing a minimum length of overlap to be considered, defined in $min_bed_length

print "\n Intersecting repeats and exons\n";

# The exon intersection
repIntersectBED($adjusted_repeat_filename,$exon_merged_bed,$transcript_bed_file,$min_bed_length,$intersect_exon_repeat_bed,$gtf_out);
repIntersectBED($adjusted_repeat_filename_no_jn,$exon_merged_bed,$transcript_bed_file,$min_bed_length,$intersect_exon_repeat_bed_no_jn,$gtf_out_exon_no_jn);   # the no junction version

splitBED($intersect_exon_repeat_bed,$intersect_exon_repeat_samesense_bed,$intersect_exon_repeat_antisense_bed);
splitBED($intersect_exon_repeat_bed_no_jn,$intersect_exon_repeat_samesense_bed_no_jn,$intersect_exon_repeat_antisense_bed_no_jn);    # the no junction version

# The intron intersection
repIntersectBED($adjusted_repeat_filename,$lncrna_intron_file,$transcript_bed_file,$min_bed_length,$intersect_intron_repeat_bed, $gtf_out_intron);
repIntersectBED($adjusted_repeat_filename_no_jn,$lncrna_intron_file,$transcript_bed_file,$min_bed_length,$intersect_intron_repeat_bed_no_jn, $gtf_out_intron_no_jn);    # the no junction version

splitBED($intersect_intron_repeat_bed,$intersect_intron_repeat_samesense_bed,$intersect_intron_repeat_antisense_bed);
splitBED($intersect_intron_repeat_bed_no_jn,$intersect_intron_repeat_samesense_bed_no_jn,$intersect_intron_repeat_antisense_bed_no_jn);    # the no junction version

# The exon intersection NONoverlap - taking those repeats that overlap exons, but then keeping only the fragments of the repeat that do NOT overlap the repeat. AKA TRIMMINGs
system("intersectBed -a $adjusted_repeat_filename -b $intersect_exon_repeat_bed -wa | subtractBed -a stdin -b $intersect_exon_repeat_bed | awk '\$3-\$2>=$min_bed_length' > $trimmings_bed_file");


# Make a count table of TE class vs overlap type
system("cat $intersect_exon_repeat_bed | awk '{print \$4,\$5}' | tr \",\" \"\t\" | awk '{print \$13,\$7}' > $repeat_class_file");

makeTable($repeat_class_file, $te_overlap_type_table);

##############################################################################################################################################
# create a GTF file containing all the exon-contained repeat elements

#print "\n Creating GTF\n";

#repeatGTF($adjusted_repeat_filename, $exon_merged_bed, $gtf_out);


##############################################################################################################################################
# calculate the correlation values

#my %correl_values=calcRcorrel($profile_exon_file, $profile_intron_file);


##############################################################################################################################################
# calculate statistics and print them to a stats file

print "\n Calculating statistics\n";

my %exon_overlap_nt={};
my %exon_overlap_count={};
my $exon_total_length=0;
my %intron_overlap_nt={};
my %intron_overlap_count={};
my $intron_total_length=0;

my %genome_repeat_nt={};
my %genome_repeat_count={};

my %exon_ss_overlap_nt={};
my %exon_ss_overlap_count={};
my %exon_as_overlap_nt={};
my %exon_as_overlap_count={};

my %intron_ss_overlap_nt={};
my %intron_ss_overlap_count={};
my %intron_as_overlap_nt={};
my %intron_as_overlap_count={};

my $exon_binomial_p;
my $intron_binomial_p;
my $exon_binomial2_p;
my $intron_binomial2_p;

%exon_overlap_nt=bedOLnt($intersect_exon_repeat_bed);
%exon_overlap_count=bedOLcount($intersect_exon_repeat_bed);
%intron_overlap_nt=bedOLnt($intersect_intron_repeat_bed);
%intron_overlap_count=bedOLcount($intersect_intron_repeat_bed);

%genome_repeat_nt=bedOLnt($adjusted_repeat_filename);
%genome_repeat_count=bedOLcount($adjusted_repeat_filename);

%exon_ss_overlap_nt=bedOLnt($intersect_exon_repeat_samesense_bed);
%exon_ss_overlap_count=bedOLcount($intersect_exon_repeat_samesense_bed);
%intron_ss_overlap_nt=bedOLnt($intersect_intron_repeat_samesense_bed);
%intron_ss_overlap_count=bedOLcount($intersect_intron_repeat_samesense_bed);

%exon_as_overlap_nt=bedOLnt($intersect_exon_repeat_antisense_bed);
%exon_as_overlap_count=bedOLcount($intersect_exon_repeat_antisense_bed);
%intron_as_overlap_nt=bedOLnt($intersect_intron_repeat_antisense_bed);
%intron_as_overlap_count=bedOLcount($intersect_intron_repeat_antisense_bed);


#counting the total number of nt in the input bed files
open (FILE, $exon_merged_bed);
while (<FILE>) {
	chomp;
	my @line=split" ";
	$exon_total_length+=($line[2]-$line[1]);
}
close (FILE);

open (FILE, $lncrna_intron_file);
while (<FILE>) {
	chomp;
	my @line=split" ";
	$intron_total_length+=($line[2]-$line[1]);
}
close (FILE);

#print this to file, starting with a header:
open (OUT, ">$stats_filename");
open (OUTCAND, ">$candidate_list");
print OUT "ID\tClass\tFamily\tintron_overlap_nt\tintron_overlap_count\tintron_ss_overlap_nt\tintron_ss_overlap_count\tintron_as_overlap_nt\tintron_as_overlap_count\tintron_total_length\texon_overlap_nt\texon_overlap_count\texon_ss_overlap_nt\texon_ss_overlap_count\texon_as_overlap_nt\texon_as_overlap_count\texon_total_length\tgenome_overlap_length\tgenome_overlap_count\n";

foreach my $key ( sort (keys %exon_overlap_nt )){
	#these commands fill in any missing values in the hashes, otherwise the output file has gaps in it in the small number of classes having no overlaps
	unless($exon_ss_overlap_nt{$key}){$exon_ss_overlap_nt{$key}=0};
	unless($exon_as_overlap_nt{$key}){$exon_as_overlap_nt{$key}=0};
	unless($intron_ss_overlap_nt{$key}){$intron_ss_overlap_nt{$key}=0};
	unless($intron_as_overlap_nt{$key}){$intron_as_overlap_nt{$key}=0};
	unless($exon_ss_overlap_count{$key}){$exon_ss_overlap_count{$key}=0};
	unless($exon_as_overlap_count{$key}){$exon_as_overlap_count{$key}=0};
	unless($intron_ss_overlap_count{$key}){$intron_ss_overlap_count{$key}=0};
	unless($intron_as_overlap_count{$key}){$intron_as_overlap_count{$key}=0};
	unless($genome_repeat_nt{$key}){$genome_repeat_nt{$key}=0};
	unless($genome_repeat_count{$key}){$genome_repeat_count{$key}=0};


	if(($exon_overlap_count{$key}>=$min_overlap_count) && ($exon_overlap_nt{$key}>=$min_overlap_nt) && $intron_overlap_nt{$key} && $intron_overlap_count{$key}){	#require minimum exon overlaps, specified at start of script
#	$exon_binomial_p=binomial($exon_ss_overlap_count{$key},$exon_overlap_count{$key},0.5);
#	$intron_binomial_p=binomial($intron_ss_overlap_count{$key},$intron_overlap_count{$key},0.5);
#	print OUT "$key\t$intron_overlap_nt{$key}\t$intron_overlap_count{$key}\t$intron_ss_overlap_nt{$key}\t$intron_ss_overlap_count{$key}\t$intron_as_overlap_nt{$key}\t$intron_as_overlap_count{$key}\t$intron_binomial_p\t$intron_total_length\t$exon_overlap_nt{$key}\t$exon_overlap_count{$key}\t$exon_ss_overlap_nt{$key}\t$exon_ss_overlap_count{$key}\t$exon_as_overlap_nt{$key}\t$exon_as_overlap_count{$key}\t$exon_binomial_p\t$exon_total_length\n";
	print OUT "$key\t$rep_class{$key}\t$rep_family{$key}\t$intron_overlap_nt{$key}\t$intron_overlap_count{$key}\t$intron_ss_overlap_nt{$key}\t$intron_ss_overlap_count{$key}\t$intron_as_overlap_nt{$key}\t$intron_as_overlap_count{$key}\t$intron_total_length\t$exon_overlap_nt{$key}\t$exon_overlap_count{$key}\t$exon_ss_overlap_nt{$key}\t$exon_ss_overlap_count{$key}\t$exon_as_overlap_nt{$key}\t$exon_as_overlap_count{$key}\t$exon_total_length\t$genome_repeat_nt{$key}\t$genome_repeat_count{$key}\n";
#	print OUT "$key\t$intron_overlap_nt{$key}\t$intron_overlap_count{$key}\t$intron_ss_overlap_nt{$key}\t$intron_ss_overlap_count{$key}\t$intron_as_overlap_nt{$key}\t$intron_as_overlap_count{$key}\t$intron_total_length\t$exon_overlap_nt{$key}\t$exon_overlap_count{$key}\t$exon_ss_overlap_nt{$key}\t$exon_ss_overlap_count{$key}\t$exon_as_overlap_nt{$key}\t$exon_as_overlap_count{$key}\t$exon_total_length\n";
	print OUTCAND "$key\n";
}
}

close (OUT);
close (OUTCAND);

##############################################################################################################################################
# create the profile output files here
# each of these profiles is filtered against the candidate list (ie only repeats that satisfy the overlap count/length criteria.

#print "\n Creating profiles\n";

#makeProfile($intersect_exon_repeat_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_exon_file";
#system($command);


#makeProfile($intersect_intron_repeat_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_intron_file";
#system($command);

#makeProfile($intersect_exon_repeat_samesense_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_exon_ss_file";
#system($command);

#makeProfile($intersect_exon_repeat_antisense_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_exon_as_file";
#system($command);

#makeProfile($intersect_intron_repeat_samesense_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_intron_ss_file";
#system($command);

#makeProfile($intersect_intron_repeat_antisense_bed, $profile_temp);
#my $command="join <\(sort $profile_temp\) <\(sort $candidate_list\) > $profile_intron_as_file";
#system($command);


############################################################
# Calculate summary statistics

my $count_genes=BEDcount($gene_bed_file);

my $length_merged_exon=BEDlength($exon_merged_bed);
my $length_stranded_intron=BEDlength($lncrna_intron_file);
my $length_repeats=BEDlength($adjusted_repeat_filename);

my $count_merged_exon=BEDcount($exon_merged_bed);
my $count_stranded_intron=BEDcount($lncrna_intron_file);
my $count_repeats=BEDcount($adjusted_repeat_filename);

my $length_overlapped_tes=BEDlength($intersect_exon_repeat_bed);
my $count_overlapped_tes=BEDcount($intersect_exon_repeat_bed);

# this command is to count the six classes of repeat intersection.
# first read them into an array
open (INBED, $intersect_exon_repeat_bed);
my @classes;
while (<INBED>) {
	chomp;
	my @line=split("\t",$_);
	push(@classes, $line[4]);
}
close(INBED);
# then count instances using a subroutine, notice you pass it an array REFERENCE.
my %intersect_class_counts=countInstances(\@classes);


open (OUT, ">$summary_stats_filename");

print OUT "\nSummary Statistics";
print OUT "\nInput Gene Annotation:\t$lncrna_gtf_file";
print OUT "\nInput Repeat Masker Annotation:\t$repeat_masker_file";
print OUT "\nTotal count of genes:\t$count_genes";

print OUT "\n\nTotal length of merged exons:\t$length_merged_exon";
print OUT "\nTotal count of merged exons:\t$count_merged_exon";
print OUT "\nTotal length of stranded introns:\t$length_merged_exon";
print OUT "\nTotal count of stranded introns:\t$count_stranded_intron";
print OUT "\nTotal length of Repeat Masker repeats:\t$length_repeats";
print OUT "\nTotal count of Repeat Masker repeats:\t$count_repeats";

print OUT "\nTotal length of intersected repeats:\t$length_overlapped_tes";
print OUT "\nTotal count of intersected repeats:\t$count_overlapped_tes";

print OUT "\n\nIntersection class counts:";
  foreach (sort keys %intersect_class_counts) {
    print OUT "\n$_ : $intersect_class_counts{$_}";
  }


print OUT "\n\n";
close (OUT);

###############################################
#calculate timing


my $end_run = time();
my $run_time = $end_run - our $start_run;
print "\nJob took $run_time seconds\n";

##############################################################################################################################################
#SUBROUTINES
##############################################################################################################################################

#make a table of counts from a two column file

sub makeTable{

	my ($inputfile, $outputfile)=@_;
	open (GIN,"$inputfile") or die "makeTable subroutine Error: Cannot open file $inputfile\n\n";
	open (GOUT, ">$outputfile") or die "makeTable subroutine Error: Cannot open file $outputfile\n\n";
#	open (GIN,"repeat_class.tmp") or die "makeTable subroutine Error: Cannot open file $inputfile\n\n";
#	open (GOUT, ">banana.txt") or die "makeTable subroutine Error: Cannot open file $outputfile\n\n";
#	print "files $inputfile, $outputfile";
	use List::MoreUtils qw/ uniq /;

	my %horse = ();
	my @barcodes = '';
	while (<GIN>) {
    		my $currLine = $_;
		chomp $currLine;
    		my ($apple, $cat) = $currLine =~ /^(.*?)[\t| ](.*?)$/;
		push(@barcodes,$apple);
	    my $temp_value = $horse{$cat}{$apple} +1;
	    $horse{$cat}{$apple} = $temp_value;
	}


	my @unique = uniq @barcodes;
	@unique = sort @unique;
	my $shit = shift @unique;

	print GOUT ".\t";

	for my $barcode (@unique){print GOUT "$barcode\t"};

	for my $fruit (sort keys %horse ) {
 	   print GOUT "\n$fruit";
	    for my $animal (@unique ) {

	       if($horse{$fruit}{$animal}){  print GOUT "\t$horse{$fruit}{$animal}"} else {print GOUT "\t0"}
		} 
	}

	close(GOUT);
	close(GIN);
}


# Send it an array REFERENCE, it gives you back a hash of counts
# eg my %intersect_class_counts=countInstances(\@classes);
sub countInstances{
	my @array = @{$_[0]};
	my %hash;
	foreach my $str (@array) {
		$hash{$str}++;
	}
	return(%hash);
}

sub BEDlength{
	my ($inbed)=@_;
	my $length;
	open (INBED, $inbed);
	while (<INBED>) {
		chomp;
		my @line=split("\t",$_);
		$length= $length+($line[2]-$line[1]);
	}

	close(INBED);
	return($length);

}

sub BEDcount{
	my ($inbed)=@_;
	my $count;
	open (INBED, $inbed);
	while (<INBED>) {
		chomp;
		my @line=split("\t",$_);
		$count++;
	}
	close(INBED);
	return($count);

}


sub splitBED{
	my ($inbed, $outssbed, $outasbed)=@_;
	open (INBED, $inbed);
	open (OUTSSBED, ">$outssbed");
	open (OUTASBED, ">$outasbed");

	while (<INBED>) {
		chomp;
		my @line=split("\t",$_);
		if($line[4]>0){print OUTSSBED "@line\n";};
		if($line[4]<0){print OUTASBED "@line\n";};
	}

	close(INBED);
	close(OUTASBED);
	close(OUTSSBED);

}

sub exonGTFtoBED{     #this accepts a GTF file, extracts exons and returns a BED where the ID column has this format ENST00000473358.1_chr1_29554_30039_+_ex1
			
	my ($gtf_file, $bed_file)=@_;
	open (GTF, $gtf_file);
	open (BED, ">$bed_file");
	while (<GTF>) {
		chomp;
		my @line=split("\t",$_);
		my @data=split(" ",$line[8]);
		my ($transid)=$data[3]=~/\"(.*?)\"/;
		my ($exon_no)=$data[17]=~/(.*?)\;/;
		my $new_id=$transid."_".$line[0]."_".$line[3]."_".$line[4]."_".$line[6];
		if($line[2] eq 'exon'){print BED "$line[0]\t$line[3]\t$line[4]\t$new_id\t1\t$line[6]\n";}
	}
close(GTF);
close(BED);
}


sub transcriptGTFtoBED{     #this accepts a GTF file, extracts exons and returns a BED where the ID column has this format ENST00000473358.1_chr1_29554_30039_+_ex1
			
	my ($gtf_file, $bed_file)=@_;
	open (GTF, $gtf_file);
	open (BED, ">$bed_file");
	while (<GTF>) {
		chomp;
		my @line=split("\t",$_);
		my @data=split(" ",$line[8]);
		my ($transid)=$data[3]=~/\"(.*?)\"/;
		my ($exon_no)=$data[17]=~/(.*?)\;/;
		my $new_id=$transid."_".$line[0]."_".$line[3]."_".$line[4]."_".$line[6];
		if($line[2] eq 'transcript'){print BED "$line[0]\t$line[3]\t$line[4]\t$new_id\t1\t$line[6]\n";}
	}
close(GTF);
close(BED);
}

sub geneGTFtoBED{     #this accepts a GTF file, extracts exons and returns a BED where the ID column has this format ENST00000473358.1_chr1_29554_30039_+_ex1
			
	my ($gtf_file, $bed_file)=@_;
	open (GTF, $gtf_file);
	open (BED, ">$bed_file");
	while (<GTF>) {
		chomp;
		my @line=split("\t",$_);
		my @data=split(" ",$line[8]);
		my ($transid)=$data[3]=~/\"(.*?)\"/;
		my ($exon_no)=$data[17]=~/(.*?)\;/;
		my $new_id=$transid."_".$line[0]."_".$line[3]."_".$line[4]."_".$line[6];
		if($line[2] eq 'gene'){print BED "$line[0]\t$line[3]\t$line[4]\t$new_id\t1\t$line[6]\n";}
	}
close(GTF);
close(BED);
}

#repIntersectBED performs an intersect between a repeat bed file, an exon bed file, filters the result and returns another bed as well as a GTF.
#pass it: repeat.bed,lncrna.exon.bed,minimum.intersect.length,output.bed.file,output.gtf.file
#important: the output BED score indicates what time of intersection it it, the sign indicates whether the transcript and repeat are on same or opposite strand.
# 1=repeat overlaps start of first exon, TSS
# 2=repeat overlaps start of internal exon (or last exon) but not end
# 3=repeat overlaps start and end of internal exon
# 4=repeat overlaps end only of internal exon (or first exon)
# 5=repeat overlap end of last exon, TTS
# 6=repeat lies within the exon
sub repIntersectBED{  
	my ($repeat_bed_file,$exon_bed_file,$transcript_bed_file,$minimum_intersect,$output_bed,$output_gtf)=@_;
	my %classifications;    #hash holding the classifier id, key=id of repeat
	open(OUT,">$output_bed");
	open(OUTGTF,">$output_gtf");	   #will write a GTF format output too
	
	#here we are using the transcript annotation to find those repeats overlapping TSS and TTS, so we overlap with transcripts. any overlap here prioritises over an overlap with part of another transcript
	my @transoverlap=`intersectBed -a $repeat_bed_file -b $transcript_bed_file -wb `;   
	foreach my $panda (@transoverlap){
		my @koala=split("\t",$panda);
		my $temp_strand=$koala[11];
		chomp $temp_strand;
		if($temp_strand eq "+"){   #check strand of the transcript

			if($koala[1]==$koala[7]){$classifications{$koala[3]}=1;};    #assign as TSS overlap
			if($koala[2]==$koala[8]){$classifications{$koala[3]}=5;};    #assign as TTS overlap
		}
		if($temp_strand eq "-"){
			if($koala[1]==$koala[7]){$classifications{$koala[3]}=5;};    
			if($koala[2]==$koala[8] && $koala[2]>$koala[7]){$classifications{$koala[3]}=1;};
		}
	}


	my @results=`intersectBed -a $repeat_bed_file -b $exon_bed_file -wb | awk '\$3-\$2>=$minimum_intersect' `;   #do the intersect of repeats with exons and save in results array
	my $resultcount=scalar(@results);

	foreach my $line (@results){    #
		chomp $line;
		my @shrimp=split("\t",$line);
		my($repeat_start,$repeat_end,$repeat_strand,$repeat_name,$repeat_family,$repeat_class,$repstart,$repend,$repleft,$hypothetical_start,$hypothetical_end)=$shrimp[3]=~/^.*?,(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?)$/;

		my $strandedness;    #this is + if theyre both on same strand, - if opposite
		if($repeat_strand eq $shrimp[11]){$strandedness="+"} else {$strandedness="-"}
		my $overlap_type;

		if($classifications{$shrimp[3]}){
			$overlap_type=$strandedness.$classifications{$shrimp[3]};
			print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n"; 
			print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";
			next;}

		if($shrimp[1]>$shrimp[7] && $shrimp[2]<$shrimp[8]){
			$overlap_type=$strandedness."6";
			print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n"; 
			print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";
			next};

		if($shrimp[1]==$shrimp[7] && $shrimp[2]==$shrimp[8]){  #if the two beds overlap exactly it means that the repeat is actually encompassing the exon
			$overlap_type=$strandedness."3";
			print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n"; 
			print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";
			next};



		if($shrimp[11] eq "+"){
			if($shrimp[1]==$shrimp[7]){
				$overlap_type=$strandedness."2";
				print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n";
				print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";};    
			if($shrimp[2]==$shrimp[8]){
				$overlap_type=$strandedness."4";
				print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n";
				print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";
				};
		}
		if($shrimp[11] eq "-"){
			if($shrimp[1]==$shrimp[7]){
				$overlap_type=$strandedness."4";
				print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n";
				print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";			
				};    
			if($shrimp[2]==$shrimp[8]){
				$overlap_type=$strandedness."2";
				print OUT "$shrimp[0]\t$shrimp[1]\t$shrimp[2]\t$shrimp[3]\t$overlap_type\t$shrimp[5]\n";
				print OUTGTF "$shrimp[0]\ttransposon.profile\texonic.repeat\t$shrimp[1]\t$shrimp[2]\t.\t$shrimp[5]\t.\ttranscript=$shrimp[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$overlap_type\n";			
				};
		}

	} 






close (OUT);
close (OUTGTF);
}










sub repeatGTF{	      #pass it the repeat BED filename and the exon BED filename and it gives back the GTF, one entry per overlap event. ie nonredundant for transcripts.
	
	my($repeat_bed_file, $exon_bed_file, $output_gtf) = @_;
#	print "\n\nsubroutine print $repeat_bed_file, $exon_bed_file, $output_gtf\n\n";
	open (OUT, ">$output_gtf");
	
	my @test=`intersectBed -a $repeat_bed_file -b $exon_bed_file -wb`;

	foreach (@test) {
		chomp $_;
		my @line=split('\t',$_);
#		my @transcripts=split(';',$line[9]);
#		foreach my $trans (@transcripts) { 
#			my($trans_id, $exon_start,$exon_end,$exon_strand)=$trans=~/^(.*?)_.*?_(.*?)_(.*?)_(.*?)$/;
			my($repeat_start,$repeat_end,$repeat_strand,$repeat_name,$repeat_family,$repeat_class,$repstart,$repend,$repleft,$hypothetical_start,$hypothetical_end)=$line[3]=~/^.*?,(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?),(.*?)$/;
			print OUT "$line[0]\ttransposon.profile\texonic.repeat\t$line[1]\t$line[2]\t.\t$line[5]\t.\ttranscript=$line[9];repeat_start=$repeat_start;repeat_end=$repeat_end;repeat_strand=$repeat_strand;repeat_name=$repeat_name;repeat_family=$repeat_family;repeat_class=$repeat_class;repstart=$repstart;repend=$repend;repleft=$repleft;repeat_hypothetical_start=$hypothetical_start;repeat_hypothetical_end=$hypothetical_end;overlap_type=$line[4]\n";
#		} 
	}
	close(OUT);

}


sub binomial {
    my $k = shift(@_); #number of successes
    my $n = shift(@_); #number of trials
    my $p = shift(@_); #prob of success
    my $prob;

    $prob = ($p**$k) * ((1 - $p)**($n - $k)) * &factorial($n) / (&factorial($k) * &factorial($n - $k));

    return $prob;
}

sub factorial {
    my $n = shift(@_);
    my $fact = 1;

    if (($n < 0) or (170 < $n)) {
	die "Factorial out of range";
    }

    for(my $i = 1; $i <= $n; $i++) {
	$fact *= $i;
    }

    return $fact;
}


sub calcRcorrel{ #give it two profile files, it calls the R correlation script on them and returns a hash with the calculated values.
	my($first_filename, $second_filename) = @_;
	my $temp_filename='temp.correl.out';  #this is defined inside the R script
	my %correl_values={};
	system ("Rscript profile.correlation.1.R $first_filename $second_filename");
	
	open (FILE, $temp_filename);
		while (<FILE>) {
			chomp;
			my @line=split" ";
			$correl_values{$line[0]}=$line[1];
	}
	return %correl_values;
}


sub bedOLnt{   #you give it a bed, it counts the number of nt overlap per repeat class, returns this as a hash.
	my %nt_overlap={};   #hash holding the number of nt of overlap for each repeat class.	
	my($bed_filename) = @_;
	open (FILE, $bed_filename);
	while (<FILE>) {
		chomp;
		my @line=split" ";
		my ($rep_class)=$line[3]=~/chr.*?\,.*?\,.*?\,.*?\,(.*?)\,/;
		$nt_overlap{$rep_class}+=$line[2]-$line[1];
	}
return %nt_overlap;
}



sub bedOLcount{   #you give it a bed, it counts the number of times you get overlap per repeat class, returns this as a hash.
	my %count_overlap={};   #hash holding the number of nt of overlap for each repeat class.	
	my($bed_filename) = @_;
	open (FILE, $bed_filename);
	while (<FILE>) {
		chomp;
		my @line=split" ";
		my ($rep_class)=$line[3]=~/chr.*?\,.*?\,.*?\,.*?\,(.*?)\,/;
		$count_overlap{$rep_class}++;
	}
#my $hash=$count_overlap{'AluSg'};
#print " $bed_filename hash $hash\n";
return %count_overlap;
}




sub makeProfile{
	my($filename, $output_profile_file) = @_;
	open (OUT, ">$output_profile_file");
	open (FILE, $filename);
		my %hoa={};
		my %total={};
	 while (<FILE>) {
 		chomp;
		my @line=split" ";
		my $n='';
		my $m='';
		my $width=10000;

		# this command gives overlap for the whole of repeats that intersect somewhere on a lncRNA
		#my ($rep_coord_start, $rep_coord_end,$repstrand,$repclass, $repstart, $repend,$repleft,$hypoth_genome_start,$hypoth_genome_end) =$line[3]=~/^chr.*,(\d*),(\d*),([+-]),(.*?),.*,.*,(.*),(.*),(.*),(.*),(.*)$/; 


		# these definitions are the original ones, where we're only drawing profiles for tha part of the 
		my ($repstrand,$repclass, $repstart, $repend,$repleft,$hypoth_genome_start,$hypoth_genome_end) =$line[3]=~/^chr.*,\d*,\d*,([+-]),(.*?),.*,.*,(.*),(.*),(.*),(.*),(.*)$/; my $rep_coord_start=$line[1]; my $rep_coord_end=$line[2];
		
		

		if($repstrand eq "+"){
			my $counter=0;
			$m=0;
			$n=$hypoth_genome_start; 
			$total{$repclass}++;

				while ($n<$hypoth_genome_end){ 
					if($n<$rep_coord_start){$n++; $counter++; }; 
					if($n>=$rep_coord_start && $n<$rep_coord_end){$hoa{$repclass}[$counter]++;  $counter++;$n++;}; 
					if($n>=$rep_coord_end){$n++; $counter++; } 
					}
				}


		if($repstrand eq "-"){
			my $counter=0;
			$n=$hypoth_genome_end; 
			$total{$repclass}++;

				while ($n>=$hypoth_genome_start){ 
					#print "counter $counter n $n hypothstart $hypoth_genome_start hypothend $hypoth_genome_end recordstart $rep_coord_start repcordend $rep_coord_end\n";
					if($n>=$rep_coord_end){$n--; $counter++; }; 
					if($n>=$rep_coord_start && $n<$rep_coord_end){$hoa{$repclass}[$counter]++;  $counter++; $n--;}; 
					if($n<$rep_coord_start){$n--; $counter++; } 
					}
				}
		}

		
			foreach my $key ( sort (keys %hoa )){
				if($hoa{$key}){			#this gets rid of annoying empty cells in hash
				my $lcounter=0;
				my $motif_length=scalar(@{$hoa{$key}});
				print OUT "\n$key\t$motif_length\t$total{$key}"; 
				foreach(@{$hoa{$key}}){
					my $value=($_/$total{$key}); 
					print OUT "\t$value";
					$lcounter++;
					}
				while($lcounter<10000){     #this fills out each line to equal width of 10000 with NAs, because R won't accept data of differing widths
					print OUT "\tNA";
					$lcounter++;
				}
				}
		}
	close(FILE);
	close(OUT);
}






