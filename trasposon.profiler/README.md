# Transposon Profiler
## A program to annotate Transposable and Repeat Elements in LncRNA
#
##### Version 17
##### Rory Johnson 
#
### Typical command:
##### perl transposon.profiler.16.pl -e input.files/gencode.v18.long_noncoding_RNAs.gtf -r input.files/repeat.masker.hg19.test -o test17.out
#####
##### Inputs are a GTF annotation and a RepeatMasker file
### 
### Explanation of overlap classes:
##### +1 "tss.plus"
##### +2 "splice_acc.plus"
##### +3 "encompass.plus"
##### +4 "splice_don.plus"
##### +5 "tts.plus"
##### +6 "inside.plus"
