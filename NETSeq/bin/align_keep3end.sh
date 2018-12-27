#!/bin/bash

###############################################################################
# NET-Seq analysis pipeline
#
# The script assumes that sample .fastq files are organized as follows:
#
# ..-data
#    |
#    |-sample1
#    | |-sample1.fastq.gz
#    |
#    |-sample2
#    | |-sample2.fastq.gz
#    .....
#    and so on.
###############################################################################

index="../ref/MG1655"
base_data="../data"
base_results="../results"

NUMPROC=30 # For snowflake, adjust as needed
READLEN=19 # Trim reads to this length

start_time=`date +%s`

for sample in $(ls $datadir)
do
    sampledir="$base_data/$sample"
    resultdir="$base_results/$sample"
    bedfile="$resultdir/$sample/${sample}_aligned.bed"
    zcat $sampledir/$sample.fastq.gz \
        | bowtie2 -p $NUMPROC --trim-to $READLEN -x $index -U - \
        | samtools view -bhS -F 4 \
        | bedtools bamtobed -i '-' \
        > $bedfile

    grep -w '+' $bedfile \
        | awk '{print $1"\t"$2"\t"($2+1)"\t"".""\t"$5"\t"$6}' \
        > $resultdir/plus_3end_2.bed
    
    grep -w '-' $bedfile \
        | awk '{print $1"\t"($3-1)"\t"($3)"\t"".""\t"$5"\t"$6}' \
        > $resultdir/minus_3end_2.bed 
    
done

echo "run time is $(expr `date +%s` - $start_time) s"

#bowtie -S -m 1 -v 0 indexes/NC_000913.2 netseq_data/clean_Trimmed_wt-mmc.fastq.gz 2> wt_mmc_NET.log|samtools view -b -h -F 4|bedtools bamtobed -i '-' > wt_mmc_NET_aligned_2.bed
#grep -w '+' wt_mmc_NET_aligned_2.bed|awk '{print $1"\t"$2"\t"($2+1)"\t"".""\t"$5"\t"$6}'> plus_3end_2.bed
#grep -w '-' wt_mmc_NET_aligned_2.bed|awk '{print $1"\t"($3-1)"\t"($3)"\t"".""\t"$5"\t"$6}'> minus_3end_2.bed && echo run time is $(expr `date +%s` - $start_time) s
