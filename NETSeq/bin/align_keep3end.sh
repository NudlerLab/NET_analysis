start_time=`date +%s`
bowtie -S -m 1 -v 0 indexes/NC_000913.2 netseq_data/clean_Trimmed_wt-mmc.fastq.gz 2> wt_mmc_NET.log|samtools view -b -h -F 4|bedtools bamtobed -i '-' > wt_mmc_NET_aligned_2.bed
grep -w '+' wt_mmc_NET_aligned_2.bed|awk '{print $1"\t"$2"\t"($2+1)"\t"".""\t"$5"\t"$6}'> plus_3end_2.bed
grep -w '-' wt_mmc_NET_aligned_2.bed|awk '{print $1"\t"($3-1)"\t"($3)"\t"".""\t"$5"\t"$6}'> minus_3end_2.bed && echo run time is $(expr `date +%s` - $start_time) s
