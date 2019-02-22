RUN=$1

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed_data/${RUN}_R1_trimmed.fastq.gz -p trimmed_data/${RUN}_R2_trimmed.fastq.gz ${RUN}/R1.fastq.gz ${RUN}/R2.fastq.gz > trimmed_data/output_files/${RUN}_cutadapt_output.txt &


