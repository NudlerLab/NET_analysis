
# coding: utf-8

# In[19]:


import sys
import os
import re
import gzip


# In[24]:


from itertools import zip_longest

def grouper(infile, n, fillvalue=None):
    args = [iter(infile)] * n
    return zip_longest(*args, fillvalue=fillvalue)

infile=os.path.join('../input_files_NETseq/wt-mmc.fastq.gz')
outfile=os.path.join('../input_files_NETseq/Trimmed_wt-mmc.fastq.gz')
with gzip.open(infile,'rt') as infile, gzip.open(outfile,'wt') as outfile:
    fastq_per_read=grouper(infile,4,',')
    for read in fastq_per_read:
        header=read[0]
        seq=read[1]
        sign=read[2]
        quality=read[3]
        trimmed_seq=seq[0:19]
        trimmed_quality=quality[0:19]
        outfile.write(header + trimmed_seq + '\n' + sign + trimmed_quality + '\n')


# In[25]:


infile=os.path.join('../input_files_NETseq/Trimmed_wt-mmc.fastq.gz')
outfile=os.path.join('../input_files_NETseq/clean_Trimmed_wt-mmc.fastq.gz')
with gzip.open(infile,'rt') as infile:
        lines = infile.readlines()
with gzip.open(outfile, 'wt') as outfile:
        lines = filter(lambda x: x.strip(), lines)
        outfile.writelines(lines)


# In[ ]:


# next, run bowtie 1.0 with 'clean_Trimmed_wt-mmc.fastq.gz', following command:
# ./bowtie -m 1 -v 0 U000963_index/U000963_index --suppress 1,3,6,7,8 U000963_reads/clean__Trimmed_wt-mmc.fastq.gz > wt_mmc_NET.map
#use the amp file as input for the next jupyter notebook

