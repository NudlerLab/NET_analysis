
# coding: utf-8

# In[1]:


import sys
import os
import numpy as np
import pandas as pd
from collections import Counter
from operator import itemgetter


# In[2]:


def rawdata(map_file):
    lis_dicts=[]
    pDict = {}
    mDict = {}
    infile = open(map_file, 'r')
    line = infile.readline()
    while line != '':
        fields = line.split()
        col1 = str(fields[0])   #strand; note: if sequencing was performed without barcode reading, the column numbering is changed
        col2 = int(fields[1])   #left-most position
        col3 = str(fields[2])   #footprint seq
        length = len(col3)      #footprint length        
        if col1 == '+': #for plus strand
            pend5 = col2 + 1                #Bowtie uses 0-based offset: transform to 1-based and subtract 1st base 
            if pend5 in pDict:
                pDict[pend5] += 1.0
            else:
                pDict[pend5] = 1.0 
        elif col1 == '-':               #for minus strand
            end3 = col2 + 1         #for minus strand, Bowtie gives leftmost position (3' end) with zero-based numbering
            mend5 = end3 + length - 1
            if mend5 in mDict:
                mDict[mend5] += 1.0
            else:
                mDict[mend5] = 1.0 
        else:
                pass
        line = infile.readline()
    for pend5 in range(1,4641653):          
        if pend5 not in pDict:
            pDict[pend5] = 0
    p_list=[(p, pDict[p]) for p in sorted(pDict)]
    for mend5 in range(1,4641653):          
        if mend5 not in mDict:
            mDict[mend5] = 0
    m_list=[(m, mDict[m]) for m in sorted(mDict)]
    infile.close()
    lis_dicts.append(p_list)
    lis_dicts.append(m_list)
    return lis_dicts

def get_local_max_signal(lis_dicts):
    lis_pause_lis=[]
    for lis in lis_dicts:
        numeric_list=[(int(pos),float(count))for pos,count in lis]
        trunc_list=numeric_list[0:100000]
        pause_list=[]
        win_side=25
        for i in range(win_side,100000-win_side): #range from list_m[0](which equals pos 1 in genome) to list_m[4639651] (not including list_m[4639651] by default)
            window=trunc_list[(i-win_side):(i+win_side+1)] #1st window: list_m[0]-list_m[51](not including list[51] by default)
            dic_freq=Counter(mem[1] for mem in window) #count the number oftimes each signal appear in window. goal: filter regions with multiple equal local max signal
            win_count_lis=[count for pos,count in window]
            win_count_mean=np.array(win_count_lis).mean()
            win_count_std=np.array(win_count_lis).std()
            win_count_max=max(window,key=itemgetter(1))[1]
            if window[win_side][1]==win_count_max and dic_freq[window[win_side][1]]==1 and win_count_max>=(win_count_mean + win_count_std*5):  #identify window max for signal of reads.
                pause_list.append(window[win_side])
            else:
                pass
        lis_pause_lis.append(pause_list)
    return lis_pause_lis


# In[3]:


map_file=os.path.join('wt_mmc_NET_NC_000913.2.map')
lis_dicts=rawdata(map_file)


# In[4]:


lis_dicts


# In[5]:


lis_pause_lis=get_local_max_signal(lis_dicts)
print (lis_pause_lis[1][0:100])

