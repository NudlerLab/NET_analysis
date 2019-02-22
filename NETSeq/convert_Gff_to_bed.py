
# coding: utf-8

# ### Convert a GFF3 downloaded from NCBI to a BED file. (Using E. coli ref genome as an example)

# In[46]:


import pandas as pd


# #### Read in the file using pandas. The first 3 lines of the file can be skipped because they do not contain any relevant gene/CDS information

# In[47]:


df=pd.read_csv('NC_000913.2.gff3.txt',sep='\t',header=None,skiprows=3) #Skip first three lines in the file
df.drop(1,axis=1,inplace=True)  #drop column 1 
df.head()


# #### The file has duplicates because it lists both genes and CDS. For this file I only want to keep the genes

# In[48]:


df_gene=df.loc[df[2] == 'gene'].copy()
df_gene.head()


# #### The info for gene names is contained in column 8, so we want to extract the gene name for each gene interval and make a new column with gene names

# In[49]:


df_gene.loc[:,('names')]=df_gene.apply(lambda x:x[8].split(';')[2].split('=')[1],axis=1).values
df_gene.head()


# #### Drop more of the unneccesary columns....

# In[50]:


df_gene.drop([5,7,8],axis=1,inplace=True)
df_gene.head()


# #### BED files are specific in their column order, so we change the order of our columns so that it is in the correct format

# In[51]:


cols=[0,3,4,2,'names',6]
df_gene=df_gene[cols]
df_gene.head()


# #### BED files are also 0 indexed where GFF are 1 indexed, so we subtract one from the first coordinate
# ### WARNING: Only run this cell once or you'll end up subtracting more than one from the coordinate

# In[52]:


df_gene[3]=df_gene[3]-1
df_gene.head()


# ### Save the file as a BED file with all gene intervals

# In[53]:


df_gene.to_csv("NC_000913.2.bed",sep='\t',header=False,index=False)

