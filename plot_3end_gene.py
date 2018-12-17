import argparse
import pandas as pd
import matplotlib
matplotlib.use('TkAgg') ##NOTE this line is important if running on MAC
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def make_counts(file_lis):
	#for f in args.infile:
	for f in file_lis:
		print('plus' in f)
		if 'plus' in f:
			p_3end=pd.read_table(f,header=None)
		else:
			print('minus' in f)
			if 'minus' in f:
				m_3end=pd.read_table(f,header=None)
			else:
				raise Exception('There must be one plus file and one minus file. Did not find file with plus and/or minus in the file name')
	
	p_counts=pd.DataFrame(p_3end[2].value_counts()).sort_index().reindex(range(1,4641653),fill_value=0)
	m_counts=pd.DataFrame(m_3end[2].value_counts()).sort_index().reindex(range(1,4641653),fill_value=0)

	return p_counts,m_counts

def get_pos(g_file,gene):
	gene_annot=pd.read_csv(g_file,sep='\t',header=None)
	gene_strand=gene_annot[gene_annot[4]==gene].iloc[0][5]
	gene_pos=((gene_annot[gene_annot[4]==gene].iloc[0][1])-12,((gene_annot[gene_annot[4]==gene].iloc[0][2])+13))

	return gene_pos,gene_strand

def make_plot(gene_name,p_df,m_df):
	
	gene_pos,gene_strand=get_pos(args.genome,gene_name)

	if gene_strand == '+':
		gene_plot=m_df.iloc[gene_pos[0]:gene_pos[1],:]
	elif gene_strand == '-':
		gene_plot=p_df.iloc[gene_pos[0]:gene_pos[1],:]
	
	x=gene_plot.index.tolist()
	y=gene_plot[2].tolist()

	if gene_strand == '-':
		y=y[::-1]
	else:
		pass

	plt.figure(figsize=(15,2))
	sns.barplot(x=x,y=y,color='blue',alpha=0.7)
	sns.despine(top=False,right=False,bottom=False,left=False)
	#plt.xticks(np.arange(min(x), max(x)+1, 100))
	#plt.ylim(0,600)
	#plt.yticks([0,300,600])
	plt.xticks([])
	plt.savefig(gene_name+'_test_3.pdf', bbox_inches='tight')
	

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Two input files required",nargs=2,type=str)
parser.add_argument("-gf", "--genome", help="bed file of genome",type=str)
parser.add_argument("-g", "--genes", help="genes of interest",nargs= '*',type=str)
args = parser.parse_args()

p_df,m_df=make_counts(args.infile)
for gene_name in args.genes:
	make_plot(gene_name,p_df,m_df)