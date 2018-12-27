import sys
import csv
import numpy as np
import pandas as pd

df_ecocyc = pd.read_table(
        'All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.csv',
        sep=r',', skipinitialspace=True)

df_genes_ecocyc = pd.DataFrame(df_ecocyc.loc[:,[
                                'Genes',
                                'Centisome-Position',
                                'Left-End-Position',
                                'Right-End-Position',
                                'Transcription-Direction']])

df_expression = pd.read_csv('Expression_Data_MG1655_sorted_RPKM.csv')
df_genes_and_expression = pd.merge(df_expression, df_genes_ecocyc,
                how='left', left_on='Gene', right_on='Genes')

dfs_list = []
for i in range(1,16):
	df_syno_ecocyc = pd.DataFrame(df_ecocyc.loc[:,[
                                        'Centisome-Position',
                                        'Left-End-Position',
                                        'Right-End-Position',
                                        'Transcription-Direction',
                                        'synonyms{0}'.format(str(i))]])
	df_genes_and_expression_add_syno = pd.merge(df_genes_and_expression, df_syno_ecocyc,
                        how='left',left_on='Gene',right_on=df_syno_ecocyc.columns[4].lstrip().rstrip())
	dfs_list.append(df_genes_and_expression_add_syno)

df_all_names_and_expression = pd.concat(dfs_list)
df_only_syno_expression = df_all_names_and_expression.drop_duplicates(['Centisome-Position_y'])
df_only_syno_expression_cleaned = df_only_syno_expression[
        ['Gene',
        'Centisome-Position_y',
        'Left-End-Position_y',
        'Right-End-Position_y',
        'Transcription-Direction_y',
        'Translation efficiency (AU)',
        'mRNA level (RPKM)']
        ].reset_index().drop(['index'],1).drop_duplicates(['Gene'])

df_only_syno_expression_cleaned2 = df_only_syno_expression_cleaned.drop(
        df_only_syno_expression_cleaned.index[df_only_syno_expression_cleaned.Gene == 'lpp'])

df_only_syno_expression_cleaned3 = df_only_syno_expression_cleaned2.drop(
        df_only_syno_expression_cleaned2.index[df_only_syno_expression_cleaned2.Gene == 'ade'])

df_only_syno_expression_del_ambig = df_only_syno_expression_cleaned3.drop(
        df_only_syno_expression_cleaned3.index[df_only_syno_expression_cleaned3.Gene == 'thiB'])

df_only_syno_expression_del_ambig.rename(columns={
        'Centisome-Position_y': 'Centisome-Position',
        'Left-End-Position_y': 'Left-End-Position',
        'Right-End-Position_y':'Right-End-Position',
        'Transcription-Direction_y':'Transcription-Direction'}, inplace=True)

df_only_syno_expression_del_ambig = df_only_syno_expression_del_ambig[
        ['Gene',
        'mRNA level (RPKM)',
        'Translation efficiency (AU)',
        'Centisome-Position',
        'Left-End-Position',
        'Right-End-Position',
        'Transcription-Direction']]

df_genes_and_expression_inner = pd.merge(df_expression, df_genes_ecocyc, how='inner', left_on='Gene', right_on='Genes')
chuncks = [df_only_syno_expression_del_ambig, df_genes_and_expression_inner]
df_expression_table_full = pd.concat(chuncks).reset_index().drop(['index','Genes'],1)

##df_expression_table_full.to_csv('expression_table.csv') #expression_table_full from now on = etf
df_leading = pd.DataFrame(df_expression_table_full[
        (((df_expression_table_full['Left-End-Position'] > 0)
        & (df_expression_table_full['Left-End-Position'] < 1588710)
        & (df_expression_table_full['Transcription-Direction'] == '+'))
        | ((df_expression_table_full['Left-End-Position'] >3923850)
        & (df_expression_table_full['Left-End-Position'] < 4639676)
        | (df_expression_table_full['Transcription-Direction'] == '+')))
        & ((df_expression_table_full['Left-End-Position'] > 1588710)
        & (df_expression_table_full['Left-End-Position'] < 3923850)
        (df_expression_table_full['Transcription-Direction'] == '-'))])

df_leading['leading'] = '+'
df_leading['lagging'] = '-'

df_lagging = pd.DataFrame(df_expression_table_full[
        (((df_expression_table_full['Left-End-Position'] > 0)
        & (df_expression_table_full['Left-End-Position'] < 1588710)
        & (df_expression_table_full['Transcription-Direction'] == '-'))
        | ((df_expression_table_full['Left-End-Position'] > 3923850)
        & (df_expression_table_full['Left-End-Position'] < 4639676)
        & (df_expression_table_full['Transcription-Direction'] == '-')))
        | ((df_expression_table_full['Left-End-Position'] > 1588710)
        & (df_expression_table_full['Left-End-Position'] < 3923850)
        & (df_expression_table_full['Transcription-Direction'] == '+'))])

df_lagging['leading'] = '-'
df_lagging['lagging'] = '+'
chuncks_collisions = [df_leading, df_lagging]
df_etf_collision = pd.concat(chuncks_collisions)
df_etf_collision.to_csv('expression_table_collisions.csv')
#This code construct a dataframe with all expressed genes from E.coli MG1655.
#'All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.csv' is a smart table imported from ecocyc.
#'Expression_Data_MG1655_sorted_RPKM.csv' is my edited version of a file with list of all expressed genes with RNA-seq RPKM and ribo-seq translation-effciency values.
#The code merges the two files to give columns of: Gene name, centisome position, left coordinate,right coordinate, plus/minus strand, RPKM, translation A.U.
#challenge was to match different names for the same genes in the two different lists.
# line 5: read_table is used instead of read_csv. Reason - wanted to remove the white sapce preceding the synonym name of each gene.
# line 6: remove the columns'common-name' and all synonyms columns. did not want crazy duplication of these columns with each merge in the loop.
# line 8: first merge based on Genes column in smart table and Gene columns in expression data. merge on expression data to retain info of all expressed genes
# line 10: starts a loop that constructs multiple data frames for each synonym column plus the rest of the columns that appear in the first merged df.
# line 11: creates a list to hold all the new dfs.
# line 12: merges synonyms1 df with the first merge df.
# line 13: append the merged df to df list
# loop: repeat process for all synonyms df up until 15.
# line 14: concat all dfs in the list.
