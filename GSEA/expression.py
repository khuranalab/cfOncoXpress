# Weiling Li (wel4007@med.cornell.edu)
# generating an expression matrix for GSEA

import numpy as np
import pandas as pd



def combinedf(path,names):
    c1=pd.read_csv(f'{path}/{names[0]}_correction.csv')[['gene_id','gene_name','correction']]
    c1.columns=['gene_id','gene_name',names[0]]
    c2=pd.read_csv(f'{path}/{names[1]}_correction.csv')[['gene_id','gene_name','correction']]
    c2.columns=['gene_id','gene_name',names[1]]
    cc=pd.merge(c1,c2, on=['gene_id','gene_name']) #
    for i in range(2,len(names)):
        c=pd.read_csv(f'{path}/{names[i]}_correction.csv')[['gene_id','gene_name','correction']]
        c.columns=['gene_id','gene_name',names[i]]
        cc = pd.merge(cc,c, on=['gene_id','gene_name'])
    
    namelist=['gene_id','gene_name']
    for i in range(0,len(names)):
        namelist.append(names[i])
    cc.columns=namelist

    return cc

cfoncoXpress_dir='csv_files/prediction'
herberts=pd.read_csv('csv_files/herberts_meta.txt',sep='\s+',header=None)[1].values
correction = combinedf(cfoncoXpress_dir, herberts)


pbmctpm=pd.read_csv('csv_files/controls34TPM.csv')
pbmclog2tpm=pbmctpm.drop(columns=['gene_id','gene_name'])
pbmclog2tpm = np.log2(pbmclog2tpm+1)
pbmclog2tpm['gene_id']=pbmctpm['gene_id']
pbmclog2tpm['gene_name']=pbmctpm['gene_name']
#pbmclog2tpm


names=pd.read_csv('csv_files/MBCnames',header=None)[0]
names26 = names[names != 'SRR12587654'].reset_index(drop=True)
mbc26=combinedf('csv_files/prediction', names26)


##########
merged_df = pd.merge(pbmclog2tpm, correction, on=['gene_id', 'gene_name'])
merged_df = pd.merge(merged_df, mbc26, on=['gene_id', 'gene_name'])
mydf=merged_df.drop(columns=['gene_id','gene_name'])
mydf.index=merged_df['gene_name']

df_min_max_scaled = mydf.copy() 
  
# apply normalization techniques 
for column in df_min_max_scaled.columns: 
    df_min_max_scaled[column] = (df_min_max_scaled[column] - df_min_max_scaled[column].min()) / (df_min_max_scaled[column].max() - df_min_max_scaled[column].min())     
  
df_min_max_scaled.insert(loc=0, column='DESCRIPTION', value=['NA']*df_min_max_scaled.shape[0])
df_min_max_scaled.insert(loc=0, column='NAME', value=df_min_max_scaled.index)
df_min_max_scaled.to_csv('pbmc.prostate.breast26.minmax.txt',index=False, sep='\t')

print(df_min_max_scaled.shape)
