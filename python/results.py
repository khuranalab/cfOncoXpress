# Weiling Li (wel4007@med.cornell.edu)
# correlation results for 13 samples with matched RNA-seq

import collections
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import scipy.spatial.distance
import scipy.signal
import time
import sys
from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error


def resmat(predmat,truthmat,samplename):
    mynewmatrix=pd.DataFrame(columns=['gene_id','gene_name', 'TPM', 'correction'])
    for i,row in predmat.iterrows():
        res=truthmat[(truthmat['IDENTIFIER']==row['gene_name']) & (truthmat['GENE_ID']==row['gene_id'])][samplename]
        if len(res)==1:
            tpm=res.values[0]
            newrow = {'gene_id':row['gene_id'], 'gene_name': row['gene_name'],'TPM':tpm, 'correction':row['correction']}
            new_df = pd.DataFrame([newrow])
            mynewmatrix = pd.concat([mynewmatrix, new_df], axis=0, ignore_index=True)        
    return mynewmatrix

mat=pd.read_csv('../csv_files/rna.csv',sep='\t')

outdir=f'csv_files/prediction' # 
for samplename in mat.columns[1:-1].values:
    print(samplename)
    pred=pd.read_csv(f'{outdir}/{samplename}_correction.csv')[['gene_id','gene_name','correction']]
    newmat=resmat(pred,mat,samplename)
    newmat.to_csv(f'{outdir}/{samplename}_correction_withTPM.csv')  # with TPM
    r, p = scipy.stats.spearmanr(newmat['TPM'],newmat['correction']) #y_true, y_pred
    print(r,p)
    mse = mean_squared_error(newmat['TPM'],newmat['correction'])
    print(mse)
    
