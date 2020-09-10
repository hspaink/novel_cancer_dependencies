import os
import sys
import numpy as np
import pandas as pd
from taigapy import TaigaClient

PEDIATRIC_CANCERS = ['Ewing_sarcoma',
                     'AML','neuroblastoma',
                     'rhabdomyosarcoma',
                     'osteosarcoma',
                     'malignant_rhabdoid_tumor',
                     'medulloblastoma',
                     't-ALL',
                     'b-ALL']

def remap_index(df, index_map, df_idx_name='index', drop_nans=True):
    df[df_idx_name] = [index_map.loc[index_map == idx].index[0] 
                       if idx in index_map.values else np.nan 
                       for idx in df.index]
    df.set_index(df_idx_name, inplace=True)
    if drop_nans:
        df = df.loc[~df.index.isna()]
    return df

def get_from_taiga(name, version, file, split_attribute=None, col=None):
    assert split_attribute in [None, 'header', 'column']
    tc = TaigaClient()
    
    data = tc.get(name=name, version=version, file=file)
    
    if split_attribute == 'header':
        data = data[[col for col in data.columns if '&' not in col]]  # In some data entries containing two genes exist
        head = pd.DataFrame([data.columns.str.split(' ',1).str[0].astype(str), 
                             data.columns.str.split(' ',1).str[1].str.strip('()').astype(int)], 
                            index=['gene','geneID'], columns=data.columns)
        data = pd.concat([head, data], join='inner').T.set_index(['geneID', 'gene']).T
    
    elif split_attribute == 'column': 
        assert col is not None
        data = pd.concat([pd.DataFrame(data[col].str.split(' ',1).tolist(), columns=['gene','geneID']),
                          data.loc[:, data.columns != col]], axis=1, join='inner').set_index('geneID')
        data.index = data.index.str.strip('()').astype(int)
        
        
    return data
        

