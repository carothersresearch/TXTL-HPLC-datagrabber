import pandas as pd
import os
from subprocess import call
import numpy as np


def get_most_updated(file_list):
        assert file_list
        if len(file_list) > 1:
            # split the file name by the period in the extension
            versions = [i.split('.')[0].split('-')[-1] for i in file_list]
            versions = [i for i in versions if '/' not in i]
            versions.sort()
            latest = [i for i in file_list if ('-' + versions[-1]) in i]
            return latest[0]
        else:
            return file_list[0]


def import_init_conc_data(EXP_DATE):
    
    # List the files in the folder
    file_names = [i for i in os.listdir(EXP_DATE) if "~" not in i]
    # load the folders:
    genex_mt=[(EXP_DATE+'/'+i) for i in file_names if 'genex_mt' in i] 
    buffers_mt=[(EXP_DATE+'/'+i) for i in file_names if 'buffers_mt' in i] 

    genemt_df = pd.read_excel(get_most_updated(genex_mt), index_col=0).dropna(axis=1, how='all')
    buffersmt_df = pd.read_excel(get_most_updated(buffers_mt), index_col=0).dropna(axis=1, how='all')

    def unpivot_df(df, type=None):
        df = df.stack().reset_index()
        df.columns=['experiment', 'component', 'concentration']
        if not type:
            df['type']=['DNA'] * len(df)
        else: 
            df['type']=['compound'] * len(df)
        return df

    input_data = pd.concat([unpivot_df(buffersmt_df, type=1),unpivot_df(genemt_df)])

    return pd.pivot_table(input_data, values='concentration', index=['experiment', 'type','component'])


def import_final_conc_data(EXP_DATE):
    
    data_file_names = os.listdir(EXP_DATE+'/data')
    metab_conc=[(EXP_DATE+'/data/'+i) for i in data_file_names if 'compiled' in i]
    metab_conc = get_most_updated(metab_conc)
    
    compiled_dfs = []
    for sheet_name, df in pd.read_excel(metab_conc, index_col=0, sheet_name=None).items():
        
        df = df.stack().reset_index()
        df.drop(['level_1'], axis=1, inplace=True)
        df.columns=['experiment','peak area']
        df['metabolite']=sheet_name
        compiled_dfs.append(df)
    all_data = pd.concat(compiled_dfs)
    all_data['count']=all_data.groupby(['metabolite', 'experiment']).transform('count')
    all_data['median'] = all_data.groupby(['experiment', 'metabolite'])['peak area'].transform('median')
    all_data['diff'] = abs(all_data['peak area']-all_data['median'])
    all_data['MAD'] = (all_data.groupby(['experiment', 'metabolite']).transform('sum')/all_data.groupby(['experiment', 'metabolite']).transform('count'))['diff']
    all_data['cutoff'] = (2*all_data['MAD']) + all_data['median']

    cleaned_df = all_data[(all_data['count']>2) & (all_data['peak area'] < all_data['cutoff']) | (all_data['count']<=2)]
    return cleaned_df[['experiment', 'metabolite', 'peak area']].groupby(['experiment', 'metabolite']).mean()

def peak_to_conc(met_data, std_curve_filepath): 
    std_curves = pd.read_csv(std_curve_filepath, index_col=0)# .dropna(axis=1, how='all')
    std_curves = std_curves.reset_index()
    std_curves.rename({'index': 'metabolite'}, axis=1, inplace=True)
    std_curves['metabolite'] = [i.lower() for i in std_curves['metabolite']]
    
    met_data = met_data.merge(std_curves, left_on='metabolite', right_on='metabolite')
    met_data['concentration'] = met_data['peak area']*met_data['slope']+ met_data['intercept']
    return met_data