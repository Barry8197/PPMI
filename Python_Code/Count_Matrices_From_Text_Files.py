#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os

dir = '/home/bryan/PPMI_IP/Python_Code/'
os.chdir(dir)

meta_df = pd.read_csv('./../Data/meta_data.11192021.csv')
IR3_path = '/data/PPMI/IR3_RNASeq/counts/'
IR3_filenames = os.listdir(IR3_path)

def random_sample(cl_event , d_status) :

    sample_list = meta_df[np.logical_and((meta_df['Genetic Status'] == cl_event),(meta_df['Disease Status'] == d_status))]['HudAlphaSampleName']

    return sample_list

def extract_counts(filename) : 
    df_tmp = pd.read_table(filename , skiprows=1)
    df_tmp = df_tmp.iloc[: , [0,-1]].set_index('Geneid')
    df_tmp.columns = [df_tmp.columns[0].split('.')[4]]
    
    return df_tmp

for status in meta_df['Disease Status'].drop_duplicates()[:-3] :
    path = dir + '/../Data_2/' + status.replace(' ','_')
    if not os.path.exists(path) :
        os.makedirs(path)
    event_list = meta_df[meta_df['Disease Status'] == status]['Genetic Status'].drop_duplicates()    
    for event in event_list[event_list.notna()] :
        path = dir + '/../Data_2/' + status.replace(' ','_') + '/' +event
        if not os.path.exists(path) :
            os.makedirs(path)
        df = None
        print('Creating Count Matrices for %s %s' % (status,event))
        for sample in random_sample(event , status) :
            for file in IR3_filenames :
                if file.split('.')[4] == sample :
                    filename = IR3_path + file
                    df_tmp = extract_counts(filename)
                    if df is None : 
                        df = df_tmp 
                    else :
                        df = df.join(df_tmp)
        df.to_csv((path + '/count_matrix.csv'))

