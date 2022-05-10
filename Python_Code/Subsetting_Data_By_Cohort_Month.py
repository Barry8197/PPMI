#!/usr/bin/env python
# coding: utf-8

# In[128]:


import pandas as pd
import numpy as np
import os

os.chdir('/home/bryan/PPMI_IP/Python_Code/')

meta_df = pd.read_csv('./../Data/meta_data.11192021.csv')
IR2_filenames = np.genfromtxt('./../Data/filenamesIR2.txt' , dtype = 'str')
IR3_filenames = np.genfromtxt('./../Data/filenamesIR3.txt' , dtype = 'str')

def random_sample(cl_event , d_status , num) :
    
    sample_list = meta_df[np.logical_and((meta_df['Clinical Event'] == cl_event),(meta_df['Disease Status'] == d_status))]['HudAlphaSampleName']
    
    return sample_list.sample(n = num , replace = True)

for status in meta_df['Disease Status'].drop_duplicates()[:-3] :
    path = './../Data_2/' + status.replace(' ','_')
    if not os.path.exists(path) :
        os.makedirs(path)
    for event in meta_df['Genetic Status'].drop_duplicates()[:-1] :
        path = './../Data_2/' + status.replace(' ','_') + '/' +event
        if not os.path.exists(path) :
            os.makedirs(path)
        for sample in random_sample(event , status , 10) : 
            for file in IR3_filenames : 
                if file.split('.')[4] == sample :
                    os.system('cp /data/PPMI/IR3_RNASeq/counts/' + file + ' ' + path + '/')
