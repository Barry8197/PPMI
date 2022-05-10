import pandas as pd
import os

df = pd.read_csv('./../Data/meta_data.11192021.csv')
filename_IR2 = os.listdir('/data/PPMI/IR2_RNASeq/counts/')
filename_IR3 = os.listdir('/data/PPMI/IR3_RNASeq/counts/')

both = 0
count_IR2 = 0
count_IR3 = 0
neither = 0
for pat in df['HudAlphaSampleName'] : 
    check_IR2 = False 
    check_IR3 = False

    for filename in filename_IR3 :
        if filename.split('.')[4] == pat :
            check_IR3 = True

    for filename in filename_IR2 : 
        if filename.split('.')[4] == pat :
            check_IR2 = True

    if (check_IR2 == True) and (check_IR3 == True) :
        both += 1
    elif check_IR2 == True :
        count_IR2 += 1
    elif check_IR3 == True :
        count_IR3 += 1
    else :
        neither += 1


print('count IR3 %f' % count_IR3)
print('count IR2 %f' % count_IR2)
print('count both %f' % both)
print('count neither %f' % neither)
