import pandas as pd
import os

filename_IR2 = os.listdir('/data/PPMI/IR2_RNASeq/counts/')

print(filename_IR2[1].split('.'))
