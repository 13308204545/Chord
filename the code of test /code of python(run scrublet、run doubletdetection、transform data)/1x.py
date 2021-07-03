#R python
import os
#os.chdir("Desktop/python test/")
#os.getcwd()
#import scrublet as scr
import sys
sys.path.append("/Users/xiongkexu/Desktop/python test")
import doubletdetection
import numpy as np
import tarfile
import matplotlib.pyplot as plt
import pandas as pd

def load_csv(path):
    data_read = pd.read_csv(path)
    list = data_read.values.tolist()
    data = np.array(list)
    print(data.shape)
    # print(data)
    return data

#os.chdir("/Users/xiongkexu/Desktop/双胞项目/data/HTOdemux8/gradient/doubletrate/0.1T/")
path = "counts2.csv"
raw_counts =load_csv(path)
raw_counts=np.delete(raw_counts,0, axis = 1)
raw_counts=raw_counts.astype(np.float)
raw_counts=raw_counts.transpose()    

clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)
#doublets = clf.fit(raw_counts)
np.savetxt('doubletdetection2.txt',clf.voting_average_)
