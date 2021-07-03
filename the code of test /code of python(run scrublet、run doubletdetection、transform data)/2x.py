#R python
import os
#os.chdir("Desktop/python test/")
#os.getcwd()
import scrublet as scr
#import solo as sl
#import doubletdetection
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

path = "counts2.csv"
raw_counts =load_csv(path)
raw_counts=np.delete(raw_counts,0, axis = 1)
raw_counts=raw_counts.astype(np.float)
raw_counts=raw_counts.transpose()   

scrub = scr.Scrublet(raw_counts)
doublet_scores=scrub.scrub_doublets()
x=scrub.calculate_doublet_scores()
np.savetxt('./scrublet2.txt',x)
#try:
#  doublet_scores, predicted_doublets = scrub.scrub_doublets()
#except:
#  np.savetxt('./scrublet2.txt',scrub.calculate_doublet_scores())
