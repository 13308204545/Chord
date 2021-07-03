#R python
import os
#os.chdir("Desktop/python test/")
#os.getcwd()
import scanpy as sc
import pandas as pd
adata = sc.read_loom("./seu2.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata.write('seu2.h5ad')

