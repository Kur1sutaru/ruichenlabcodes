#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scvelo as scv
import numpy as np
import os
import pandas as pd
import anndata
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams


# 1) To read the .h5ad obj with models already calculated - obj_adata_kns_graph_wumap_velGraph_2_SSmodel.h5ad
#===============================
# user's parameters
#===============================
finh5ad="/storage/chentemp/u247700/kns/velocity/degsvelo/obj_adata_kns_graph_wumap_velGraph_2_SSmodel.h5ad"
prefix="kns"
output_path="/storage/chentemp/u247700/kns/velocity/degsvelo/"
#===============================


# create out_dir
os.makedirs(output_path, exist_ok=True)
os.chdir(output_path)


adata = anndata.read_h5ad("/storage/chentemp/u247700/kns/velocity/degsvelo/obj_adata_kns_graph_wumap_velGraph_2_SSmodel.h5ad")

print("input was loaded")



# 2) Using dynamical model to determine the transcritpion rate, splicing rate, and degradation rate

scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
pl.savefig('transcritpionsplicingdegradationrate.png')


scv.get_df(adata, 'fit*', dropna=True).head()
adata.to_csv('resultstfspldeg.csv')

# 3) degs ranked
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv('resultsdegs.csv')


print("process was done")
