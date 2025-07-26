import cellrank as cr
import pandas as pd

adata=cr.read('/pathto/hAG_cnts.h5ad')

obs=adata.obs
obsm=adata.obsm
obs.to_csv('obs.csv')
mtx=adata.X
mtx=pd.DataFrame(mtx)
mtx.to_csv('mtx.csv')

genes=adata.var_names
genes=pd.DataFrame(genes)
genes.to_csv('genes.csv')

cells=adata.obs_names
cells=pd.DataFrame(cells)
cells.to_csv('cells.csv')
