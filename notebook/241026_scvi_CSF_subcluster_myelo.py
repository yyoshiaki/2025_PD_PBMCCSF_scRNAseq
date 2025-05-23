#!/usr/bin/env python
# coding: utf-8

# @Author: yyasumizu
# @Date: 2024-09-29
# @Last Modified time: 2024-09-29

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import matplotlib.pyplot as plt
import seaborn as sns
import scvi

import statsmodels.api as sm
import statsmodels.stats.multitest as multi

import re
import os

seed = 0
np.random.seed(seed)

import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis', transparent=False, frameon=False)  # low dpi (dots per inch) yields small inline figures

import matplotlib as mpl
# 2 lines below solved the facecolor problem.
# mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
sc.settings.autosave = True
sc.logging.print_header()

version = f'241026_scvi_CSF_subcluster_myelo'

f_adata = '/home/yy693/pi_hafler/ASAP/scanpy/241026_scvi_CSF_doublet0.15_addribo_removelowribo_chem_l1d10/241026_scvi_CSF_doublet0.15_addribo_removelowribo_chem_l1d10.celltypist.scrubletall.h5ad'
# f_cl1 = '/vast/palmer/home.mccleary/yy693/cellxgene-data/yy693/Yale/ASAP/241010_scvi_CSF_doublet0.15_addribo_removelowribo.cxg.celltypist.scrubletall_annotations/241015.csv'
f_cl1='/vast/palmer/pi/hafler/yy693/ASAP/scanpy/241026_scvi_CSF_doublet0.15_addribo_removelowribo_chem_l1d10/241026.csv'
cluster = 'Myeloid'

results_file = f'../scanpy/{version}/{version}.h5ad'
results_file_cxg = f'../scanpy/{version}/{version}.cxg.h5ad'

sc.settings.figdir = '../scanpy/{}/graph'.format(version)
sc.settings.cachedir = '../scanpy/{}/cache'.format(version)
# %config InlineBackend.figure_format = 'retina'

os.makedirs(sc.settings.figdir, exist_ok=True)
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)


adata = sc.read(f_adata)
df_clu = pd.read_csv(f_cl1, skiprows=2, index_col=0)
adata.obs = pd.merge(adata.obs, df_clu, left_index=True, right_index=True, how='left')
adata = adata[adata.obs['cluster_L1'] == cluster].copy()

print(adata)

patterns_to_exclude = ["^IGKV", "^IGLV", "^IGHV", "^IGLC", 
                       "^TRAV", "^TRAJ", 
                       "^TRBV", "^TRBD", "^TRBJ", 
                       "^TRGV", "^TRGJ",
                       "^TRDV", "^TRDD", "^TRDJ", 
                       "^MT-", "^RPL", "^RPS"]

sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key = 'SampleID')
adata.var.highly_variable = \
    ((~adata.var.index.str.contains('|'.join(patterns_to_exclude), flags=re.IGNORECASE)) &
                 (adata.var.highly_variable))


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["SampleID", 'chemistry'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
)

model = scvi.model.SCVI(adata)

print(model)

model.train()

model_dir = os.path.join(f'../scanpy/{version}', "scvi_model")
model.save(model_dir, overwrite=True)

model = scvi.model.SCVI.load(model_dir, adata=adata)

SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
print(latent.shape)

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)

# sc.tl.umap(adata, min_dist=0.3, spread=2)
sc.tl.umap(adata, min_dist=0.5, spread=1)
SCVI_CLUSTERS_KEY = "leiden_scVI"
# sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=1)
for res in [0.7, 1, 1.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )

adata.write(results_file)
adata_cxg = sc.AnnData(adata.raw.X)
adata_cxg.obs = adata.obs
adata_cxg.obsm = adata.obsm
adata_cxg.var.index = adata.raw.var.index
adata_cxg.write(results_file_cxg)

sc.pl.umap(adata,
           color=['CD3E', 'CD4', 'CD8A', 'CCR7', 'CD28', 'FOXP3', 'CD74',
            'CD14', 'FCGR3A', 'LAMP3',
                       'MS4A1', 'IGHG1', 'HBB', 'PF4',
                 "n_genes_by_counts", "total_counts", "pct_counts_mt"],
           ncols=3, save='QC')

for res in [0.7, 1, 1.5, 2.0]:
    sc.tl.rank_genes_groups(adata, groupby=f"leiden_res_{res:4.2f}", method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(
        adata, groupby=f"leiden_res_{res:4.2f}", standard_scale="var", n_genes=5,
        save=f"leiden_res_{res:4.2f}"
    )
    sc.pl.umap(adata, color=f"leiden_res_{res:4.2f}", legend_loc="on data", save=f"leiden_res_{res:4.2f}")