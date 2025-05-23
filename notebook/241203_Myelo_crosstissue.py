#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
# import scrublet as scr
import scirpy as ir
# import muon as mu
#from vpolo.alevin import parser # to parse alevin output
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
import statsmodels.stats.multitest as multi

import re

seed = 0
np.random.seed(seed)

import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


# In[2]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis', transparent=False, frameon=False)  # low dpi (dots per inch) yields small inline figures

import matplotlib as mpl
# 2 lines below solved the facecolor problem.
# mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
sc.settings.autosave = True
sc.logging.print_header()

version = '241204_Myelo_crosstissue'
# input_table = '../data/231009_PBMC_HH514.csv'

results_file = '../scanpy/{}/res.h5ad'.format(version)
# results_file_mu = '../scanpy/{}/res.h5mu'.format(version)
results_file_cellxgene = '../scanpy/{}/res.cxg.h5ad'.format(version)

import os
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)

sc.settings.figdir = '../scanpy/{}/graph'.format(version)
sc.settings.cachedir = '../scanpy/{}/cache'.format(version)
# %config InlineBackend.figure_format = 'retina'

import os
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)
os.makedirs(sc.settings.figdir, exist_ok=True)


# In[3]:


dict_adata = {}


# ## PBMC Myelo (ours)

a = sc.read('/vast/palmer/pi/hafler/yy693/ASAP/scanpy/241021_BLD_integrate/241021_BLD_integrate.h5ad')

# d_clu = pd.read_csv('/vast/palmer/pi/hafler/yy693/ASAP/scanpy/241010_scvi_BLD_subcluster_myelo/241015.csv', index_col=0, skiprows=2)
# a.obs = pd.merge(a.obs, d_clu, left_index=True, right_index=True, how='left')
a = a[a.obs.cluster_L2.isin(['CD14 Mono', 'CD16 Mono', 'CD16 Mono C1', 'cDC2',
       'pDC', 'cDC1'])]

a.obs['cluster'] = a.obs['cluster_L2']
a.obs['sample'] = a.obs['SampleID']
a = sc.AnnData(a.layers['counts'], obs=a.obs[['cluster', 'sample']], var=a.var)

a.obs['organ'] = 'PBMC'
a.obs['project'] = 'PBMC ours'

dict_adata['PBMC_ours'] = a


# ## CSF Myelo (ours)
a = sc.read('/vast/palmer/pi/hafler/yy693/ASAP/scanpy/241027_CSF_integrate/241027_CSF_integrate.h5ad')
a = a[a.obs.cluster_L2.isin(['CD14 Mono', 'CD16 Mono', 'CSF Mac', 'cDC2',
       'pDC', 'cDC1', 'migDC'])]

a.obs['organ'] = 'CSF'
a.obs['project'] = 'CSF ours'
a.obs['cluster'] = a.obs['cluster_L2']
a.obs['sample'] = a.obs['SampleID']
a = sc.AnnData(a.layers['counts'], obs=a.obs[['cluster', 'organ', 'project', 'sample']], var=a.var)

dict_adata['CSF_ours'] = a

# ## Cross-tissue myeloid cell
a = sc.read('../data/crosstissue/CountAdded_PIP_myeloid_object_for_cellxgene.h5ad')

a.obs['cluster'] = a.obs['Manually_curated_celltype']
a.obs['project'] = 'Crosstissue'
a.obs['sample'] = a.obs['Donor']
a.obs['organ'] = a.obs['Organ']
a = sc.AnnData(a.layers['counts'], obs=a.obs[['cluster', 'project', 'sample', 'organ']], var=a.var)

dict_adata['Crosstissue'] = a


# ## Microglia Olah et al., 2020
d = pd.read_csv('../data/Microglia_Olah_et_al_2020/SupplementaryData14.csv', index_col=0)
d_obs = pd.read_csv('../data/Microglia_Olah_et_al_2020/SupplementaryData15.csv')
d_barcode = pd.read_csv('../data/Microglia_Olah_et_al_2020/cell_barcodes.csv', index_col=0)
d_obs = pd.merge(d_barcode, d_obs, how='left', left_on='cellid', right_on='sample_id').set_index('cellid')

a = sc.AnnData(d.T, obs=d_obs)

a.obs['organ'] = 'Brain'
a.obs['project'] = 'MicroOlah'
a.obs['sample'] = a.obs['donor']
a.obs['cluster'] = a.obs['cluster_label']

dict_adata['Microglia_Olah'] = a


# ## merge, downstream
adata = sc.concat(dict_adata, index_unique="-")
adata.obs['cluster'] = adata.obs.project.astype(str) + '-' + adata.obs.cluster.astype(str)


sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_cells(adata, min_genes=1000)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], 
    percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], groupby='sample', rotation=90,
             jitter=0.4, multi_panel=True, save='QC.pdf')

adata = adata[adata.obs.n_genes_by_counts < 4000, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]

# filter out BCR,TCR,MT,Ribosomal genes
patterns_to_exclude = ["^IGKV", "^IGLV", "^IGHV", "^IGLC", 
                       "^TRAV", "^TRAJ", 
                       "^TRBV", "^TRBD", "^TRBJ", 
                       "^TRGV", "^TRGJ",
                       "^TRDV", "^TRDD", "^TRDJ", 
                       "^MT-", "^RPL", "^RPS"]

mask = ~adata.var.index.str.contains('|'.join(patterns_to_exclude), flags=re.IGNORECASE)

adata = adata[:, mask]

# adata.write(results_file_human_count)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    # subset=True,
    subset=False,
    layer="counts",
    flavor="seurat_v3",
    span=1,
    # flavor="seurat",
    batch_key="sample",
)

sc.pl.highly_variable_genes(adata, log=True)
# adata = adata[:, adata.var.highly_variable]

cell_cycle_genes = [x.strip() for x in open('../data/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.pp.scale(adata, max_value=10)
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'S_score', 'G2M_score'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

sce.pp.harmony_integrate(adata, 'sample')
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40, use_rep='X_pca_harmony')
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=4,
    values_to_plot="logfoldchanges",
    min_logfoldchange=2,
    vmax=7,
    vmin=-7,
    cmap="bwr",
    save='leiden_QC'
)
adata.write(results_file)

adata_cxg = sc.AnnData(adata.raw.X)
adata_cxg.obs = adata.obs
adata_cxg.var = adata.var
adata_cxg.obsm = adata.obsm
adata_cxg.write(results_file_cellxgene)

sc.pl.umap(adata, color=['CD3E', 'CD4', 'FOXP3', 'CCR7', 'CD8A', 'NKG7', 'MS4A1', 'CD74', 'HBB', 'PF4', # platelet
                         'GNLY', 'ITGAX', 'LAMP3', 'IRF8', 'XCR1', 'FCER1A'], save='markers')

sc.pl.umap(adata, color=['cluster', 'organ', 'project'], save='QC')
