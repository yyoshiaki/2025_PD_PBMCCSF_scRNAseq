# %%
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import decoupler as dc
# import scrublet as scr
# import muon as mu
#from vpolo.alevin import parser # to parse alevin output
import matplotlib.pyplot as plt
import seaborn as sns

import scvi
import torch

import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from adjustText import adjust_text
import re

seed = 0
np.random.seed(seed)


import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# %%
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

# %%
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis', transparent=False, frameon=False)  # low dpi (dots per inch) yields small inline figures

import matplotlib as mpl
# 2 lines below solved the facecolor problem.
# mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.facecolor'] = (1,1,1,1)
sc.settings.autosave = True
sc.logging.print_header()

version = '250417_MS_scanvi'
# input_table = '../data/231009_PBMC_HH514.csv'

results_file_query = '../scanpy/{}/merged.CSF_MSandHA_integrated_60704.scanvi.h5ad'.format(version)

save_dir = '../scanpy/{}/scvi'.format(version)

import os
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)

sc.settings.figdir = '../scanpy/{}/graph'.format(version)
sc.settings.cachedir = '../scanpy/{}/cache'.format(version)
# %config InlineBackend.figure_format = 'retina'

import os
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)
os.makedirs(sc.settings.figdir, exist_ok=True)

# %%
# adata_query = sc.read_h5ad('/vast/palmer/pi/hafler/yy693/Ocrevus/data/merged.CSF_MSandHA_integrated_60704.h5ad')

adata_query = sc.read('../data/CSF_MS.allinhousedata.raw.h5ad')
adata_query.layers['counts'] = adata_query.X
adata_query.obs['disease'] = adata_query.obs['disease'].replace({'MS':'MS', 'HC':'Healthy'})
adata_query.obs.patient = adata_query.obs.patient.astype(str)

# %%
# adata_ref = sc.read('/home/yy693/pi_hafler/ASAP/scanpy/241028_CSF_downstream/res.h5ad')
adata_ref = sc.read('/home/yy693/pi_hafler/ASAP/scanpy/250407_CSF_downstream/res.h5ad')

# retain only overlapped genes
genes = list(set(adata_ref.var.index) & set(adata_query.var.index))
adata_ref = adata_ref[:, genes].copy()
adata_query = adata_query[:, genes].copy()

# calcurate percent ribo
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata_query.var["mt"] = adata_query.var_names.str.startswith("MT-")
# ribosomal genes
adata_query.var["ribo"] = adata_query.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata_query.var["hb"] = adata_query.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    adata_query,
    qc_vars=["mt", "ribo"],
    percent_top=None,
    log1p=False,
    layer="counts",
    inplace=True,
)

adata_query.obs['SampleID'] = adata_query.obs['patient']
# adata_query.obs['chemistry'] = 'SC5P-R2'


adata_query.X = adata_query.layers['counts'].copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata_query, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata_query)

# %%
# scvi.model.SCVI.setup_anndata(
#     adata_ref,
#     layer="counts",
#     categorical_covariate_keys=["SampleID", 'chemistry'],
#     continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
# )
scvi.model.SCVI.setup_anndata(
    adata_ref,
    layer="counts",
    batch_key="SampleID",
    # continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
)



# %%
scvi_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)
scvi_ref.train()

# %%
SCVI_LATENT_KEY = "X_scVI"

adata_ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)

# %%
sc.pl.umap(
    adata_ref,
    color=["SampleID", "cluster_L2"],
    frameon=False,
    ncols=1,
)

# %%
scvi_ref_path = os.path.join(save_dir, "adata_scvi_ref")
scvi_ref.save(scvi_ref_path, overwrite=True)


# %%
SCANVI_LABELS_KEY = "cluster_L2_scanvi"
adata_ref.obs[SCANVI_LABELS_KEY] = adata_ref.obs["cluster_L2"].values

# %%
# unlabeled category does not exist in adata.obs[labels_key]
# so all cells are treated as labeled
scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_ref,
    unlabeled_category="Unknown",
    labels_key=SCANVI_LABELS_KEY,
)

# %%
scanvi_ref.train(max_epochs=20, n_samples_per_label=100)

# %%
SCANVI_LATENT_KEY = "X_scANVI"

adata_ref.obsm[SCANVI_LATENT_KEY] = scanvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)

# %%
sc.pl.umap(
    adata_ref,
    color=["SampleID", "cluster_L2"],
    frameon=False,
    ncols=1,
)

# %%
scanvi_ref_path = os.path.join(save_dir, "adata_scanvi_ref")
scanvi_ref.save(scanvi_ref_path, overwrite=True)

# %%
# again a no-op in this tutorial, but good practice to use
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_ref_path)

# %%
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, scanvi_ref_path)

# %%
scanvi_query.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10,
)

# %%
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()



sc.pp.neighbors(adata_query, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)

# %%
adata_query.write(results_file_query)

# %%
sc.pl.umap(adata_query, color=['MS4A1', 'TCL1A', 'BANK1', 'CD27', 'MZB1', 
'CD3E', 'CD4', 'SELL', 'FAS', 'CD28', 'CCL5', 'FOXP3', 'CD8A', 'MKI67',
'NKG7', 'FCGR3A', 'NCAM1', 'CST3', 'CD14', 'CDKN1C',
'C1QC', 'LAMP3', 'ITM2C', 'CLEC9A', 'FCER1A', 'FN1', 'predictions_scanvi'], save='query_l2markers')

sc.pl.dotplot(adata_query, var_names=['MS4A1', 'TCL1A', 'BANK1', 'CD27', 'MZB1', 
'CD3E', 'CD4', 'SELL', 'FAS', 'CD28', 'CCL5', 'FOXP3', 'CD8A', 'MKI67',
'NKG7', 'FCGR3A', 'NCAM1', 'CST3', 'CD14', 'CDKN1C',
'C1QC', 'LAMP3', 'ITM2C', 'CLEC9A', 'FCER1A', 'FN1'], 
    groupby='predictions_scanvi', cmap='Purples', save='l2markers')

# # %%
# df = adata_query.obs.groupby(["cluster_label", SCANVI_PREDICTIONS_KEY]).size().unstack(fill_value=0)
# norm_df = df / df.sum(axis=0)

# plt.figure(figsize=(8, 8))
# _ = plt.pcolor(norm_df)
# _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
# _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
# plt.xlabel("Predicted")
# plt.ylabel("Observed")
# plt.savefig('../scanpy/{}/graph/scanvi_query_predicted_observed.png'.format(version), dpi=300)