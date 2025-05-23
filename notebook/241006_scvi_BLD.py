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

version = f'241006_scvi_BLD_doublet0.15_addribo_removelowribo_chem'

f_file_list = '/home/yy693/pi_hafler/ASAP/data/241010_cellranger_files_GEX.csv'
f_sample = '/home/yy693/pi_hafler/ASAP/data/ASAP BLD and CSF Single Cell Analysis 2024 September v1.xlsx'
f_cli = '/home/yy693/pi_hafler/ASAP/data/ASAPFinalClinicalData_20041002_2000029032KooMainPro-ASAPINearf.delspace.csv'

results_file_raw = f'../scanpy/{version}/{version}.raw.h5ad'
results_file_celltypist = f'../scanpy/{version}/{version}.celltypist.scrubletall.h5ad'
results_file_cxg_celltypist = f'../scanpy/{version}/{version}.cxg.celltypist.scrubletall.h5ad'


sc.settings.figdir = '../scanpy/{}/graph'.format(version)
sc.settings.cachedir = '../scanpy/{}/cache'.format(version)
# %config InlineBackend.figure_format = 'retina'

os.makedirs(sc.settings.figdir, exist_ok=True)
os.makedirs('../scanpy/{}'.format(version), exist_ok=True)


# threshold for the number of cells per sample
th_sample_cell = 50


df_sample = pd.read_excel(f_sample)
df_sample =  df_sample[['SampleID', 'Diagnosis', 'Gender', 'Age']]
df_sample.columns = ['DonorID', 'Diagnosis', 'Gender', 'Age']
df_sample['Diagnosis'] = df_sample['Diagnosis'].replace('PD+RBD', 'PDRBD')
# df_sample.head()

df_file = pd.read_csv(f_file_list)
df_file = df_file[df_file['SampleType'] == 'BLD']
# df_file.head()

# list_blacklist = ['/home/yy693/pi_hafler/ASAP/output/cellranger/20220502_lz438-1_5p/YPD036CSF_HHT_GEX/outs']
list_blacklist = ['/home/yy693/pi_hafler/ASAP/output/cellranger/20220610_lz438_5p/YPD040CSF_HHT_GEX/outs',
'/home/yy693/pi_hafler/ASAP/output/cellranger/20220610_lz438_5p/YPD040BLD_HHT_GEX/outs']
df_file = df_file[~df_file['path'].isin(list_blacklist)]

assert df_file.SampleID.duplicated().sum() == 0

list_cols = ['cellranger_sample_id', 'SampleID', 
       'SampleType', 'DonorID', 'Estimated Number of Cells',
       'Mean Reads per Cell', 'Median Genes per Cell', 'Number of Reads',
       'Valid Barcodes', 'Sequencing Saturation', 'Q30 Bases in Barcode',
       'Q30 Bases in RNA Read', 'Q30 Bases in UMI', 'Reads Mapped to Genome',
       'Reads Mapped Confidently to Genome',
       'Reads Mapped Confidently to Intergenic Regions',
       'Reads Mapped Confidently to Intronic Regions',
       'Reads Mapped Confidently to Exonic Regions',
       'Reads Mapped Confidently to Transcriptome',
       'Reads Mapped Antisense to Gene', 'Fraction Reads in Cells',
       'Total Genes Detected', 'Median UMI Counts per Cell',
       'Q30 Bases in RNA Read 2', 'chemistry']


dict_adata = {}
for pos,row in df_file.iterrows():
    a = sc.read_10x_h5(f"{row['path']}/filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()
    for col in list_cols:
        a.obs[col] = row[col]

    # sc.pp.scrublet(a, expected_doublet_rate=doublet_rate(a.shape[0], row['chemistry']))

    dict_adata[row['SampleID']] = a

adata = sc.concat(dict_adata, join='outer', index_unique='-')
adata.obs = pd.merge(adata.obs.reset_index(), df_sample, left_on='DonorID', right_on='DonorID', how='left').set_index('index')

df_cli = pd.read_csv(f_cli, index_col=0)
df_cli = df_cli.loc[:,~df_cli.columns.isin(adata.obs.columns)]
adata.obs = pd.merge(adata.obs.reset_index(), df_cli, how='left', left_on='DonorID', right_on='SubjectID').set_index('index')

sc.pp.scrublet(adata, sim_doublet_ratio=0.05, expected_doublet_rate=0.064) # estimated by cell number / sample number
adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)

adata.write(results_file_raw)

# adata = adata[adata.obs['predicted_doublet'] == False].copy()
# adata = adata[adata.obs['doublet_score'] < 0.18].copy()
adata = adata[adata.obs['doublet_score'] < 0.15].copy()

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", 'pct_counts_ribo'],
    jitter=0.4,
    multi_panel=True,
)

with plt.rc_context({"figure.figsize": (48, 5)}):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", 'pct_counts_ribo'],
        jitter=False,
        multi_panel=True,
        groupby='DonorID',
        rotation=90,
        save='sample.pdf'
    )

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=50)

print(f"Sample removed, # cells after QC <= {th_sample_cell}", 
list(adata.obs['SampleID'].value_counts()[(adata.obs['SampleID'].value_counts() < th_sample_cell)].index))
adata = adata[adata.obs['SampleID'].isin(
    adata.obs['SampleID'].value_counts()[(adata.obs['SampleID'].value_counts() >= th_sample_cell)].index)]

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

filter_genes = ['PPBP', 'GNG11', 'PF4', 'HBB']
threshold_exp = 1

# sc.pl.umap(adata,color=filter_genes, save='filtergenes.pdf')

for g in filter_genes:
    adata = adata[adata.raw[:, g].X <= threshold_exp]

adata = adata[adata.obs['pct_counts_ribo'] > 5].copy()

patterns_to_exclude = ["^IGKV", "^IGLV", "^IGHV", "^IGLC", 
                       "^TRAV", "^TRAJ", 
                       "^TRBV", "^TRBD", "^TRBJ", 
                       "^TRGV", "^TRGJ",
                       "^TRDV", "^TRDD", "^TRDJ", 
                       "^MT-", "^RPL", "^RPS"]

sc.pp.highly_variable_genes(adata, n_top_genes=4000, batch_key = 'SampleID')
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

sc.tl.umap(adata, min_dist=0.3, spread=2)
SCVI_CLUSTERS_KEY = "leiden_scVI"
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=1)

sc.pl.umap(adata,
           color=['CD3E', 'CD4', 'CD8A', 'CCR7', 'CD28', 'FOXP3', 'CD74', 'CD14', 
                       'MS4A1', 'IGHG1', 'HBB', 'PF4',
                 "n_genes_by_counts", "total_counts", "pct_counts_mt"],
           ncols=3, save='QC')


import celltypist
from celltypist import models

model = 'Immune_All_Low'
predictions = celltypist.annotate(adata, 
                                  model = f'/home/yy693/pi_hafler/ASAP/data/celltypist/{model}.pkl',
                                  majority_voting = True)

df_pred = predictions.predicted_labels
df_pred.columns = [f'predicted_labels_{model}', 
                   f'over_clustering_{model}', 
                    f'majority_voting_{model}']

adata.obs = pd.merge(adata.obs, df_pred, how='left', left_index=True, right_index=True)

model = 'COVID19_HumanChallenge_Blood'
predictions = celltypist.annotate(adata, 
                                  model = f'/home/yy693/pi_hafler/ASAP/data/celltypist/{model}.pkl',
                                  majority_voting = True)

df_pred = predictions.predicted_labels
df_pred.columns = [f'predicted_labels_{model}', 
                   f'over_clustering_{model}', 
                    f'majority_voting_{model}']

adata.obs = pd.merge(adata.obs, df_pred, how='left', left_index=True, right_index=True)

adata.write(results_file_celltypist)
adata_cxg = sc.AnnData(adata.raw.X)
adata_cxg.obs = adata.obs
adata_cxg.obsm = adata.obsm
adata_cxg.var.index = adata.raw.var.index
adata_cxg.write(results_file_cxg_celltypist)