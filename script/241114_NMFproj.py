import os
import glob
import pandas as pd
import scanpy as sc
from NMFproj import *

# Set the working directory
os.chdir('/home/yy693/pi_hafler/ASAP')

# Parameters setting
dir_output = "./output/241114_NMFproj"

fixed_W = pd.read_csv('/home/yy693/pi_hafler/yyoshiaki-git/NMFprojection/data/NMF.W.CD4T.csv.gz', index_col=0)

# Get the list of samples in the directory
list_h5ad = glob.glob('output/CD4T_screfmapping/*/*/*.h5ad')

# NMF process
for f_h5ad in list_h5ad:
    # print(sample)
    print('Processing', f_h5ad)
    # prefix = f"./output/{project_name}/{project}/{sample}"
    
    adata = sc.read(f_h5ad)

    X_norm, X_trunc, df_H, fixed_W_trunc = NMFproj(adata.to_df().T, fixed_W, return_truncated=True, normalized=False)
    df_H.to_csv('{}_projection.csv'.format(f_h5ad.replace('.h5ad', '')))
    # if os.path.exists(f"{prefix}_CD4T_AssayData.h5ad"):
    # cmd = [
    #     "NMFproj",
    #     f"{prefix}_CD4T_X.tsv",
    #     "--outputprefix", f"{prefix}_CD4T",
    #     "/home/yy693/pi_hafler/yyoshiaki-git/NMFprojection/data/NMF.W.CD4T.csv.gz"
    # ]
    
    # subprocess.run(cmd)
    
