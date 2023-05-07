#! /usr/bin/env python
# -*- coding:utf-8 -*-
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from sklearn.metrics import adjusted_rand_score


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    for res in sorted(list(np.arange(0.02, 2, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        if count_unique_leiden == fixed_clus_count:
            return res


def mclust_R(x, n_clusters, model='EEE', random_seed=2020):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    # os.environ['R_HOME'] = '/GPUFS/sysu_ydyang_10/.conda/envs/r-base/lib/R'
    # os.environ['R_USER'] = '/GPUFS/sysu_ydyang_10/.conda/envs/r-base/lib/python3.9/site-packages/rpy2'

    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(x), n_clusters, model)
    mclust_res = np.array(res[-2]).astype(int) - 1
    
    return mclust_res


def eval_mclust_ari(labels, z, n_clusters):
    raw_preds = mclust_R(z, n_clusters)
    
    preds = raw_preds[labels != -1]
    labels = labels[labels != -1]
    
    return  raw_preds
