import numpy as np
from sklearn.kernel_approximation import Nystroem
from sklearn.metrics.pairwise import rbf_kernel
import torch
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import gzip
from torch.utils.data import Subset

def guassian_kernel(source, target, kernel_mul=2.0, kernel_num=5, fix_sigma=None):
    n_samples = int(len(source))+int(len(target))
    total = torch.cat([source, target], dim=0)

    total0 = total.unsqueeze(0).expand(int(total.size(0)), \
                                       int(total.size(0)), \
                                       int(total.size(1)))
    total1 = total.unsqueeze(1).expand(int(total.size(0)), \
                                       int(total.size(0)), \
                                       int(total.size(1)))
    L2_distance = ((total0-total1)**2).sum(2)


    if fix_sigma:
        bandwidth = fix_sigma
    else:
        bandwidth = torch.sum(L2_distance.data) / (n_samples**2-n_samples)
    bandwidth /= kernel_mul ** (kernel_num // 2)
    bandwidth_list = [bandwidth * (kernel_mul**i) for i in range(kernel_num)]

    kernel_val = [torch.exp(-L2_distance / bandwidth_temp) for \
                  bandwidth_temp in bandwidth_list]

    return sum(kernel_val)



def mmd(source, target, kernel_mul=2.0, kernel_num=5, fix_sigma=None):
    batch_size_source = int(source.size()[0])
    batch_size_target = int(target.size()[0])

    if batch_size_source == 0 or batch_size_target == 0:
        return torch.tensor(0.0)
    else:

        kernels = guassian_kernel(source, target, kernel_mul=kernel_mul, kernel_num=kernel_num, fix_sigma=fix_sigma)

        XX = kernels[:batch_size_source, :batch_size_source]  # Source<->Source
        YY = kernels[batch_size_source:, batch_size_source:]  # Target<->Target
        XY = kernels[:batch_size_source, batch_size_source:]  # Source<->Target
        YX = kernels[batch_size_source:, :batch_size_source]  # Target<->Source

        M = torch.zeros((batch_size_source + batch_size_target, batch_size_source + batch_size_target))
        M[:batch_size_source, :batch_size_source] = 1.0 / (batch_size_source * batch_size_source) if batch_size_source > 0 else 0
        M[batch_size_source:, batch_size_source:] = 1.0 / (batch_size_target * batch_size_target) if batch_size_target > 0 else 0
        M[:batch_size_source, batch_size_source:] = -1.0 / (batch_size_source * batch_size_target) if batch_size_source > 0 and batch_size_target > 0 else 0
        M[batch_size_source:, :batch_size_source] = -1.0 / (batch_size_target * batch_size_source) if batch_size_source > 0 and batch_size_target > 0 else 0

        loss = torch.sum(M * kernels)
        return loss



def calculate_mmd_matrix(domain_rna, domain_meta):
    domain_rna_umap = pd.DataFrame(domain_rna.obsm['X_umap'], index=domain_rna.obs.index)
    domain_meta_umap = pd.DataFrame(domain_meta.obsm['X_umap'], index=domain_meta.obs.index)

    rna_groups = domain_rna.obs['bayes'].unique()
    meta_groups = domain_meta.obs['louvain'].unique()
    rna_indices = {group: domain_rna.obs.index[domain_rna.obs['bayes'] == group] for group in rna_groups}
    meta_indices = {group: domain_meta.obs.index[domain_meta.obs['louvain'] == group] for group in meta_groups}

    num_rna_groups = len(rna_groups)
    num_meta_groups = len(meta_groups)

    matrix = [[0] * num_meta_groups for _ in range(num_rna_groups)]

    for rna_group in rna_groups:
        for meta_group in meta_groups:
            rna_data = domain_rna_umap.loc[rna_indices[rna_group], :].values
            meta_data = domain_meta_umap.loc[meta_indices[meta_group], :].values
            mmddt = mmd(torch.tensor(rna_data), torch.tensor(meta_data))
            rna_idx = list(rna_groups).index(rna_group)
            meta_idx = list(meta_groups).index(meta_group)
            matrix[rna_idx][meta_idx] = mmddt.item()

    df_matrix = pd.DataFrame(matrix, index=[f'RNA Group {group}' for group in rna_groups], columns=[f'Meta Group {group}' for group in meta_groups])
    return df_matrix

