import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import networkx as nx
import scglue
from itertools import chain


def create_anndata_from_file(file_path):
    """
    Function to read a tab-separated values file and create an AnnData object
    with spatial coordinates and data.

    Parameters:
    - file_path: str, the path to the input file.

    Returns:
    - adata: AnnData object containing the data and spatial coordinates.
    """
    # Read the file
    data_df = pd.read_csv(file_path, sep='\t')  # Assuming the data is tab-separated

    # Get the first and second rows of data (excluding the index)
    first_row = data_df.iloc[:, 2]  # Assuming the first row's index is 0, skipping the first column
    second_row = data_df.iloc[:, 3]  # Same as above, assuming the second row's index is 1

    # Convert the first and second rows to spatial coordinates
    X_spatial = np.stack((first_row.to_numpy(), second_row.to_numpy()), axis=-1)

    # Create an index
    index = first_row.astype(str) + 'x' + second_row.astype(str)
    data_df.index = index

    # Drop the first four columns as they have been used for indexing
    data_df = data_df.drop(data_df.columns[0:4], axis=1)
    # data_df = data_df.applymap(lambda x: float(x / 10))
    data_df = data_df.astype(int)
    # Create an AnnData object
    adata = ad.AnnData(X=data_df)

    # Add the DataFrame's row index as the AnnData's obs (observations) index
    adata.obs = pd.DataFrame(index=data_df.index)

    # Add the DataFrame's column index as the AnnData's var (variables) index
    adata.var = pd.DataFrame(index=data_df.columns)

    # Add spatial coordinates
    adata.obsm['X_spatial'] = X_spatial

    return adata

def create_sma_anndata_from_file(file_path):
    """
    Function to read a tab-separated values file and create an AnnData object
    with spatial coordinates and data.

    Parameters:
    - file_path: str, the path to the input file.

    Returns:
    - adata: AnnData object containing the data and spatial coordinates.
    """
    # Read the file
    data_df = pd.read_csv(file_path, sep=',')  # Assuming the data is tab-separated

    # Get the first and second rows of data (excluding the index)
    first_row = data_df.iloc[:, 0]  # Assuming the first row's index is 0, skipping the first column
    second_row = data_df.iloc[:, 1]  # Same as above, assuming the second row's index is 1

    # Convert the first and second rows to spatial coordinates
    X_spatial = np.stack((first_row.to_numpy(), second_row.to_numpy()), axis=-1)

    # Create an index
    index = first_row.astype(str) + 'x' + second_row.astype(str)
    data_df.index = index

    # Drop the first four columns as they have been used for indexing
    data_df = data_df.drop(data_df.columns[0:4], axis=1)
    # data_df = data_df.applymap(lambda x: float(x / 10))
    data_df = data_df.astype(int)
    # Create an AnnData object
    adata = ad.AnnData(X=data_df)

    # Add the DataFrame's row index as the AnnData's obs (observations) index
    adata.obs = pd.DataFrame(index=data_df.index)

    # Add the DataFrame's column index as the AnnData's var (variables) index
    adata.var = pd.DataFrame(index=data_df.columns)

    # Add spatial coordinates
    adata.obsm['X_spatial'] = X_spatial

    return adata
def process_sma_mt_data(file_path):
    """
    Process single cell RNA sequencing data.

    Parameters:
    file_path (str): The path to the file containing the data.

    Returns:
    adata (anndata.AnnData): The processed AnnData object.
    """
    adata = create_sma_anndata_from_file(file_path)
    adata.var_names_make_unique()

    adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)

    sc.pp.normalize_total(adata, inplace=True)
    #sc.pp.log1p(adata)

    sc.pp.scale(adata)

    sc.pp.pca(adata, n_comps=100)

    sc.pp.neighbors(adata)

    sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=0.4,
        random_state=0,
        n_iterations=2,
        directed=False,
    )

    return adata
def process_mt_data(file_path,scale=True):
    """
    Process single cell RNA sequencing data.

    Parameters:
    file_path (str): The path to the file containing the data.

    Returns:
    adata (anndata.AnnData): The processed AnnData object.
    """

    adata = create_anndata_from_file(file_path)
    adata.var_names_make_unique()

    adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=2000)

    sc.pp.normalize_total(adata, inplace=True)
    #sc.pp.log1p(adata)
    if scale:
        sc.pp.scale(adata)

    sc.pp.pca(adata, n_comps=100)

    sc.pp.neighbors(adata)

    sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=0.4,
        random_state=0,
        n_iterations=2,
        directed=False,
    )

    return adata



def process_st_data(visium_path,scale=True):
    """
    Process Visium spatial transcriptomics data.

    Parameters:
    visium_path (str): The path to the directory containing the Visium output files.

    Returns:
    adata (anndata.AnnData): The processed AnnData object.
    """
    adata = sc.read_visium(visium_path)
    adata.var_names_make_unique()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata.layers["counts"] = adata.X.copy()

    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.histplot(
        adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
        kde=False,
        bins=40,
        ax=axs[1],
    )
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.histplot(
        adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
        kde=False,
        bins=60,
        ax=axs[3],
    )
    plt.show()


    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)

    if scale:
        sc.pp.scale(adata)

    sc.pp.pca(adata, n_comps=100)

    sc.pp.neighbors(adata)

    sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=0.5,
        random_state=0,
        n_iterations=2,
        directed=True,
    )

    return adata

