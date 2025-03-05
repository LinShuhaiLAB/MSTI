import pandas as pd
import re
import anndata as ad
from collections import defaultdict
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.sparse import csr_matrix
from tqdm import tqdm


def extract_molecule_meta(equation):

    molecule_ids = re.findall(r'\bC\d+\b', equation)

    unique_molecule_ids = list(dict.fromkeys(molecule_ids))

    return ', '.join(unique_molecule_ids)


def extract_molecule_rna(enzyme_series):

    return enzyme_series.str.replace('///', ', ').str.strip()



def map_kegg_to_entries(dataframe):

    reverse_mapping = defaultdict(list)


    for _, row in dataframe.iterrows():
        entry = row['ENTRY']
        for kegg in row['KEGG']:
            reverse_mapping[kegg].append(entry)


    df = pd.DataFrame([(kegg, entries) for kegg, entries in reverse_mapping.items()], columns=['KEGG', 'ENTRIES'])

    return df


def process_kegg_reaction(kegg_reaction, sample_mt, kegg_rna):



    kegg_reaction['KEGG'] = kegg_reaction['EQUATION'].apply(extract_molecule_meta)
    kegg_reaction['EC'] = extract_molecule_rna(kegg_reaction['ENZYME'])


    kegg_reaction = kegg_reaction[['ENTRY', 'KEGG', 'EC']]


    kegg_reaction['KEGG'] = kegg_reaction['KEGG'].str.split(',')
    kegg_reaction = kegg_reaction.explode('KEGG', ignore_index=True)


    kegg_reaction['EC'] = kegg_reaction['EC'].str.split(',')
    kegg_reaction = kegg_reaction.explode('EC', ignore_index=True)


    kegg_reaction = kegg_reaction.merge(sample_mt[['KEGG', 'Precursor (Da)']], on='KEGG', how='left')
    kegg_reaction = kegg_reaction.merge(kegg_rna[['EC', 'gene_name']], on='EC', how='left')


    kegg_reaction = kegg_reaction.dropna(axis=0)

    return kegg_reaction


def create_graph_from_data(rna, meta, kegg_reaction):

    G = nx.Graph()


    rna_nodes_to_add = rna.var.index.tolist()
    meta_nodes_to_add = meta.var.index.tolist()
    nodes_to_add = rna_nodes_to_add + meta_nodes_to_add
    G.add_nodes_from(nodes_to_add)

    if kegg_reaction is not None:
        for index, row in kegg_reaction.iterrows():
            precursor = str(row['Precursor (Da)']).strip()
            gene_name = str(row['gene_name']).strip()

            if G.has_node(precursor) and G.has_node(gene_name):
                G.add_edge(precursor, gene_name, weight=1, sign=1)

    for node in G.nodes:
        G.add_edge(node, node, weight=1, sign=1)

    return G


def create_basic_graph(rna, meta):

    G = nx.Graph()


    rna_nodes_to_add = rna.var.index.tolist()
    meta_nodes_to_add = meta.var.index.tolist()
    nodes_to_add = rna_nodes_to_add + meta_nodes_to_add
    G.add_nodes_from(nodes_to_add)

    for node in G.nodes:
        G.add_edge(node, node, weight=1, sign=1)

    return G



def add_edges_with_weight_sign(G, mtadata,target_meta_list, stadata,target_rna_list,top_n=3):
    if not isinstance(G, nx.Graph):
        raise ValueError("G must be a networkx Graph object")
    for target_meta, target_rna in zip(target_meta_list, target_rna_list):

        meta_list = get_top_similar_genes(mtadata, target_meta, top_n)
        rna_list = get_top_similar_genes(stadata, target_rna, top_n)

        for meta in meta_list:
            for rna in rna_list:
                if not G.has_edge(meta, rna):
                    G.add_edge(meta, rna, weight=1, sign=1)


    return G
def add_edges(G, meta_list, rna_list):


    for meta in meta_list:
        for rna in rna_list:
            if not G.has_edge(meta, rna):
                G.add_edge(meta, rna, weight=1, sign=1)


    return G

def get_top_similar_genes(adata, target_gene, top_n):





    target_vector = adata[:, target_gene].X.toarray().flatten()

    target_mean = np.mean(target_vector)
    target_std = np.std(target_vector)

    similarity_scores = {}

    for gene in tqdm(adata.var.index, desc="Processing"):
        if gene == target_gene:
            continue

        gene_vector = adata[:, gene].X.toarray().flatten()

        if np.std(gene_vector) == 0:
            continue

        correlation, _ = pearsonr(target_vector, gene_vector)
        similarity_scores[gene] = correlation

    top_genes = sorted(similarity_scores, key=similarity_scores.get, reverse=True)[:top_n]

    # print(f"Top {top_n} genes similar to {target_gene}:")
    # for i, gene in enumerate(top_genes):
    #     print(f"{i+1}. {gene}: {similarity_scores[gene]:.4f}")

    return top_genes







def MSTI_Guidancegraph(rna, meta, target_gene=None, target_meta=None, top_n=3):


    print('Creating graph...')
    G = create_basic_graph(rna, meta)

    print('Adding edges...')
    G = add_edges_with_weight_sign(G, mtadata=meta,target_meta_list=target_meta, stadata=rna,target_rna_list=target_gene,top_n=top_n)


    print('Graph completed!')
    return G

def MSTI_kegg_Guidancegraph(rna, meta, kegg_reaction,sample_mt, kegg_rna,target_gene=None, target_meta=None, top_n=3,):



    if target_gene is not None and target_meta is not None:
        kegg_reaction = process_kegg_reaction(kegg_reaction, sample_mt, kegg_rna)
        print('Creating KEGG reaction graph...')
        G = create_graph_from_data(rna, meta, kegg_reaction)
        print('KEGG reaction graph completed!')

        print('Adding edges...')
        G = add_edges_with_weight_sign(G, mtadata=meta,target_meta_list=target_meta, stadata=rna,target_rna_list=target_gene,top_n=top_n)

    else:
        kegg_reaction = process_kegg_reaction(kegg_reaction, sample_mt, kegg_rna)
        print('Creating KEGG reaction graph...')
        G = create_graph_from_data(rna, meta, kegg_reaction)
        print('KEGG reaction graph completed!')


    print('Graph completed!')
    return G