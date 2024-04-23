"""Module with the main classes/functions for the base calculations of the 3D
analysis of transcriptoms."""

import numpy as np
import pandas as pd
from scipy import spatial

# import pprint  # for testing only!
from great_3d import timing


@timing
def calculate_distance(df, distance_metric="seuclidean"):
    """Take a genes position dataFrame, return the distance matrix as a pandas data frame.
    - Args:
    - `df`: The genes position dataFrame.
    - `distance_metric`: The numpy distance metric to use (default: "seuclidean").

    - Return:
    - The distance matrix as a pandas data frame.
    """
    geneNames = df.index
    da = df[["X", "Y", "Z"]].to_numpy()
    ndarray = spatial.distance.pdist(da, distance_metric, V=None)  # Look at other available implementations from numpy
    matrix_uni = spatial.distance.squareform(ndarray)
    return pd.DataFrame(matrix_uni, index=geneNames, columns=geneNames)


@timing
def sorting_distances(dist_df):
    """Take a genes distance matrix (Pandas DataFrame), return a dictionary of
    gene names sorted by distance with each gene as a key.
    - Args:
    - `dist_df`: The genes distance matrix (Pandas DataFrame).

    - Return:
    - A dictionary of gene names sorted by distance.
    """
    sorted_indices = np.argsort(dist_df.values, axis=1)
    sorted_distances = np.sort(dist_df.values, axis=1)
    gene_names = dist_df.index
    return {
        gene_names[i]: pd.Series(sorted_distances[i], index=gene_names[sorted_indices[i]])
        for i in range(len(gene_names))
    }
    # The above implementation might be more complicated but it is 2x times faster
    # return {
    #     gene_name: dist_df.loc[gene_name, :].sort_values()
    #     for gene_name in dist_df.index
    # }


@timing
def sum_correlation(dists_sorted, geDF, no_genes, correlation_type):
    """Take the dictionnary of the closest genes, the gene expression file, the
    number of genes we need to compute correlation and the correlation type.

    - Args:
    - `dists_sorted`: The dictionary of the closest genes.
    - `geDF`: The gene expression pandas data frame.
    - `no_genes`: The number of genes to compute correlation.
    - `correlation_type`: The correlation type.

    - Return:
    - A complex dictionary containing the sum of correlations, the sum of
    absolute correlations, and the neighboring genes for each gene.
    """
    # FIXME If a gene has zero expression we give it zero correlation immediately
    geNumpy = geDF.to_numpy()
    corrMnp = np.corrcoef(geNumpy)
    corrMat = pd.DataFrame(corrMnp, index=geDF.index, columns=geDF.index)
    correlation_sums = {}
    # Compute correlation matrix once
    for gene_ref, sorted_genes in dists_sorted.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(sorted_genes[1: no_genes+1].index)
        # Select all proximal correlations from the coorelation matrix.
        correlation = corrMat.loc[gene_ref, selected_genes]
        # correlation = np.nan_to_num(correlation)  # Correction in case of NAN values
        sum_correlationA = sum(abs(correlation))
        sum_correlation = sum(correlation)
        correlation_sums[gene_ref] = [sum_correlation, sum_correlationA, selected_genes]
    return correlation_sums
