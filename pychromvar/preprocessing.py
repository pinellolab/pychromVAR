import numpy as np
import scipy as sp
import re

from typing import Union
from anndata import AnnData
from mudata import MuData
from pysam import Fastafile
from tqdm import tqdm
from pynndescent import NNDescent

def get_bg_peaks(data: Union[AnnData, MuData], niterations=50, n_jobs=-1):
    """
    Find background peaks based on GC bias.

    Args:
        data (Union[AnnData, MuData]):
            AnnData object with peak counts or MuData object with 'atac' modality.
        niterations (int, optional): 
            Number of background peaks to sample. Defaults to 50.
        n_jobs:

    Raises:
        TypeError: _description_
    """
    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError(
            "Expected AnnData or MuData object with 'atac' modality")

    # check if the object contains bias in Anndata.varm
    assert "gc_bias" in adata.var.columns, "Cannot find gc bias in the input object, please first run add_gc_bias!"

    reads_per_peak = np.log10(np.sum(adata.X, axis=0))

    # here if reads_per_peak is a numpy matrix, convert it to array
    if isinstance(reads_per_peak, np.matrix):
        reads_per_peak = np.squeeze(np.asarray(reads_per_peak))

    mat = np.array([reads_per_peak, adata.var['gc_bias'].values])
    chol_cov_mat = np.linalg.cholesky(np.cov(mat))
    trans_norm_mat = sp.linalg.solve_triangular(
        a=chol_cov_mat, b=mat, lower=True).transpose()

    index = NNDescent(trans_norm_mat, metric="euclidean", n_neighbors=niterations, n_jobs=n_jobs)
    knn_idx, _ = index.query(trans_norm_mat, niterations)

    adata.varm['bg_peaks'] = knn_idx

    return None


def add_peak_seq(data: Union[AnnData, MuData], genome_file: str, delimiter="-"):
    """
    Add the DNA sequence of each peak to data object. 
    The sequences will be used in GC bias estimation and motif binding sites matching.

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        genome_file (str): 
            Filename of genome reference
        delimiter (str, optional): 
            Delimiter that separates peaks. Defaults to "-".
    """

    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError("Expected AnnData or MuData object with 'atac' modality")

    fasta = Fastafile(genome_file)
    adata.uns['peak_seq'] = [None] * adata.n_vars

    for i in tqdm(range(adata.n_vars)):
        peak = re.split(delimiter, adata.var_names[i])
        chrom, start, end = peak[0], int(peak[1]), int(peak[2])
        adata.uns['peak_seq'][i] = fasta.fetch(chrom, start, end).upper()

    return None


def add_gc_bias(data: Union[AnnData, MuData]):
    """
    Compute GC bias for each peak.

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
    Returns:
        _type_: _description_
    """

    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError("Expected AnnData or MuData object with 'atac' modality")

    assert "peak_seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"

    bias = np.zeros(adata.n_vars)

    for i in tqdm(range(adata.n_vars)):
        seq = adata.uns['peak_seq'][i]

        freq_a = seq.count("A")
        freq_c = seq.count("C")
        freq_g = seq.count("G")
        freq_t = seq.count("T")

        if freq_a + freq_c + freq_g + freq_t == 0:
            bias[i] = 0.5
        else:
            bias[i] = (freq_g + freq_c) / (freq_a + freq_c + freq_g + freq_t)

    adata.var['gc_bias'] = bias

    return None