from typing import Union
from pysam import Fastafile
from anndata import AnnData
from mudata import MuData
import numpy as np
from tqdm import tqdm

def add_motif_to_anndata(data: Union[AnnData, MuData], motifs):
    """
    Add motif information to anndata object.

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        motifs: 
            List of motifs
    """

    data.uns['motif_id'] = [None] * len(motifs)
    data.uns['motif_name'] = [None] * len(motifs)

    for i in range(len(motifs)):
        data.uns['motif_id'][i] = motifs[i].matrix_id
        data.uns['motif_name'][i] = motifs[i].name

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
    adata.uns['seq'] = [None] * adata.n_vars

    for i in tqdm(range(adata.n_vars)):
        peak = adata.var_names[i].split(delimiter)
        chrom, start, end = peak[0], int(peak[1]), int(peak[2])
        adata.uns['seq'][i] = fasta.fetch(chrom, start, end).upper()

    return None

def compute_expectation(count: np.array) -> np.array:
    """
    Compute expetation accessibility per peak and per cell by assuming identical 
    read probability per peak for each cell with a sequencing depth matched to that cell
    observed sequencing depth.

    Args:
        count (_type_): _description_
    """

    a = np.sum(count, axis=0)
    a = np.expand_dims(a, axis=0)
    a /= np.sum(count)

    b = np.sum(count, axis=1)
    b = np.expand_dims(b, axis=1)

    exp = np.dot(b, a)

    return exp