from typing import Union
from pysam import Fastafile
from anndata import AnnData
from mudata import MuData
import numpy as np
from scipy import sparse
from tqdm import tqdm
import re

# TODO
def get_genome(genome:str="hg38"):
    """_summary_

    Args:
        genome (str, optional): _description_. Defaults to "hg38".
    """

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