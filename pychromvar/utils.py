from typing import Union
from pysam import Fastafile
from anndata import AnnData
from mudata import MuData

from tqdm import tqdm

def add_peak_seq(data: Union[AnnData, MuData], genome_file: str, delimiter="-"):
    """
    Add DNA sequence to data object

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
        adata.uns['seq'][i] = fasta.fetch(chrom, start, end)

    return None

