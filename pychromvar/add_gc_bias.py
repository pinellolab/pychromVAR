from typing import Union
from anndata import AnnData
from mudata import MuData

from tqdm import tqdm
import numpy as np

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

    assert "seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"

    bias = np.zeros(adata.n_vars)

    for i in tqdm(range(adata.n_vars)):
        seq = adata.uns['seq'][i]

        freq_a = seq.count("A")
        freq_c = seq.count("C")
        freq_g = seq.count("G")
        freq_t = seq.count("T")

        if freq_a + freq_c + freq_g + freq_t == 0:
            bias[i] = 0.5
        else:
            bias[i] = (freq_g + freq_c) / (freq_a + freq_c + freq_g + freq_t)

    adata.varm['gc_bias'] = bias

    return None