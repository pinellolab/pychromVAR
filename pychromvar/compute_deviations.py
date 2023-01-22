from typing import Union
from anndata import AnnData
from mudata import MuData
import numpy as np

from .utils import compute_expectation

def compute_deviations(data: Union[AnnData, MuData]):
    """
    Compute raw and bias-corrected deviations

    Args:
        data (Union[AnnData, MuData]): _description_
    """

    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError(
            "Expected AnnData or MuData object with 'atac' modality")

    # check if the object contains bias in Anndata.varm
    assert "bg_peaks" in adata.varm_keys(
    ), "Cannot find background peaks in the input object, please first run get_bg_peaks!"

    expectation = compute_expectation(count=adata.X)
    motif_match = adata.varm['motif_match'].transpose()

    observed = np.dot(motif_match, adata.X.transpose())
    expected= np.dot(motif_match, expectation.transpose())

    obs_dev = ((observed - expected) / expected).transpose()

    # compute background deviations for bias-correction
    n_bg_peaks = adata.varm['bg_peaks'].shape[1]
    bg_dev = np.zeros(shape=(n_bg_peaks, adata.n_obs, len(adata.uns['motif_name'])), dtype=np.float32)

    for i in range(n_bg_peaks):
        bg_peak_idx = adata.varm['bg_peaks'][:, i]
        bg_motif_match = adata.varm['motif_match'][bg_peak_idx, :].transpose()

        bg_observed = np.dot(bg_motif_match, adata.X.transpose())
        bg_expected= np.dot(bg_motif_match, expectation.transpose())

        bg_dev[i, :, :] = ((bg_observed - bg_expected) / bg_expected).transpose()

    mean_bg_dev = np.mean(bg_dev, axis=0)
    std_bg_dev = np.std(bg_dev, axis=0)

    adata.obsm['chromvar'] = (obs_dev - mean_bg_dev) / std_bg_dev
    
    return None