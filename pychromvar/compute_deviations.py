from typing import Union
from anndata import AnnData
from mudata import MuData
import numpy as np
from tqdm import tqdm
from scipy import sparse
from multiprocessing import Pool, cpu_count

def compute_deviations(data: Union[AnnData, MuData], n_jobs=-1):
    """
    Compute raw and bias-corrected deviations.

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        n_jobs:
            Number of cpus used for motif matching. If set to -1, all cpus will be used. Default: 1
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

    count = adata.X
    if not isinstance(count, sparse.csr_matrix):
        count = sparse.csr_matrix(count)

    print("Computing expectation reads per cell and peak...")
    a = np.sum(count, axis=0)
    a /= np.sum(count)
    b = np.sum(count, axis=1)
    expectation = sparse.csr_matrix(np.dot(b, a))

    print("Computing observed deviations...")
    obs_dev = _compute_dev((count, expectation, adata.varm['motif_match']))[0]

    # compute background deviations for bias-correction
    print("Computing background deviations...")
    n_bg_peaks = adata.varm['bg_peaks'].shape[1]
    bg_dev = np.zeros(shape=(n_bg_peaks, adata.n_obs, len(adata.uns['motif_name'])), dtype=np.float32)

    if n_jobs == -1:
        n_jobs = cpu_count()

    if n_jobs == 1:
        for i in tqdm(range(n_bg_peaks)):
            bg_peak_idx = adata.varm['bg_peaks'][:, i]
            motif_match = adata.varm['motif_match'][bg_peak_idx, :]
            bg_dev[i, :, :] = _compute_dev((count, expectation, motif_match))[0]

    elif n_jobs > 1:
        # prepare arguments for multiprocessing
        arguments_list = list()
        for i in range(n_bg_peaks):
            bg_peak_idx = adata.varm['bg_peaks'][:, i]
            motif_match = adata.varm['motif_match'][bg_peak_idx, :]
            arguments = (count, expectation, motif_match)
            arguments_list.append(arguments)

        # run the function with multiple cpus
        with Pool(processes=n_jobs) as pool:
            all_results = pool.map(_compute_dev, arguments_list)

        # parse the results
        for i in range(n_bg_peaks):
            bg_dev[i, :, :] = all_results[i]

    mean_bg_dev = np.mean(bg_dev, axis=0)
    std_bg_dev = np.std(bg_dev, axis=0)

    adata.obsm['chromvar_dev'] = obs_dev - mean_bg_dev
    adata.obsm['chromvar_z'] = adata.obsm['chromvar_dev'] / std_bg_dev
    
    return None


def _compute_dev(arguments):
    count, expectation, motif_match = arguments

    observed = np.dot(count, motif_match)
    expected= np.dot(expectation, motif_match)

    return (observed - expected) / expected
