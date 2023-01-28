from typing import Union
from anndata import AnnData
from mudata import MuData
import numpy as np
from scipy import sparse
from multiprocessing import Pool, cpu_count
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def compute_deviations(data: Union[AnnData, MuData], n_jobs=-1):
    """
    Compute raw and bias-corrected deviations.

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        n_jobs:
            Number of cpus used for motif matching. If set to -1, all cpus will be used. Default: -1
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

    if isinstance(adata.X, sparse.csr_matrix):
        adata.X = np.array(adata.X.todense())

    logging.info('computing expectation reads per cell and peak...')
    expectation = compute_expectation(count=adata.X)

    logging.info('computing observed motif deviations...')
    motif_match = adata.varm['motif_match'].transpose()

    obs_dev = _compute_dev(
        (motif_match, adata.X.transpose(), expectation.transpose())).transpose()

    # compute background deviations for bias-correction
    logging.info('computing background deviations...')

    n_bg_peaks = adata.varm['bg_peaks'].shape[1]
    bg_dev = np.zeros(shape=(n_bg_peaks, adata.n_obs, len(
        adata.uns['motif_name'])), dtype=np.float32)

    if n_jobs == -1:
        n_jobs = cpu_count()

    if n_jobs == 1:
        for i in range(n_bg_peaks):
            bg_peak_idx = adata.varm['bg_peaks'][:, i]
            bg_motif_match = adata.varm['motif_match'][bg_peak_idx, :].transpose()
            bg_dev[i, :, :] = _compute_dev((bg_motif_match, adata.X.transpose(),
                                            expectation.transpose())).transpose()

    elif n_jobs > 1:
        # prepare arguments for multiprocessing
        arguments_list = list()
        for i in range(n_bg_peaks):
            bg_peak_idx = adata.varm['bg_peaks'][:, i]
            bg_motif_match = adata.varm['motif_match'][bg_peak_idx, :].transpose()
            arguments = (bg_motif_match, adata.X.transpose(), expectation.transpose())
            arguments_list.append(arguments)

        # run the function with multiple cpus
        with Pool(processes=n_jobs) as pool:
            all_results = pool.map(_compute_dev, arguments_list)

        # parse the results
        for i in range(n_bg_peaks):
            bg_dev[i, :, :] = all_results[i].transpose()

    mean_bg_dev = np.mean(bg_dev, axis=0)
    std_bg_dev = np.std(bg_dev, axis=0)

    dev = (obs_dev - mean_bg_dev) / std_bg_dev
    dev = np.nan_to_num(dev, 0)

    dev = AnnData(dev, dtype=np.float32)
    dev.obs_names = adata.obs_names
    dev.var_names = adata.uns['motif_name']

    return dev


def _compute_dev(arguments):
    motif_match, count, expectation = arguments

    observed = np.dot(motif_match, count)
    expected = np.dot(motif_match, expectation)

    return ((observed - expected) / expected)


def compute_expectation(count: np.array) -> np.array:
    """
    Compute expetation accessibility per peak and per cell by assuming identical 
    read probability per peak for each cell with a sequencing depth matched to that cell
    observed sequencing depth.
    Args:
        count (_type_): _description_
    """

    a = np.sum(count, axis=0, keepdims=True)
    a /= np.sum(count)

    b = np.sum(count, axis=1, keepdims=True)

    exp = np.dot(b, a)

    return exp
