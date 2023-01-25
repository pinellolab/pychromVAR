import os
import numpy as np
import scipy as sp
from joblib import Parallel, delayed
from multiprocessing import Pool, cpu_count

from typing import Union, Literal, get_args
from tqdm import tqdm
import MOODS.scan
import MOODS.tools
from anndata import AnnData
from mudata import MuData

from .utils import add_motif_to_anndata

_BACKGROUND = Literal["subject", "genome", "even"]


def match_motif(data: Union[AnnData, MuData], motifs, pseudocounts=0.0001, p_value=5e-05,
                background: _BACKGROUND = "even", genome_file: str = None, n_jobs=1):
    """
    Perform motif matching to predict binding sites using MOODS. 
    This function wraps 

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        motifs: 
            List of motifs
        pseudocounts:
            Pseudocounts for each nucleotide. Default value is 0.0001
        p_value:
            P-value threshold for motif matching. Default: 5e-05
        background:
            Background distribution of nucleotides for computing thresholds from p-value. 
            Three options are available: "subject" to use the subject sequences, "genome" to use the
            whole genome (need to provide a genome file), or even using 0.25 for each base.
            Default: "subject".
        genome_file:
            If background is set to genome, a genome file must be provided. Default: None
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

    assert "seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"

    options = get_args(_BACKGROUND)
    assert background in options, f"'{background}' is not in {options}"

    if background == "genome":
        assert os.path.exists(genome_file), f"{genome_file} does not exist!"

    add_motif_to_anndata(data=adata, motifs=motifs)

    # compute background distribution
    seq = ""
    if background == "subject":
        for i in range(adata.n_vars):
            seq += adata.uns['seq'][i]
        bg = MOODS.tools.bg_from_sequence_dna(seq, 0)
    elif background == "genome":
        # TODO
        bg = MOODS.tools.flat_bg(4)
    else:
        bg = MOODS.tools.flat_bg(4)

    # prepare motif data
    n_motifs = len(motifs)

    matrices = [None] * 2 * n_motifs
    thresholds = [None] * 2 * n_motifs
    for i, motif in enumerate(motifs):
        counts = (tuple(motif.counts['A']),
                  tuple(motif.counts['C']),
                  tuple(motif.counts['G']),
                  tuple(motif.counts['T']))

        matrices[i] = MOODS.tools.log_odds(counts, bg, pseudocounts)
        matrices[i+n_motifs] = MOODS.tools.reverse_complement(matrices[i])

        thresholds[i] = MOODS.tools.threshold_from_p(matrices[i], bg, p_value)
        thresholds[i+n_motifs] = thresholds[i]

    # create scanner
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices=matrices, bg=bg, thresholds=thresholds)
    motif_match_res = np.zeros(shape=(adata.n_vars, n_motifs), dtype=np.uint8)

    if n_jobs == -1:
        n_jobs = cpu_count()

    if n_jobs == 1:
        for i in tqdm(range(adata.n_vars)):
            results = scanner.scan(adata.uns['seq'][i])
            for j in range(n_motifs):
                if len(results[j]) > 0 or len(results[j+n_motifs]) > 0:
                    motif_match_res[i, j] = 1

    elif n_jobs > 1:
        # prepare arguments for multiprocessing
        
        # seqs = [seq for seq in adata.uns['seq']]
        # arguments_list = list()
        #for i in range(adata.n_vars):
            #arguments = (scanner, adata.uns['seq'][i])
        #    arguments = (adata.uns['seq'][i])
        #    arguments_list.append(arguments)
        all_results = Parallel(n_jobs=n_jobs)(delayed(scanner.scan)(seq) for seq in adata.uns['seq'])

        #with Pool(processes=n_jobs) as pool:
        #    all_results = pool.map(scanner.scan, seqs)

        # parse outputs
        for i, results in enumerate(all_results):
            for j in range(n_motifs):
                if len(results[j]) > 0 or len(results[j+n_motifs]) > 0:
                    motif_match_res[i, j] = 1

    adata.varm['motif_match'] = sp.sparse.csr_matrix(motif_match_res)

    return None

def _match_motif(scanner, seq):
    results = scanner.scan(seq)

    return results

