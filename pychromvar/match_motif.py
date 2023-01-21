from typing import Union
from tqdm import tqdm
import MOODS.scan
import MOODS.tools
import MOODS.parsers
from anndata import AnnData
from mudata import MuData
from pysam import Fastafile
import numpy as np


def match_motif(data: Union[AnnData, MuData], motifs,
            pseudocounts=0.0001, pvalue=0.0001):
    """
    Perform motif matching to predict binding sites

    Args:
        data (Union[AnnData, MuData]): 
            AnnData object with peak counts or MuData object with 'atac' modality.
        motifs: 
            List of motifs
        pseudocounts:
            Pseudocounts
    """

    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError(
            "Expected AnnData or MuData object with 'atac' modality")

    assert "seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"

    # prepare motif data
    matrices = [None] * len(motifs)
    adata.uns['motif_id'] = [None] * len(motifs)
    adata.uns['motif_name'] = [None] * len(motifs)

    for i in range(len(motifs)):
        motif = motifs[i]
        adata.uns['motif_id'][i] = motif.matrix_id
        adata.uns['motif_name'][i] = motif.name
        motif.pseudocounts = {
            'A': pseudocounts, 'C': pseudocounts, 'G': pseudocounts, 'T': pseudocounts}
        pssm = motif.pssm
        matrices[i] = (tuple(pssm['A']),
               tuple(pssm['C']),
               tuple(pssm['G']),
               tuple(pssm['T']))

    matrices_reverse = [MOODS.tools.reverse_complement(matrix) for matrix in matrices]

    # create scanner
    scanner_forward = MOODS.scan.Scanner(7)
    scanner_reverse = MOODS.scan.Scanner(7)

    # compute threshold
    bg = MOODS.tools.flat_bg(4)
    thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]

    scanner_forward.set_motifs(matrices=matrices, bg=bg, thresholds=thresholds)
    scanner_reverse.set_motifs(matrices=matrices_reverse, bg=bg, thresholds=thresholds)

    motif_match_res = np.zeros(shape=(adata.n_vars, len(motifs)))
    for i in tqdm(range(adata.n_vars)):
        results_forward = scanner_forward.scan(adata.uns['seq'][i])
        results_reverse = scanner_reverse.scan(adata.uns['seq'][i])

        # check the matching result for each motif
        for j, res in enumerate(results_forward):
            motif_match_res[i, j] = len(res)

        for j, res in enumerate(results_reverse):
            motif_match_res[i, j] = len(res)

    motif_match_res[motif_match_res > 0] = 1
    adata.varm['motif_match'] = motif_match_res

    return None