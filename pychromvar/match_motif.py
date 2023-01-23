import os
import numpy as np

from typing import Union, Literal, get_args
from tqdm import tqdm
import MOODS.scan
import MOODS.tools
import MOODS.parsers
from anndata import AnnData
from mudata import MuData

from .utils import add_motif_to_anndata

_BACKGROUND = Literal["subject", "genome", "even"]


def match_motif(data: Union[AnnData, MuData], motifs, pseudocounts=0.0001, p_value=5e-05,
                background: _BACKGROUND = "subject", genome_file: str = None):
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
        bg = MOODS.tools.bg_from_sequence_dna(seq, pseudocounts)
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
    motif_match_res = np.zeros(shape=(adata.n_vars, n_motifs))

    # for each peak
    for i in tqdm(range(adata.n_vars)):
        results = scanner.scan(adata.uns['seq'][i])
        # for each motif
        for j in range(n_motifs):
            if len(results[j]) > 0:
                motif_match_res[i, j] = 1
            elif len(results[j+n_motifs]) > 0:
                motif_match_res[i, j] = 1

    adata.varm['motif_match'] = motif_match_res

    return None
