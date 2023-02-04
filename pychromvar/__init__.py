__version__ = "0.0.4"
__version_info__ = tuple([int(num) for num in __version__.split('.')])  # noqa: F401

from .preprocessing import get_bg_peaks, add_gc_bias, add_peak_seq
from .match_motif import match_motif
from .compute_deviations import compute_deviations, compute_expectation
from .get_genome import get_genome

