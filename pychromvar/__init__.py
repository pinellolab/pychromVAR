__version__ = "0.0.2"
__version_info__ = tuple([int(num) for num in __version__.split('.')])  # noqa: F401

from .add_gc_bias import add_gc_bias
from .get_bg_peaks import get_bg_peaks
from .match_motif import match_motif
from .compute_deviations import compute_deviations
from .utils import add_peak_seq

