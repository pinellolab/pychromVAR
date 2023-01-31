__version__ = "0.0.3"
__version_info__ = tuple([int(num) for num in __version__.split('.')])  # noqa: F401

from .preprocessing import *
from .match_motif import *
from .compute_deviations import *
from .get_genome import *

