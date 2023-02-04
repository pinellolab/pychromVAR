import numpy as np
from anndata import AnnData

from pychromvar.preprocessing import add_gc_bias

def test_gc_bias():
    data = AnnData(np.array([[1, 2,], [3, 4]]), dtype=np.float32)
    data.uns["peak_seq"] = ["AAAAAAAA", "CCCCGGGG"]
    add_gc_bias(data=data)

    assert np.array_equal(data.var['gc_bias'].values, np.array([0.0, 1.0]))