import numpy as np
from pychromvar.compute_deviations import *

def test_compute_expectation():
    count = np.array([[1, 0, 1, ], [0, 1, 1]], dtype=np.float32)
    exp = compute_expectation(count)

    # make sure output has same dimensionas
    assert exp.shape == count.shape

    # check the output
    assert np.array_equal(exp, np.array([[0.5, 0.5, 1], [0.5, 0.5, 1]]))
