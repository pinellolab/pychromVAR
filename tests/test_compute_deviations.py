import numpy as np
from pychromvar.compute_deviations import compute_expectation, _compute_deviations

def test_compute_expectation():
    count = np.array([[1, 0, 1, ], [0, 1, 1]], dtype=np.float32)
    exp = compute_expectation(count)

    # make sure output has same dimensionas
    assert exp.shape == count.shape

    # check the output
    assert np.array_equal(exp, np.array([[0.5, 0.5, 1], [0.5, 0.5, 1]]))

def test_compute_deviations():
    count = np.array([[1, 0, 1, ], [0, 1, 1]], dtype=np.float32)
    exp = compute_expectation(count)

    motif_match = np.array([[1, 1], [0, 1], [1, 0]], dtype=np.int8)

    dev = _compute_deviations((motif_match, count, exp)).transpose()

    # make sure output has same dimensionas
    assert dev.shape == (count.shape[1], motif_match.shape[0])

    # check the output
    assert np.array_equal(dev, np.array([[0.0, -1, 1], [0, 1, -1], [0, 0, 0]]))