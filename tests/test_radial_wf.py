import numpy as np
from numerov.cy.core import radial_wf as radial_wf_cy
from numerov.core import radial_wf as radial_wf_py
from sympy.physics.hydrogen import R_nl

STEP = 0.001
MAX_DIFF = 1e-7

def test_cy(n=10, l=5):
    """ test that python and cython wf are equivilent
    """
    r_py, y_py = radial_wf_py(n, l, step=STEP)
    r_cy, y_cy = radial_wf_cy(n, l, step=STEP)
    assert (r_py == r_cy).all()
    assert (y_py == y_cy).all()


def test_sympy(n=10, l=5):
    """ test that sympy wf and numerov are approximately equal
    """
    r_cy, y_cy = radial_wf_cy(n, l, step=STEP)
    y_sympy = np.array([R_nl(n, l, r).evalf() for r in r_cy])
    max_diff = np.max(np.abs(y_cy - y_sympy))
    assert max_diff < MAX_DIFF