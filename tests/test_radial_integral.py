from numerov import radial_integral
from sympy.physics.hydrogen import R_nl
from sympy import integrate, oo, var

STEP = 0.0001
FRAC_DIFF = 1e-10

def test_radial_integral(n1=12, l1=5, n2=15, l2=4):
    """ compare radial integral with sympy
    """
    var("r")
    integral_sympy = integrate(R_nl(n1, l1, r) * r**3 * R_nl(n2, l2, r), (r, 0, oo)).evalf()
    integral_numerov = radial_integral(n1, l1, n2, l2, step=STEP)
    frac_diff = abs(integral_numerov - integral_sympy) / integral_sympy
    assert frac_diff < FRAC_DIFF
