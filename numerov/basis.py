def generate_basis(nvals):
    """ generate a | n, l ⟩ basis set for nvals
        
        Nb. quantum numbers are stored as float to allow for quantum defects
    """
    for n in nvals:
        for l in range(n):
            yield (float(n), float(l))