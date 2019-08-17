#! python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 21:48:03 2019

@author: Adam
"""

def generate_basis(nvals):
    """ generate a | n, l ) basis set for nvals

        Nb. quantum numbers are stored as float to allow for quantum defects
    """
    for n in nvals:
        for l in range(n):
            yield (float(n), float(l))
