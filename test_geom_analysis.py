"""
Tests for geom_analysis.py
"""

import pytest
import geom_analysis as ga

def test_calculate_distance():
    coord1 = [0, 0, 2]
    coord2 = [0, 0, 0]
    
    observed = ga.calculate_distance(coord1,coord2)

    assert observed == 2.0

def test_bond_check():
    bond1 = 1.4
    bond2 = 1.7

    observed1 = ga.bond_check(bond1) 
    observed2 = ga.bond_check(bond2)

    assert observed1 == True
    assert observed2 == False

def test_bond_check_error():
    bond_length = "a"

    with pytest.raises(TypeError):
        observed = ga.bond_check(bond_length)
 
