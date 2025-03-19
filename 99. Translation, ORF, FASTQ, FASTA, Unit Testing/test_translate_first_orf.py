#!user/bin/env python3
"""Test behavior of translate_first_orf.py"""

from translate_first_orf import find_first_orf 
from translate_first_orf import translate_first_orf 
from Bio.Seq import Seq


def test_short_orf():
    """Identify short ORF"""
    assert find_first_orf("AUGCCCUAG") == "AUGCCCUAG", "expect three codon ORF"

def test_orf_in_orf():
    """Identify first ORF when two present"""
    assert find_first_orf("AUGCUGUAACUGUAG") == "AUGCUGUAA", "expect first complete ORF"

def test_missing_stop_codon():
    """Identify no ORF when missing stop codon"""
    assert find_first_orf("AUGCUG") == "", "expect no ORF in AUGCUG - lacks stop codon"

def test_out_of_frame_stop():
    """Identify no ORF when stop codon is out of frame"""
    assert find_first_orf("AUGAUAA") == "", "expect no ORF in AUGAUAA - stop codon out of frame"

def test_dna_sequence():
    """Identify protein sequence within DNA"""
    assert translate_first_orf(Seq("AAATGCCCTAG")) == "MP*", "expect MP protein within sequence"
    
test_short_orf()
test_orf_in_orf()
test_missing_stop_codon()
test_out_of_frame_stop()
test_dna_sequence()
