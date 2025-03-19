#!user/bin/env python3
"""Test behavior of identify_sequence_refactored.py"""

from identify_sequence_refactored import identify_sequence # import the function to test (make sure you rename your script to "indentify_sequence_refactored.py")

def test_dna_sequence():
    """Identify a DNA sequence"""
    assert identify_sequence("ATGATGA") == "nucleic acid", "expect ATGATGA identifies as nucleic acid"

def test_dna_lowercase_sequence():
    """Identify a DNA sequence in lowercase"""
    assert identify_sequence("atgatga") == "nucleic acid", "except atgatga identifies as nucleic acid"

def test_rna_sequence():
    """Identify an RNA sequence"""
    assert identify_sequence("AUGAUGA") == "rna", "except identifies as RNA"

def test_aminoacid_sequence():
    """Identify an amino acid sequence"""
    assert identify_sequence("VMPMIN") == "amino acid", "except VMPMIN identifies as amino acids"

def test_nonsequence():
    """TODO: write what you think should happen with a non-sequence"""
    assert identify_sequence("ZZXyDI") == "invalid sequence", "neither nucleic acid or amino acid"


test_dna_sequence()
test_dna_lowercase_sequence()
test_rna_sequence()
test_aminoacid_sequence()
test_nonsequence()