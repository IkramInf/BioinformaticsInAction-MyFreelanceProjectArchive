#!user/bin/env python3
"""Test behavior of interleaved.py"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.getcwd()))
from interleaved import interleave
from Bio import SeqIO

def test_interleaved_list():
    """Interleave two lists"""
    list1 = ["A", "B", "C"]
    list2 = ["1", "2", "3"]
    expected = ["A", "1", "B", "2", "C", "3"]
    assert interleave(list1, list2) == expected, "expect two lists to be interleaved"


def test_interleaved_SeqRecords():
    """Interleave two iterators of SeqRecords.

    Because SeqRecord comparisons are not supported, this test gets LONG.
    """
    file1 = SeqIO.parse("first3reads_Aip02.R1.fastq", "fastq")
    file2 = SeqIO.parse("first3reads_Aip02.R2.fastq", "fastq")

    expected = []
    for record in SeqIO.parse("first3reads_Aip02.interleave_manual.fastq", "fastq"):
        expected.append(record)

    result = interleave(file1, file2)

    # lists are the same size
    assert len(result) == len(expected), "expect the two lists to have the same number of elements"
    assert result[1].id == expected[1].id, "expect the same indexed sequence to have the same ID"
    assert result[2].id == expected[2].id, "expect the next indexed sequence to also be the same"
    
test_interleaved_list()
test_interleaved_SeqRecords()