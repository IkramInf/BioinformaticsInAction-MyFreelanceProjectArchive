import os
import pytest
from assignment5 import io_utils
#from assignment5.io_utils import get_filehandle, is_gene_file_valid

_DIRECTORY_FOR_UNIGENE = "/data/PROGRAMMING/assignment5"
#_DIRECTORY_FOR_UNIGENE = "./data_directories"

def test_get_filehandle_reading():
    # does it open a file for reading
    # test
    test = io_utils.get_filehandle(_DIRECTORY_FOR_UNIGENE+"/Homo_sapiens/TGM1.unigene", "r")
    assert hasattr(test, "read") is True, "Not able to open for reading"


def test_get_filehandle_writing():
    # does it open a file for writing
    # test
    test = io_utils.get_filehandle('test.txt', mode="w")
    assert hasattr(test, "write") is True, "Not able to open for writing"
    test.close()
    os.remove('test.txt')
    
def test_is_gene_file_valid():
    # is it a valid gene file location
    test = io_utils.is_gene_file_valid(_DIRECTORY_FOR_UNIGENE+"/Homo_sapiens/TGM1.unigene")
    assert type(test) == bool, "The path doesn't exist"