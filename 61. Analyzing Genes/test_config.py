import pytest
from assignment5 import config

_DIRECTORY_FOR_UNIGENE = "/data/PROGRAMMING/assignment5"
#_DIRECTORY_FOR_UNIGENE = "./data_directories"
_FILE_ENDING_FOR_UNIGENE = "unigene"

def test_get_directory_for_unigene():
    test = config.get_directory_for_unigene()
    assert test == _DIRECTORY_FOR_UNIGENE, "Incorrect directory"
    
    
def test_get_extension_for_unigene():
    test = config.get_extension_for_unigene()
    assert test == _FILE_ENDING_FOR_UNIGENE, "Incorrect extension"

    
def test_get_keywords_for_hosts():
    test = config.get_keywords_for_hosts()
    assert hasattr(test, 'keys'), "Not a valid mapping dictionary of hosts"
    
def test_get_error_string_4_EXCEPTION_TYPE():
    try:
        test = open("check_exception.txt", "r")
        test.close()
    except Exception:
        test = config.get_error_string_4_EXCEPTION_TYPE()
    
    assert test, "Don't a error string"