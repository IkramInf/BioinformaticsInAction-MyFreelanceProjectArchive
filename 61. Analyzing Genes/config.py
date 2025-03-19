import sys

_DIRECTORY_FOR_UNIGENE = "/data/PROGRAMMING/assignment5"
#_DIRECTORY_FOR_UNIGENE = "./data_directories"
_FILE_ENDING_FOR_UNIGENE = "unigene"

def get_directory_for_unigene():
    """
    Return directory for unigene file
    """
    return _DIRECTORY_FOR_UNIGENE

def get_extension_for_unigene():
    """
    Return extension for unigene file
    """    
    return _FILE_ENDING_FOR_UNIGENE

def get_keywords_for_hosts():
    """
    Return a mapped dictionary of common names with scientific names.
    """
    
    homo_sapiens = "Homo_sapiens"
    bos_taurus = "Bos_taurus"
    equus_caballus = "Equus_caballus"
    mus_musculus = "Mus_musculus"
    ovis_aries = "Ovis_aries"
    rattus_norvegicus = "Rattus_norvegicus"

    host_keywords = {"Homo sapiens" : homo_sapiens,
                    "Human" : homo_sapiens,
                    "Humans" : homo_sapiens,
                    "Bos taurus" : bos_taurus,
                    "Cow" : bos_taurus,
                    "Cows" : bos_taurus,
                    "Equus caballus" : equus_caballus,
                    "Horse" : equus_caballus,
                    "Horses" : equus_caballus,
                    "Mus musculus" : mus_musculus,
                    "Mice" : mus_musculus,
                    "Mouse" : mus_musculus,
                    "Ovis aries" : ovis_aries,
                    "Sheep" : ovis_aries,
                    "Sheeps" : ovis_aries,
                    "Rattus norvegicus" : rattus_norvegicus,
                    "Rat" : rattus_norvegicus,
                    "Rats" : rattus_norvegicus}
    return host_keywords

def get_error_string_4_EXCEPTION_TYPE():
    """
    Return error string
    """
    
    return sys.exc_info()[1]