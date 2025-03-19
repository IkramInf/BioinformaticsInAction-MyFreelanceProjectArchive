import os

def get_filehandle(filename, mode):
    """
    Provide utilities to open a file in reading, writing mode with proper error handling.
    """
    try:
        f = open(filename, mode)
        return f
    except (OSError, ValueError) as error:
        print(f"IOError: Could not open the file: {filename} for type '{mode}'", filename=sys.stderr)
        raise
        
def is_gene_file_valid(filename):
    """
    This function will check to make sure the given file name exists, if it does it return True else it will return False.
    
    Parameters:
        filename : Name of file to check for existence
    Returns:
        Returns True if exists otherwise False
    """
    
    isExist = os.path.exists(filename)
    
    if isExist:
        return True
    else:
        return False