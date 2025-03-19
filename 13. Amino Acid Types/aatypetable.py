import os

def AAtypetable(filelist, outputfile):
    # store non-existant, invalid and errored filenames into invalidFilename list
    invalidFilename = []
    
    # write all amino acids together in a set
    AA = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'}
    
    # run a for loop to iterate over file lists
    for filename in filelist:
        # return true, if file exists
        if os.path.isfile(filename):
            # return 0, if file is empty
            if os.stat(filename).st_size == 0:
                invalidFilename.append(filename)
            else:
                # if amino acids contain any invalid characters
                if set(readAAsequence(filename)).difference(AA):
                    invalidFilename.append(filename)     
        else:
            # if file is not exist
            invalidFilename.append(filename)
            
    # separate valid files from file lists
    validFile = [file for file in filelist if file not in invalidFilename]
    
    # write in a file
    with open(outputfile, "w") as writer:
        writer.write(f"# Filename    Polar    Small    Hydro\n")
        # to iterate over valid files
        for file in validFile:
            # compute polar, small, hydro fractions by the previous AAtypes function
            polarSmallHydro = AAtypes(readAAsequence(file))
            writer.write(f"{file}    {polarSmallHydro[0]}    {polarSmallHydro[1]}    {polarSmallHydro[2]}\n")
                
    return invalidFilename