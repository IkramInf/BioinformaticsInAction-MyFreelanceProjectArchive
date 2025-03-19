def openReadingFrame(string):
    position = []
    for i in range(len(string)) :
        # Check for a start amino acid
        if string[i] == "M":
            # Look at amino acid up to a *
            for j in range(i, len(string)):
                # check for a stop amino acid
                if string[j] == "*":
                    position.extend([i, j])
    
    # if position has start and end index
    if len(position) != 0:
        # take the start to end portion from string and return
        return string[position[0]:position[1]]
    else:
        # if position has no index, return empty string
        return ""