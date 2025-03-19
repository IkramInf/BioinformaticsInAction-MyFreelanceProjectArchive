# function to parse fasta sequence
def parser(filename):
    with open(filename, "r") as f:
        # read the file as a list
        file = f.readlines()
        # find out description line indices, started with '>'
        desp_index = [ind for ind in range(len(file)) if file[ind].startswith(">")]
        # separate descriptions
        descriptions = [line.replace("\n", "")[1:] for line in file if line.startswith(">")]

        # separate sequences
        sequences = []
        for i in range(len(desp_index)):
            if desp_index[i] == desp_index[-1]:
                sequences.append("".join(file[desp_index[i]+1:]).replace("\n", ""))
            else:
                sequences.append("".join(file[desp_index[i]+1:desp_index[i+1]]).replace("\n", ""))
    # store descriptions and sequences in a dictionary            
    container = {d : s for d, s in zip(descriptions, sequences)}
    return container

if __name__ == "__main__":
    container = parser("project2_data/Example-DNA1-Project-Part2.txt")
    print(container)