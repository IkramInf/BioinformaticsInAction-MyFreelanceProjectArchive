# import required libraries
import sys

# check for correct number of command line arguments
# and the existence of each file
if len(sys.argv) == 6:
    SIF, GO1, GO2, MAP1, MAP2 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    try:
        f1 = open(SIF, "r")
        f1.close()
        f2 = open(GO1, "r")
        f2.close()
        f3 = open(GO2, "r")
        f3.close()
        f4 = open(MAP1, "r")
        f4.close()
        f5 = open(MAP2, "r")
        f5.close()
    except FileNotFoundError:
        sys.exit("Please check the path of each file and put all of the paths correctly.")        
else:
    sys.exit("NotEnoughArgumentsError!!! Please provide five comand line arguments (<SIF> <GO1> <GO2> <MAP1> <MAP2>).")


def get_mapping(map_file):
    # Open the file.
    f = open(map_file, "r")

    # Result is a list of dictionaries.
    mapping_list = []

    # Skip the header on the first line.
    header = f.readline()

    sp_ens, tr_ens, ID_ens = {}, {}, {}
    for line in f:
        # TODO: PUT YOUR CODE HERE
        items = line.split("\t")
        n_cols = len(items)
        
        if n_cols == 4:
            sp_ens[items[1].strip()] = items[0]
            tr_ens[items[2].strip()] = items[0]
            ID_ens[items[3].strip()] = items[0]
        elif n_cols == 3:
            sp_ens[items[1].strip()] = items[0]
            tr_ens[items[2].strip()] = items[0]
        else:
            ID_ens[items[1].strip()] = items[0]
    
    mapping_list.append(sp_ens)
    mapping_list.append(tr_ens)
    mapping_list.append(ID_ens)
    # Remember to close the file after we're done.
    f.close()

    return mapping_list


def get_go_terms(mapping_list, go_file):
    # Open the file.
    f = open(go_file, "r")

    # This will be the dictionary that this function returns.
    # Entries will have as a key an Ensembl ID and the value will
    # be a set of GO terms.
    go_dict = dict()

    for line in f:
        # TODO: PUT YOUR CODE HERE
        if not line.strip().startswith("!"):
            items = line.strip().split("\t")
            uniprot_id, go_term = items[1], items[4]
            
            if uniprot_id in mapping_list[0].keys():
                go_dict.setdefault(mapping_list[0][uniprot_id], set()).add(go_term)
            if uniprot_id in mapping_list[1].keys():
                go_dict.setdefault(mapping_list[1][uniprot_id], set()).add(go_term)
            if uniprot_id in mapping_list[2].keys():
                go_dict.setdefault(mapping_list[2][uniprot_id], set()).add(go_term)

    # Remember to close the file after we're done.
    f.close()

    return go_dict


def compute_score(alignment_file, go_one_dict, go_two_dict):
    # Open the file.
    f = open(alignment_file, "r")

    # Keep track of the number of proteins we can't map to GO terms
    # and the score.
    unmappable_one = 0
    unmappable_two = 0
    score = 0.0

    for line in f:
        # TODO: PUT YOUR CODE HERE
        ens1, ens2 = line.strip().split()
        
        if ens1 in go_one_dict.keys():
            go1 = go_one_dict[ens1]
        else:
            unmappable_one += 1
            
        if ens2 in go_two_dict.keys():
            go2 = go_two_dict[ens2]
        else:
            unmappable_two += 1
            
        score += len(go1.intersection(go2)) / len(go1.union(go2))

    # Remember to close the file after we're done.
    f.close()

    # Return the statistics and the score back so the main code
    # can print it out.
    return unmappable_one, unmappable_two, score


def main():
    # TODO: PUT YOUR CODE HERE
    mapping_list1 = get_mapping(MAP1)
    mapping_list2 = get_mapping(MAP2)
    go_one_dict = get_go_terms(mapping_list1, GO1)
    go_two_dict = get_go_terms(mapping_list2, GO2)
    print(compute_score(SIF, go_one_dict, go_two_dict))


if __name__ == '__main__':
    main()

