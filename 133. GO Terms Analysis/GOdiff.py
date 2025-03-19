"""
    Combine annotation and diff exp file and print common genes along with description and p-value
        a) read annotation file
        b) read  differential expression data file
        c) write output to a file
        d) print one gene of same group of genes with lowest p-value
"""

def read_annotation(filename):
    """
    Read columns 4, 5, 6 from annotation file
    """
    # store annotation data into dict; gene name as key and description as value
    annotation_dict = {}
    # store feature data into dict; feature as key and list of genes as value
    feature_dict = {}
    # open the annotation file as read mode
    with open(filename, "r") as ifile:
        # iterate over each line of the file
        annotations = ifile.readlines()
        for annotation in annotations[1:]:  # skip header
            # split tab-separated line
            annotation = annotation.strip().split("\t")
            # separating columns 4, 5 and 6
            gene_name, description, feature = annotation[3:]
            # store data into dict
            annotation_dict[gene_name] = description
            feature_dict.setdefault(feature, []).append(gene_name)

    # return the results
    return (annotation_dict, feature_dict)

def read_diff_exp(filename):
    """
    Read columns 1 and 5 from diff exp data file
    """
    # store diff exp data into dict; gene name as key and p-value as value
    diffexp_dict = {}
    # open the diff exp file as read mode
    with open(filename, "r") as ifile:
        # iterate over each line of the file
        diff_expressions = ifile.readlines()
        for diff_expression in diff_expressions:
            # split tab-separated line
            diff_expression = diff_expression.strip().split("\t")
            # separating columns 1 and 5
            gene_name, p_value = diff_expression[0], diff_expression[4]
            # store data into dict
            diffexp_dict[gene_name] = p_value

    # return the results
    return diffexp_dict

def write_output_to_file(output_filename, annotation_dict, diffexp_dict, feature_dict):
    """
    Write output to a file
    """
    # open a output file in write mode
    with open(output_filename, "w") as ofile:
        # write header of the file
        ofile.write("gene name\tdescription\tp-value\n")
        # iterate over key, value in annotation_dict
        for gene_name, description in annotation_dict.items():
            # check that gene name is in diffexp_dict
            if gene_name in diffexp_dict:
                # write data into file
                ofile.write(f"{gene_name}\t{description}\t{diffexp_dict[gene_name]}\n")

    # print the one gene of same group of genes on terminal
    print("Same Group of Genes:")
    # iterate over each feature types
    for genes in feature_dict.values():
        # find out lowest p-value
        temp = {gene: float(diffexp_dict[gene]) for gene in genes if gene in diffexp_dict}
        if temp:
            temp = sorted(temp.items(), key=lambda x: x[1])[0]
            print(f"{temp[0]}\t{temp[1]}")

def main():
    """
    Execute all three functions inside the main function
    """
    # read annotation file and store results data into dict
    annotation_dict, feature_dict = read_annotation("GO0003723.genelist")
    # read diff exp file and store results into dict
    diffexp_dict = read_diff_exp("diffexp.tsv")
    # write common genes into a file
    write_output_to_file("GOdiff_studentid.out", annotation_dict, diffexp_dict, feature_dict)

### Execute the Main Function
main()
