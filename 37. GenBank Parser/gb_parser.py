# import required libraries
import sys
import os
import re


# assign commandline argument into filename
try:
    filename = sys.argv[1]
except Exception:
    sys.exit(f"Run the script followed by a commandline argument as 'python {sys.argv[0]} genbank_filename'")


#f1, f2, f3 = "CFTR_DNA.gb", "CFTR_mRNA.gb", "CFTR_protein.gp"

# create a custom exceptions
class InvalidGenbankFile(Exception):
    pass

class CaseError(Exception):
    pass


# class to read genbank file
class GenbankParser:
    
    def __init__(self, filename):
        self.filename = filename
        
    def isValid(self, seq):
        valid_characters = {'a', 'g', 'c', 't'}
        seq = set(seq.lower())
        diff = seq.difference(valid_characters)
        if diff:
            return False
        else:
            return True
    
    def read(self):
        try:
            with open(self.filename, "r") as f:
                lines = f.readlines()
                origin_index = [i for i, line in enumerate(lines) if line.strip().startswith("ORIGIN")]
                seq = "".join(lines[origin_index[0]+1:])
                seq = "".join(re.findall("[a-zA-Z].*?", seq))
                
                try:
                    if len(seq) != 0 and self.isValid(seq):

                        feature_index = [i for i, line in enumerate(lines) if line.strip().startswith("FEATURES")]
                        features = lines[feature_index[0]+1:origin_index[0]]

                        def_index = [i for i, line in enumerate(lines) if line.strip().startswith("DEFINITION")]
                        acc_index = [i for i, line in enumerate(lines) if line.strip().startswith("ACCESSION")]
                        definition = " ".join(lines[def_index[0]].strip().split()[1:]) + " " + \
                                            "".join(lines[def_index[0]+1:acc_index[0]]).strip()

                        return definition, features, seq, origin_index

                    else:
                        raise InvalidGenbankFile
                except InvalidGenbankFile:
                    sys.exit(f"The file named '{self.filename}' is not a valid genbank file or doesn't contain dna sequence.")
                    
        except FileNotFoundError:
            sys.exit(f"The file named '{self.filename}' does not exist in the directory.")


# class to extract feature from genbank file
class Feature:
    
    def __init__(self, obj):
        self.definition, self.features, self.seq, self.origin_index = obj.read()
            
    def complement(self, seq):
        seq = seq[::-1]
        complement_pair = {'a':'t', 'c':'g', 'g':'c', 't':'a',
                           'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        complement = "".join(complement_pair.get(base) for base in seq)
        return complement
    
    def extract(self, case="lower"):
        try:
            if case.lower() == "lower":
                self.seq = self.seq.lower()
            elif case.lower() == "upper":
                self.seq = self.seq.upper()
            else:
                raise CaseError
        except CaseError:
            sys.exit("Only 'lower' and 'upper' are valid case")
            
        output_filename = "".join(filename.split('.')[0:-1]) + "_features.txt"
        output = open(output_filename, "w")
        featureLists = ['source', 'mRNA', 'CDS', 'misc_feature',  'exon', 'ncRNA', 'gene']
        
        feature_first_lines = [i for i, line in enumerate(self.features) \
                               if line.strip().split()[0] in featureLists]
        index = feature_first_lines + self.origin_index
        
        feature_second_lines=[[index[idx]+i for i,line in enumerate(self.features[index[idx]:index[idx+1]]) \
                 if line.strip().startswith("/")][0] for idx in range(len(index)-1)]

        loc_string = ["".join(self.features[p:q]) for p, q in zip(feature_first_lines, feature_second_lines)]

        locations = [re.findall("(\d+)[.>]*?(\d+)", loc) for loc in loc_string]

        is_complement = ["complement" in self.features[k] for k in feature_first_lines]
        
        output.write(f"{self.definition}\n\n")
    
        for x, y, z, c in zip(feature_first_lines, feature_second_lines, locations, is_complement):
            output.write(f">{self.features[x].strip().split()[0]} ")
            output.write(f"{self.features[y].strip()}\n")
            sequence = "".join([self.seq[int(l[0])-1:int(l[1])] for l in z])
            if c == True:
                output.write(f"{self.complement(sequence)}\n\n")
            else:
                output.write(f"{sequence}\n\n")
        output.close()


if __name__ == "__main__":
    reader = GenbankParser(filename)
    feat = Feature(reader)
    
    # you may pass an argument for case from ['lower', 'upper']; default is case='lower'
    feat.extract(case = "lower")
