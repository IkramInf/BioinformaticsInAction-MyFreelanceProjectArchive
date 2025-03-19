import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser(description='Extracting Fasta Files...')
parser.add_argument('-f', '--folder', help='Input folder name.')
parser.add_argument('-n', '--num', type=int, help='start index')
parser.add_argument('-t', '--title', default="Need title here", help='custom title')
args = parser.parse_args()

base_folder_path, Num, custom_title = args.folder, args.num, args.title

# function to create proper name
def make_name(title):
    title = title.strip().replace(",", "").replace("/", "").replace("|", "").split()
    #return title[0] + " " + "_".join([t for t in title[1:] if t.isalnum()])
    return title[0]

# function to extract fasta    
def extract_fasta(base_folder_path):
    
    if not base_folder_path.strip().endswith("/"):
        base_folder_path = base_folder_path + "/"
    
    # iterate over all files in a folder
    for filename in os.listdir(base_folder_path):
        if filename.strip().endswith(".txt"):
            with open(base_folder_path+filename, "r") as f:
                lines = f.readlines()
                # separate index with '>'
                index = [d for d in range(len(lines)) if lines[d].strip().startswith(">")]

                # extract title and sequence
                for n, i in enumerate(range(len(index))):
                    if index[i] == index[-1]:
                        title = lines[index[i]][1:].replace("\n", "")
                        sequence = "".join(lines[index[i]+1:]).replace("\n", "")
                    else:
                        title = lines[index[i]][1:].replace("\n", "")
                        sequence = "".join(lines[index[i]+1:index[i+1]]).replace("\n", "")

                    # create directory to save files
                    directory = base_folder_path + filename.strip().split(".txt")[0]
                    
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                        
                    # write to individual file
                    temp1 = str(Num+n+1).rjust(5,'0')+" "+title
                    temp = str(Num+n+1).rjust(5,'0')+" "+make_name(title)
                    input_title = temp + " " + custom_title
                    outputfile = directory+"/"+input_title+".txt"
                    with open(outputfile, "w") as writer:
                        writer.write(f">{temp1}\n")
                        writer.write(sequence)    

if __name__ == "__main__":
    
    # calling the function with folder name
    extract_fasta(base_folder_path)        

    print("Successfully extracted files!!!")