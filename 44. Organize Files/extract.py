import sys
import os

if len(sys.argv) == 2:
    if os.path.isdir(sys.argv[1]):
        base_folder_path = sys.argv[1]
    else:
        sys.exit("Exiting!!! The folder does not exist. Please enter a valid folder location.")
else:
    sys.exit("Exiting!!! Please enter a folder location as 'python code_location folder_location'.")

# function to extract fasta    
def extract_fasta(base_folder_path):
    # iterate over all files in a folder
    for filename in os.listdir(base_folder_path):
        if filename.strip().endswith(".txt"):
            with open(filename, "r") as f:
                lines = f.readlines()
                # separate index with '>'
                index = [d for d in range(len(lines)) if lines[d].strip().startswith(">")]

                # extract title and sequence
                for i in range(len(index)):
                    if index[i] == index[-1]:
                        title = lines[index[i]][1:].replace("\n", "")
                        sequence = "".join(lines[index[i]+1:]).replace("\n", "")
                    else:
                        title = lines[index[i]][1:].replace("\n", "")
                        sequence = "".join(lines[index[i]+1:index[i+1]]).replace("\n", "")

                    # create directory to save files
                    if base_folder_path.strip().endswith("/"):
                        directory = base_folder_path + filename.strip().split(".txt")[0]
                    else:
                        directory = base_folder_path + "/" + filename.strip().split(".txt")[0]
                        
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                        
                    # write to individual file
                    with open(directory + "/" + title.split()[0], "w") as writer:
                        writer.write(f">{title}\n")
                        writer.write(sequence)    

if __name__ == "__main__":
    
    # calling the function with folder name
    extract_fasta(base_folder_path)        

    print("Successfully extracted files!!!")