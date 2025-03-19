# import required libraries
import os
import re
import argparse

parser = argparse.ArgumentParser(description='Extract genbank file...')
parser.add_argument('-i', '--input', required=True, help='Genbank filename.')
parser.add_argument('-n', '--num', type=int, default=1, help='Start index of filename')
args = parser.parse_args()
Num = args.num

# read informations into a list
infos = []
with open(args.input, "r") as f:
    reads = f.read().split("//")[0:-1]
    for record in reads:
        m = re.search("DEFINITION(.*?)TITLE", record, re.S)
        #s, e = m.start(), m.end()
        definition = m.group(1).strip().replace("\n", " ").replace("  ","")
        fname = definition.replace("/","_").replace("|","_").replace(",","_").replace(" ","_")

        recs = record.split("\n")
        acc = [line.strip().split("ACCESSION")[1].strip() for line in recs if line.strip().startswith("ACCESSION")]
        seq = "".join([recs[i:] for i, line in enumerate(recs) if line.strip().startswith("ORIGIN")][0][1:])
        seq = re.sub("[^a-zA-Z]", "", seq)
        infos.append((fname, definition, acc, seq))

# write to individual files
if not os.path.exists("OUTPUT"):
    os.mkdir("OUTPUT")

for n, info in enumerate(infos):
    fname, definition, acc, seq = info
    with open("OUTPUT/"+str(Num+n+1).rjust(5,'0')+" "+fname+".txt", "w") as ofile:
        ofile.write(f">{acc[0]} {definition}\n{seq}\n")

