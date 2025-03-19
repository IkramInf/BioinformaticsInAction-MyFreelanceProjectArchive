#!/usr/bin/env python
# coding: utf-8
# import required libraries
import os
import gzip
import wget
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline

# sets up local BLAST databases
# runs BLAST searches
# analyses BLAST results
class BLASTAnalysis:
    def __init__(self, urls):
        # take url of three choosen species
        self.urls = urls
        
    def getFasta(self):
        # download and unpack protein fasta file
        unpacked_fnames = []
        for url in self.urls:
            if os.path.exists(url.split("/")[-1]):
                filename = url.split("/")[-1]
            else:
                filename = wget.download(url)
            unpacked_fname = filename.split(".")[0].strip()
            unpacked_fnames.append(unpacked_fname)
            with gzip.open(filename, "rt") as ifile, open(unpacked_fname, "w") as ofile:
                ofile.write(ifile.read())
        return unpacked_fnames
    
    def prepareBlastSearch(self):
        # prepare database
        unpacked_fnames = self.getFasta()
        if not os.path.exists("DB/"):
            os.mkdir("DB/")
        db_names = []
        for fname in unpacked_fnames:
            title = fname.split(".fa")[0]
            db = "DB/" + title
            command = f"makeblastdb -in {fname} -input_type fasta -dbtype prot -title {title} -parse_seqids -out {db}"
            p = subprocess.run(command, shell=True)
            db_names.append(db)
        return (unpacked_fnames, db_names)
            
    def runBlast(self):
        # runs each with each to produce blast results
        fnames = self.prepareBlastSearch()
        names = {}
        for i, f in enumerate(fnames[0]):
            for j, d in enumerate(fnames[1]):
                if i != j:
                    outf = f"{f}{i+1}{j+1}.tbl"
                    cmd = NcbiblastpCommandline(query=f, db=d, evalue=0.001, outfmt=6, out=outf, max_target_seqs=1)
                    qresult = cmd()
                    names[outf] = (f, d)
        return names
    
    def analyseBlast(self):
        # analyses blast results
        names = self.runBlast()
        for tab, qs in names.items():
            print(f"Query: {qs[0]}\nDatabase: {qs[1]}\n{'-'*30}")
            with open(tab, "r") as f:
                results = [line.strip().split("\t") for line in f.readlines()]
            # expect value
            e = min(results, key=lambda x: float(x[-2].split("e")[-1]))
            print(f"Best expect value (e-value): {e[-2]}\n\tQuery seq id: {e[0]}\tSubject seq id: {e[1]}")
            print(f"\tLength: {e[3]}\t\tbit score: {e[-2]}")
            v = int(e[3]) - (int(e[4])+int(e[5]))
            print(f"\tIdentifies: {v}/{e[3]} ({e[2]}%)\tMismatches: {e[4]}/{e[3]}\tGaps: {e[5]}/{e[3]}")
            # bit score
            b = max(results, key=lambda x: float(x[-1]))  
            print(f"Highest alignment score (bit score): {b[-1]}\n\tQuery seq id: {b[0]}\tSubject seq id: {b[1]}")
            print(f"\tLength: {b[3]}\t\te-value: {b[-1]}")
            v = int(b[3]) - (int(b[4])+int(b[5]))
            print(f"\tIdentifies: {v}/{b[3]} ({b[2]}%)\tMismatches: {b[4]}/{b[3]}\tGaps: {b[5]}/{b[3]}")
            print(f"Total hits found: {len(results)}\n")


if __name__ == "__main__":
    urls = ["https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/brugia_malayi/PRJNA10729/brugia_malayi.PRJNA10729.WBPS17.protein.fa.gz",
            "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/brugia_pahangi/PRJEB497/brugia_pahangi.PRJEB497.WBPS17.protein.fa.gz",
            "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/brugia_timori/PRJEB4663/brugia_timori.PRJEB4663.WBPS17.protein.fa.gz"]

    # create object for class BLASTAnalysis
    havingblast = BLASTAnalysis(urls)
    # run it to get output
    havingblast.analyseBlast()
    