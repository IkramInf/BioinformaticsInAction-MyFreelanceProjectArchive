import os
import sys
import io
import gzip
import glob
import time
import pickle
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import LabelEncoder
import openpyxl
from openpyxl import Workbook
from openpyxl import load_workbook

#FUCTIONS CREATED

def copy_and_paste_vcf_data_1file(vcf_path):
    txt_with_alldata = open('PATH', 'a')
    with gzip.open(vcf_path, "rt") as vcf_f:
        for line in vcf_f:
            if line[0] != '#':
                chromosome = str(line.split("\t")[0])
                position = str(line.split("\t")[1])
                ref = str(line.split("\t")[3])
                alt = str(line.split("\t")[4])
                txt_with_alldata.write(chromosome + "\t" + position + "\t" + ref + "\t" + alt + "\n")
        txt_with_alldata.close()

def extract_specific_chromosome(vcf_path, specific_chromosome):
    txt_with_data_chromosome = open('PATH', 'a')
    with open(vcf_path, "rt") as vcf_f:
        for line in vcf_f:
            if line[0] == str(specific_chromosome):
                txt_with_data_chromosome.write(line + '\n')

def chromosome_data_merger(chromosomevalue, vcf_path, file):
    file_to_open = open(file, 'a')
    with gzip.open(vcf_path, "rt") as vcf_f:
        for line in vcf_f:
            if line[0] != '#':
                chromosome = str(line.split('\t')[0])
                if chromosome == str(chromosomevalue):
                    file_to_open.write(str(chromosome) + '\t' + str(line.split('\t')[1]) + '\n')

def get_sample_data(vcf_path, windows_range):
    ancestor1 = 0
    ancestor2 = 0
    chromosome = ''
    min = 0
    min2 = 0
    max = 0
    list_graph_ancestor1 = []
    list_graph_ancestor2 = []
    list_graph_x_fragments = []
    x_ranges = []
    with gzip.open(vcf_path, "rt") as vcf_f:
        for line in vcf_f:
            if line[0] != '#':
                chromosome = str(line.split("\t")[0])
                break
        for line in vcf_f:
            if line[0] != '#':
                info_field_line = line.split("\t")[9]
                info_field_line_array = info_field_line.split("|")
                for i in info_field_line_array:
                    max = int(line.split("\t")[1])
                    if int(info_field_line_array[0].rstrip()) == 1:
                        ancestor1 += 1
                    elif int(info_field_line_array[1].rstrip()) == 1:
                        ancestor2 += 1
                    else:
                        pass
                if min + windows_range > max:
                    continue
                else:
                    list_graph_ancestor1.append(ancestor1)
                    list_graph_ancestor2.append(ancestor2)
                    ancestor1 = 0
                    ancestor2 = 0
                    min = max
                    list_graph_x_fragments.append(max)

#CODE

#ANALYZE ALL FILES

root_dir = "/mnt/iribhm/people/sagnetti"
patient = "3c86ba21"

chromosome_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

#location_indels = [fn for fn in glob.glob(root_dir + '/SNV/snv_mnv/') if patient in fn]
#location_snvs = [fn for fn in glob.glob(root_dir + '/SNV/indel/') if patient in fn]

#print(location_indels)
#print(location_snvs)
#print(location_indels_graylist)
#print(location_snvs_graylist)

location_indels = root_dir + '/SNV/indel/3c86ba21-7b11-4ec7-9d20-a2325197c676.consensus.20161006.somatic.indel.vcf.gz'
location_snvs = root_dir + '/SNV/snv_mnv/3c86ba21-7b11-4ec7-9d20-a2325197c676.consensus.20160830.somatic.snv_mnv.vcf.gz'

indelstxt = [fn for fn in glob.glob(root_dir + '/SNV/graylist/indel/*')]
snvstxt = [fn for fn in glob.glob(root_dir + '/SNV/graylist/snv_mnv/*')]

#location_sample = root_dir + "/TEST/paca/*/*/*/*/*/*/*/*/*.vcf.gz"
#location_samples = root_dir + "/TEST/paca/*/*/*/*/*/*/*/*/*.vcf.gz"
print(indelstxt)
print(snvstxt)

file_names = glob.glob(os.getcwd()+"/TEST/paca/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
#print("Total number of files in PACA directory are : ", len(file_names))

#Part to locate samples' data

location_samples = []

for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    ID = code.split("-")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit() and flag != str(23):
        #fname = "buca/" + ch + "_" + str(i) + ".csv"
        #fname = "buca/" + ch + "_" + code + "_" + str(i) + ".csv"
        #df = vcf_parser(filename)
        #df.to_csv(fname, index=False)
        if ID == "3c86ba21":
            location_samples.append(filename)

#Part to create a file with all indels merged from 1 specific chromosome (1 for each chromosome)
for chromosome in chromosome_list:
    for root, directories, files in os.walk(location_indels):
        for file in files:
            print(file)
            filetopen = os.path.join(root, file)
            chromosome_data_merger(chromosome, filetopen, indelstxt[0])

#Part to create a file with all snvs merged from 1 specific chromosome (1 for each chromosome)
for chromosome in chromosome_list:
    for root, directories, files in os.walk(location_snvs):
        for file in files:
            print(file)
            filetopen = os.path.join(root, file)
            chromosome_data_merger(chromosome, filetopen, snvstxt[0])

#Part to create a file with all indels of the graylist merged from 1 specific chromosome (1 for each chromosome)
"""for chromosome in chromosome_list:
    for root, directories, files in os.walk(location_indels_graylist):
        for file in files:
            print(file)
            filetopen = os.path.join(root, file)
            chromosome_data_merger(chromosome, filetopen, indelstxt[0])"""

#Part to create a file with all snvs of the graylist merged from 1 specific chromosome (1 for each chromosome)
"""for chromosome in chromosome_list:
    for root, directories, files in os.walk(location_snvs_graylist):
        for file in files:
            print(file)
            filetopen = os.path.join(root, file)
            chromosome_data_merger(chromosome, filetopen, snvstxt[0])"""

#Part to check 1vs1 each position of each chromosome

#Comparision with indels

counter = 0
coincidences = []
coincidences_txt = root_dir + "coincidences.txt"

for file in location_samples:
    for line in gzip.open(file, 'rt'):
        for line2 in gzip.open(location_indels, 'rt'):
            if line.split('\t')[0] == line2.split('\t')[0]:
                print('Coincidence Found: ' + line.split('\t')[0] + ' --- ' + line2.split('\t')[0] + ' in position ' + str(line2.split('\t')[1]))
                counter += 1
                coincidences.append(line.split('\t')[0] + ' : ' + line.split('\t')[0])
    print('You have found:  ' + str(counter) + ' coincidences in position ' + str(line2.split('\t')[1]))

#Comparision with snvs

counter = 0

for file in location_samples:
    for line in gzip.open(file,'rt'):
        for line2 in gzip.open(location_indels, 'rt'):
            if line[0] == '#':
                pass
            elif line.split('\t')[0] != line2.split('\t')[0]:
                pass
            elif line.split('\t')[0] == line2.split('\t')[0] and line.split('\t')[1] == line2.split('\t')[1]:
                counter += 1
                coincidences.append(line.split('\t')[0] + ' : ' + line.split('\t')[1] + ' : ' + line.split('\t')[3] + ' : ' + line.split('\t')[4] + ' : ' + line.split('\t')[9] )
                print('You have found:  ' + str(counter) + ' coincidences in position ' + str(line2.split('\t')[1]))

print('IN TOTAL YOU HAVE FOUND:  ' + counter + '  COINCIDENCES.')
print(coincidences)

with open(coincidences_txt, 'w') as resumefile:
    for coincidence in coincidences:
        resumefile.write(str(coincidence) + '\n')

#Part to do the % of ancestors

with open(coincidences_txt, 'r') as resumefile:
    ancestor1 = 0
    ancestor2 = 0
    chromosome = ""
    value1 = 0
    value2 = 0
    for line in resumefile.readlines():
        chromosome = line.split('\t')[0]
        values = line.split('\t')[4]
        value1 = values.split('|')[0]
        value2 = values.split('|')[1]
        if value1 == 1:
            ancestor1 += 1
        elif value2 == 1:
            ancestor2 +=1
    print('% of Ancestor 1 in Chromosome ' + str(chromosome) + 'is: ' + str(ancestor1 / (ancestor1 + ancestor2) * 100))
    print('% of Ancestor 2 in Chromosome ' + str(chromosome) + 'is: ' + str(ancestor2 / (ancestor1 + ancestor2) * 100))