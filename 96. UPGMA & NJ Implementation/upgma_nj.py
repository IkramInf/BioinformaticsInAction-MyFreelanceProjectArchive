#!/usr/bin/env python
# coding: utf-8

# import required libraries
import os
import sys
from collections import defaultdict

# ==============================================================================================================
# Unweighted Pair Group Method with Arithmetic Mean (UPGMA) Algorithm
# ==============================================================================================================

# dict of dict to hold initial distance matrix value
hashOriginal = defaultdict(dict)
# dict of dict to hold matrix value as {sequence1_code : {sequence2_code : distance}, ...}
hashClusters = defaultdict(dict)
# dict table to hold Newick format of cluster
hashNewick = defaultdict()
# dict to yield back original index of codes
codehash = defaultdict()
# hold the distance matrix
distmatrix = []

# input filename
# check for existance of file path
# exit program if file doesn't exist otherwise read the file
fileinput = input("Enter the name of the distance matrix file: ")  # DM-p127.txt
if os.path.exists(fileinput):
    with open(fileinput, "r") as f:
        numOTU = int(f.readline().strip())  # number of OTUs
        codes = f.readline().strip().split("\t")  # contains character codes for OTUs
        distmatrix = list(map(lambda x:x.split("\t"), f.read().strip().splitlines()))
        distmatrix = [list(map(int, mat)) for mat in distmatrix]

        for i in range(numOTU):
            codehash[codes[i]] = i
            hashNewick[codes[i]] = codes[i]
            hashClusters[codes[i]][0] = codes[i]
else:
    sys.exit("Error opening file!")

for i in range(numOTU):
    for j in range(numOTU):
        hashOriginal[codes[i]][codes[j]] = distmatrix[i][j]
    
print("\n")
tmp = '\t'.join(codes)
print(f"\t{tmp}")
for i, code in enumerate(codes):
    tmp = '\t'.join(list(map(str, distmatrix[i])))
    print(f"{code}:\t{tmp}")
print("\n")

print(f"{'='*30}\nResults for UPGMA\n{'='*30}\n")
MAXDIST = 999999
numClusters = numOTU

while numClusters > 1:
    # Calculate distances and find smallest
    smallest, smallestI, smallestJ = MAXDIST, 0, 0
    arrayClusters = list(hashClusters.keys())
    
    for i in range (numClusters-1):
        for j in range(i+1, numClusters):
            # Debugging:
            print(f"Cluster {i}: {arrayClusters[i]} Cluster {j}: {arrayClusters[j]}")
            tempdist = 0;
            # The following code determines distance between
            # clusters as the Unweighted Average of all distances.
            for k in range(len(arrayClusters[i])):
                for m in range(len(arrayClusters[j])):
                    tempdist = tempdist + hashOriginal[hashClusters[arrayClusters[i]][k]][hashClusters[arrayClusters[j]][m]]
                    
            tempdist = tempdist / (len(arrayClusters[i]) * len(arrayClusters[j]));
            if tempdist < smallest:
                smallest = tempdist
                smallestI = i
                smallestJ = j

    # merge smallestI and smallestJ into clusters and add
    # to hash tables
    clusterI = arrayClusters[smallestI]
    clusterJ = arrayClusters[smallestJ]
    merge =  clusterI + clusterJ
    print(f"Merging Clusters: {clusterI} and {clusterJ} with distance {smallest}\n")
    i=0;
    for j in range(len(clusterI)):
        hashClusters[merge][i] = hashClusters[clusterI][j]
        i += 1
        
    for j in range(len(clusterJ)):
        hashClusters[merge][i] = hashClusters[clusterJ][j]
        i += 1

    hashNewick[merge] = f"({hashNewick[clusterI]},{hashNewick[clusterJ]})"
    # eliminate old clusters, decrement numClusters
    del hashClusters[clusterI]
    del hashClusters[clusterJ]
    del hashNewick[clusterI]
    del hashNewick[clusterJ]
    numClusters -= 1
    
#print hash structures
arrayClusters = list(hashClusters.keys())  #should only be one cluster left
print(f"Newick format: {hashNewick[arrayClusters[0]]}\n\n")



# ==============================================================================================================
# Neighbour Joining (NJ) Algorithm
# ==============================================================================================================

# dict of dict to hold initial distance matrix value
hashOriginal = defaultdict(dict)
# dict of dict to hold matrix value as {sequence1_code : {sequence2_code : distance}, ...}
hashClusters = defaultdict(dict)
# dict table to hold Newick format of cluster
hashNewick = defaultdict()
hashNewick1 = defaultdict()
# hold the distance matrix
distmatrix = []


# check for existance of file path
# exit program if file doesn't exist otherwise read the file
if os.path.exists(fileinput):
    with open(fileinput, "r") as f:
        numOTU = int(f.readline().strip())  # number of OTUs
        codes = f.readline().strip().split("\t")  # contains character codes for OTUs
        distmatrix = list(map(lambda x:x.split("\t"), f.read().strip().splitlines()))
        distmatrix = [list(map(int, mat)) for mat in distmatrix]

        for i in range(numOTU):
            hashNewick[codes[i]] = codes[i]
            hashNewick1[codes[i]] = codes[i]
            hashOriginal[codes[i]] = codes[i]
else:
    sys.exit("Error opening file!")

for i in range(numOTU):
    for j in range(numOTU):
        hashClusters[codes[i]][codes[j]] = distmatrix[i][j]


print(f"{'='*30}\nResults for NJ\n{'='*30}\n")

# create transition matrix from dict of dict
def make_td(d):
    for k in dict(sorted(d.items(), key=lambda l:l[0].lower())):
        print("\t".join(list(map(lambda x: str(round(x, 2)), dict(sorted(d[k].items(), key=lambda l:l[0].lower())).values()))))

step = 1
MAXDIST = 999999
numClusters = numOTU
merged = []
#arrayClusters = list(hashOriginal.keys())

while numClusters > 2:
    print(f"Neighbor Joining Step {step}")
    arrayClusters = list(hashClusters.keys())
    # Calculate r values
    arrayRValues = [0.0]*numClusters
    for i in range(numClusters):
        temp = 0.0
        for j in range(numClusters):
            #since value can be stored as hash{x}{y} or hash{y}{x}
            if i != j:
                temp += hashClusters[arrayClusters[i]][arrayClusters[j]]
            
        # with divide by 2 modification due to distances being calculated twice
        arrayRValues[i] = temp / (numClusters - 2)
    
    print("Average Distance Matrix")
    for i in range(numClusters):
        print(f"r({arrayClusters[i]}): {arrayRValues[i]}", end=" ")
   
    # Determine merging clusters using transformed distances.
    smallest, smallestI, smallestJ = MAXDIST, 0, 0
    arr = defaultdict(dict)
    for k in arrayClusters:
        arr[k][k] = 0.0
    for i in range(numClusters-1):
        for j in range(i+1, numClusters):
            # calculate transformed distance
            tempTD = hashClusters[arrayClusters[i]][arrayClusters[j]] - arrayRValues[i] - arrayRValues[j]
            # check if smallest transformed distance
            #print(f"TD({arrayClusters[i]},{arrayClusters[j]}) = {tempTD}", end=",  ")
            arr[arrayClusters[i]][arrayClusters[j]] = tempTD
            arr[arrayClusters[j]][arrayClusters[i]] = tempTD
            if tempTD < smallest:
                smallest = tempTD
                smallestI = i
                smallestJ = j
                
    print("\n\nTransition Matrix")
    make_td(arr)

    # merge smallestI and smallestJ into clusters and add
    # to hash tables
    clusterI = arrayClusters[smallestI]
    clusterJ = arrayClusters[smallestJ]
    merge = clusterI + clusterJ
    # calculate branch lengths
    branch1 = (hashClusters[clusterJ][clusterI] + arrayRValues[smallestI] - arrayRValues[smallestJ]) / 2.0
    branch2 = (hashClusters[clusterJ][clusterI] + arrayRValues[smallestJ] - arrayRValues[smallestI]) / 2.0

    #print(f"\n\nMerging Clusters: {clusterI} and {clusterJ}")
    #print(f"    Distance between {clusterI} and ancestral node = {branch1}")
    #print(f"    Distance between {clusterJ} and ancestral node = {branch2}\n")
    print(f"\nMerged {clusterI} & {clusterJ}, Transition Distance = {arr[clusterI][clusterJ]}")

    for i in range(numClusters):
        if (arrayClusters[i] not in merged) and (arrayClusters[i] != clusterI) and (arrayClusters[i] != clusterJ):
            #print(arrayClusters[p], clusterI, clusterJ)
            d1 = hashClusters[arrayClusters[i]][clusterI]
            d2 = hashClusters[arrayClusters[i]][clusterJ]
            # calculate new distance value
            hashClusters[merge][arrayClusters[i]] = (d1 + d2 - hashClusters[clusterI][clusterJ]) / 2.0
            hashClusters[arrayClusters[i]][merge] = hashClusters[merge][arrayClusters[i]]
            
    hashNewick[merge] = f"({hashNewick[clusterI]},{hashNewick[clusterJ]})"
    hashNewick1[merge] = f"(({hashNewick1[clusterI]}:{branch1}, {hashNewick1[clusterJ]}:{branch2}):{float(hashClusters[clusterI][clusterJ])})"
            
    merged.extend((clusterI, clusterJ))
    # eliminate hash keys $clusterI and $clusterJ, decrement numClusters
    del hashClusters[clusterI]
    del hashClusters[clusterJ]
    del hashNewick[clusterI]
    del hashNewick[clusterJ]
    del hashNewick1[clusterI]
    del hashNewick1[clusterJ]
    
    print("\nUpdated Distance Matrix")
    arr = defaultdict(dict)
    keys = set(hashClusters.keys())
    for k1 in keys:
        for k2 in keys:
            if k1 == k2:
                arr[k1][k2] = 0.0
            else:
                arr[k1][k2] = hashClusters[k1][k2]
                
    make_td(arr)
   
    step += 1
    numClusters -= 1


#print information
#print(f"Distance between remaining clusters: {hashClusters[arrayClusters[0]][arrayClusters[1]]} {hashClusters[arrayClusters[1]][arrayClusters[0]]}")
# should have two clusters remaining
print("\nNewick format:", end= " ")
keys = list(hashNewick.keys())
print(f"({', '.join(list(hashNewick.values()))})")

print("Tree with Branch Lengths:", end=" ")
print(f"{', '.join(list(hashNewick1.values()))}")

print(f"Final Avg Distance: {hashClusters[keys[0]][keys[1]] + hashClusters[keys[1]][keys[0]]}\n")