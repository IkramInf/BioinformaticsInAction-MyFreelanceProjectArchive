
import math
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']  ############# This list is the codon labels to the coding and non-coding models

def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")   #64x64 matrix of probabilities for coding triplet mutations
    noncodingMatrix = getProbs("./noncodingModel.tab") #64x64 matrix of probabilities for non - coding triplet mutations
    # In these matrices the rows are ancestor codon and the colums are spacii codons

    id2ancestorSeq = getSeq("./Ancestor.fa") # reads the ancestor sequences
    id2spaciiSeq = getSeq("./Spacii.fa") # reads the ancestor sequences
    allID = list(id2ancestorSeq.keys())  # the two sequences above share dictionary indexes
    iscoding = {"likely coding" : [], "likely non-coding" : []}
    likelihood = {}

    for ID in allID:
        cScore = 0 #variable to contain the score of id2ancestorSeq[ID] and id2spaciiSeq[ID] with the coding model
        nScore = 0 #variable to contain the score of id2ancestorSeq[ID] and id2spaciiSeq[ID] with the non-coding model

        ####YOUR CODE GOES HERE. Compute the log Odds for both models. Which sequences came from each model?
        ####In other words, which genes are still coding for protein?

        ####Hint 1: Don't multiply probabilities for each and THEN take the log.
        ####Take the log of each probability and THEN sum. Also the logarithm is math.log(x)

        ####Hint 2: You can find the index of an element using the .index() function
        ####For example a =["a","b","c","d"]
        ####a.index("c") returns 2

        ####Hint 3: The modelCodons list represents the row and column labels.
        ####In the coding matrix,
        ####codingMatrix[modelCodons.index("TTT")][modelCodons.index("TTA")]
        ####is the probability of transition from TTT in ancestor to TTA in M. Spacii
        
        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]
        
        for i in range(0, len(ancestorSeq), 3):
            ancestorCodon = ancestorSeq[i:i+3]
            spaciiCodon = spaciiSeq[i:i+3]

            if (ancestorCodon in modelCodons) and (spaciiCodon in modelCodons):
                cScore += math.log(codingMatrix[modelCodons.index(ancestorCodon)][modelCodons.index(spaciiCodon)])
                nScore += math.log(noncodingMatrix[modelCodons.index(ancestorCodon)][modelCodons.index(spaciiCodon)])

        if cScore > nScore:
            #print(f"{ID} is coding {cScore} {nScore}")
            iscoding["likely coding"].append(ID)
            likelihood[ID] = cScore
        else:
            #print(f"{ID} is NOT coding {cScore} {nScore}")
            iscoding["likely non-coding"].append(ID)
            likelihood[ID] = nScore
            
    for label in iscoding:
        sequences = ", ".join(iscoding[label])
        print(f"List of {label} sequence IDs:\n{sequences}\n")
            
    return (iscoding, likelihood)

def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    return pMatrix

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = (line[1:].strip().split("|")[0])
            id2seq[currkey] = ""
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    f.close()
    return id2seq

iscoding, likelihood = scoreModels()  # call the scoreModels function

id2spaciiSeq = getSeq("./Spacii.fa")
true_labels = []
pred_labels = []
for label in id2spaciiSeq:
    if "_c_" in label:
        true_labels.append(1)
    else:
        true_labels.append(0)
    pred_labels.append(likelihood[label])

# Fit and transform the data
Min, Max = min(pred_labels), max(pred_labels)
pred_labels = [(x - Min) / (Max - Min) for x in pred_labels]
print(true_labels, pred_labels)

# Create the ROC curve
fpr, tpr, thresholds = roc_curve(true_labels, pred_labels)
# Calculate the AUC
roc_auc = auc(fpr, tpr)

# Plot the ROC curve
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, label='ROC curve (AUC = %0.2f)' % roc_auc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve of Ancestor Vs Spacii')
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("ROC Curve of Ancestor Vs Spacii.png", dpi=300)
plt.show()
