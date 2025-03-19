# import numpy
import numpy as np
from math import exp

# calculate euclidean distance between two points
def euclideanDistance(p, q):
    return sum([(x-y)**2 for x, y in zip(p, q)])**0.5

# function for Lloyd k-means clustering
def kMeansClustering(data, k):
    # length of data
    N = len(data)
    # choose first k data points as centers to initialize
    centers = data[:k]
    # calculate closest centers for each data points
    closest_centers = [min([(euclideanDistance(point, c), c) for c in centers])[1] for point in data]
    # calculate error
    error = (1/N)*sum([euclideanDistance(p, cc)**2 for p, cc in zip(data, closest_centers)])
    # let, error distortion as 1
    error_distortion = 1
    
    # if error distortion greater than 0.001, loop will continue
    while error_distortion > 0.001:
        # dictionary to store all clusters
        cluster = {}
        for point in data:
            closest = min([(euclideanDistance(point, c), c) for c in centers])[1]
            cluster.setdefault(tuple(closest),[]).append(point)
        
        # calculate new centers and replace into 'centers' variable
        for index, c in enumerate(centers):
            centers[index] = tuple([sum(dp)/len(dp) for dp in zip(*cluster[tuple(c)])])
        
        # calculate error for new centers
        p_error = error
        new_closest = [min([(euclideanDistance(point, c), c) for c in centers])[1] for point in data]
        error = (1/N)*sum([euclideanDistance(p, cc)**2 for p, cc in zip(data, new_closest)])
        # calculate error distortion
        error_distortion = abs(error - p_error)
    
    return centers

# Calculate Centers to Soft Clusters (E-step)
def eStep(data, beta, centers):
    # apply the formula of responsibility matrix
    distances = [[exp(-beta * euclideanDistance(p, c)) for p in data] for c in centers]
    summation = np.sum(distances, axis=0)
    responsibility_matrix = distances/summation
    return responsibility_matrix

# Calculate Soft Clusters to Centers (M-step)
def mStep(data, respMatrix, m):
    N = len(data)
    # calculate new centers based on data and responsibility matrix
    centers = [[sum([data[k][j]*respMatrix[i][k] for k in range(N)])/sum([respMatrix[i][k] \
                            for k in range(N)]) for j in range(m)] for i in range(len(respMatrix))]   
    return centers

# function for implementation of soft k-means clustering
def softKMeansClustering(data, k, beta, centers):
    # iterate over 100 iterations
    for i in range(100):
        # centers to soft clusters
        responsibility_matrix = eStep(data, beta, centers)
        # soft clusters to centers
        centers = mStep(data, responsibility_matrix, m)    
    return centers

        
if __name__ == "__main__":
    with open("rosalind_ba8d.txt", "r") as f:
        lines = f.read().strip().split()
        k = int(lines[0])
        m = int(lines[1])
        beta = float(lines[2])
        data = np.array(lines[3:]).astype(float)
        data = data.reshape(len(data)//m, m)
        data = list(map(list, data))
    
    # calling the k-means clustering function
    centers = kMeansClustering(data, k)
    # calling the soft k-means clustering function 
    soft_centers = softKMeansClustering(data, k, beta, centers)

    soft_centers = np.round(soft_centers, 3).astype(str)
    for center in soft_centers:
        print(" ".join(center))     

        