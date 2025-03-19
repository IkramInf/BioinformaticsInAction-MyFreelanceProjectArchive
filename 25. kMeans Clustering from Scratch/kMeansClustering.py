# import numpy
import numpy as np

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
        
if __name__ == "__main__":
    with open("rosalind_ba8c.txt", "r") as f:
        lines = f.read().strip().split()
        k = int(lines[0])
        m = int(lines[1])
        data = np.array(lines[2:]).astype(float)
        data = data.reshape(len(data)//m, m)
        data = list(map(list, data))
    # calling the function
    centers = kMeansClustering(data, k)
    
    # format result as output format
    centers = np.round(centers, 3).astype(str)
    for center in centers:
        print(" ".join(center))