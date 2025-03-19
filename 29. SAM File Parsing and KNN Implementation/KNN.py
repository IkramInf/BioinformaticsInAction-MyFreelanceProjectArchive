# import libraries
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score

# Implementation of k-nearest neighbors algorithm
class KNearestNeighbors():
    def __init__(self, k):
        self.k = k
    
    def fit(self, X_train, y_train):
        self.X, self.y = np.array(X_train), np.array(y_train)
        
    def getNeighbors(self, test):
        distances = np.linalg.norm(self.X - test, axis=1)
        return distances.argsort()[:self.k]
    
    def predict(self, X_test):
        if isinstance(X_test, pd.DataFrame):
            X_test = X_test.to_numpy()
        else:
            X_test = np.array(X_test)
            
        prediction = []
        for test in X_test:
            neighbors_index = self.getNeighbors(test)
            neighbors = self.y[neighbors_index]
            values, counts = np.unique(neighbors, return_counts=True)
            prediction.append(values[np.argmax(counts)])
        return prediction
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--train', help="Enter train Filename...")
    parser.add_argument('-t','--test', help="Enter test Filename...")
    parser.add_argument('-k','--k', type=int, help="Enter #neighbors...")
    parser.add_argument('-o','--output', help="Enter Output Filename...")
    args = parser.parse_args()
    
    
    # read the train and test file
    train = pd.read_csv(args.train, delimiter = "\t")
    test = pd.read_csv(args.test, delimiter = "\t")

    # separate training and target data from train
    X = train.drop(['tissue'], axis=1)
    y = train['tissue']

    # Split the data as, 75% train and 25% test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=1)

    # apply standard scalar processing on data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    
    # fit and predict
    knn = KNearestNeighbors(k = args.k)
    knn.fit(X_train, y_train)
    y_pred = knn.predict(X_test)
    
    print(f"Sklearn KNN Accuracy: {accuracy_score(y_test, y_pred)}")
    
    # predict on test data and write to a file
    prediction = knn.predict(test)
    
    with open(args.output, "w") as f:
        for i in prediction:
            f.write(f"{i}\n")
    
    
