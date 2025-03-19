# import required libraries
import pandas as pd
from lazypredict.Supervised import LazyClassifier
from sklearn.model_selection import train_test_split
import d1, d2


# load the data
positive = d1.parser("project3_data/pos-Project-Part3-Data.txt")
negative = d1.parser("project3_data/neg-Project-Part3-Data.txt")


# calculate features and make csv file for both positive and negative data
d2.calculate("project3_data/pos-Project-Part3-Data.txt", 2)
d2.calculate("project3_data/neg-Project-Part3-Data.txt", 2)

# read positive csv file
positive = pd.read_csv("pos-Project-Part3-Data.csv", header=None)
# add labels after the last column
positive[positive.shape[1]] = [1]*positive.shape[0]

# read negative csv file
negative = pd.read_csv("neg-Project-Part3-Data.csv", header=None)
# add labels after the last column
negative[negative.shape[1]] = [0]*negative.shape[0]

# merge the positive and negative data
data = pd.concat([positive, negative], axis=0)

# separate train and target data
X = data.iloc[:,:-1]
y = data.iloc[:,-1]

# split the data into train and test portion
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.25, random_state=1)


# train and prediction
clf = LazyClassifier(verbose=0, ignore_warnings=True, custom_metric=None)
models, predictions = clf.fit(X_train, X_test, y_train, y_test)

print(models)


# the 5 top best models
print("\nThe top 5 best models are :\n", models.head())





