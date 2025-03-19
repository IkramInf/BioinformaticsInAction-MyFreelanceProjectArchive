# Required Libraries
import numpy as np
import pandas as pd
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

def x_scale(x, p=7.5):
    '''
    function for scaling x1
    argument:
        x: the input variable
        p: the scaling factor (default: 7.5)
    returns:
        the scaled variable
    '''
    return 1/p * np.log(1 + x * (np.exp(p) - 1))

def y_scale(y):
    '''
    function for scaling y1 and y2
    argument:
        y: the input variable
    returns:
        the scaled variable
    '''
    return np.log(1 + y) if y >= 0 else -np.log(1 - y)

# Read the csv file
df = pd.read_csv("dataset.csv.gz", compression="gzip")

# Scale x1 and y1
df['x1'] = df['x1'].apply(x_scale)
df['y1'] = df['y1'].apply(y_scale)

# Define features and target
features = ['x1', 'x2', 'x3', 'x4']
target = 'y1'
X = df[features]
y = df[target]

# Split data into train, test and validation set
# training : 70%, validation : 15%, testing : 15%
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.30, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_val, y_val, test_size=0.50, random_state=42)

print("Shape of Training Data: ", X_train.shape)
print("Shape of Validation Data: ", X_val.shape)
print("Shape of Testing Data: ", X_test.shape)

# Train the XGBoost model and set hyperparameters
model = XGBRegressor(n_estimators=1000,
                     max_depth=7,
                     eta=0.1,
                     subsample=0.7,
                     colsample_bytree=1.0,
                     eval_metric=mean_absolute_error)

model.fit(X_train, y_train, eval_set=[(X_val, y_val)])

# Predict on validation set
y_pred = model.predict(X_val)

# Calculate RMSE on validation set
val_rmse = np.sqrt(mean_squared_error(y_val, y_pred))
print("Final Validation Root Mean Squared Error (RMSE):", val_rmse)

# Evaluate the model on the test set
y_test_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_test_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_test, y_test_pred)
r_squared = r2_score(y_test, y_test_pred)

print("\nDifferent Evaluation score for Testing Data:")
print("Mean Squared Error (MSE):", round(mse, 3))
print("Root Mean Squared Error (RMSE):", round(rmse, 3))
print("Mean Absolute Error (MAE):", round(mae, 3))
print("R-squared (R2) Score:", round(r_squared, 3))
