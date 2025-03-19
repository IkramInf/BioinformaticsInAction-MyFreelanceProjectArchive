# required libraries
import os
import math
import pandas as pd

# Read data
clinical_df = pd.read_csv("new_data_clinical_patient.txt", na_values='[Not Available]',
                         dtype={'Overall Survival (Months)': float, 'Disease Free (Months)': float})
# drop rows with NaN values in column Disease Free Status
clinical_df = clinical_df.dropna(subset=['Disease Free Status'])
clinical_df

# Create two df p1_data & p2_data
p1_data = clinical_df[clinical_df['Disease Free Status'] == "1:Recurred/Progressed"]
p1_data = p1_data.sort_values("Disease Free (Months)")

p2_data = clinical_df[clinical_df['Disease Free Status'] == "0:DiseaseFree"]
p2_data = p2_data.sort_values("Disease Free (Months)")

# create directory to store arrays
#os.makedirs("p1_data", exist_ok=True)
#os.makedirs("p2_data", exist_ok=True)

# Get p1 arrays
p1_arrays = []
p1_n = math.ceil(p1_data['Disease Free (Months)'].max())

for i in range(p1_n + 1):
    p1_mn = p1_data[p1_data['Disease Free (Months)'] <= i]
    p1_arrays.append(p1_mn.values)
    # Save DataFrame to a text file (CSV format)
    with open(f"p1_m{i}.txt", "w") as ifile:
        ifile.write(f"Sample Size: {len(p1_mn.values)}\n")
        p1_mn.to_csv(ifile, index=False)

# Get p2 arrays
p2_arrays = []
p2_n = math.floor(p2_data['Disease Free (Months)'].max())

for i in range(p2_n + 1):
    p2_mn = p2_data[p2_data['Disease Free (Months)'] >= i]
    p2_arrays.append(p2_mn.values)
    # Save DataFrame to a text file (CSV format)
    with open(f"p2_m{i}.txt", "w") as ifile:
        ifile.write(f"Sample Size: {len(p2_mn.values)}\n")
        p2_mn.to_csv(ifile, index=False)

# Display arrays on terminal
num = int(input("Enter a Number: "))
print("Printing arrays of p1_data:")
for i, array in enumerate(p1_arrays):
    Len = len(array)
    if Len >= num:
        print(f"p1_m{i}: {Len}")
        
print("\nPrinting arrays of p2_data:")
for i, array in enumerate(p2_arrays):
    Len = len(array)
    if Len >= num:
        print(f"p2_m{i}: {Len}")
