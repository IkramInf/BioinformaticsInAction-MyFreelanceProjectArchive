import os
import sys
import requests
import io
import gzip
import glob
import time
import asyncio
import nest_asyncio
nest_asyncio.apply()
import argparse
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from pyensembl import EnsemblRelease
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_extraction.text import TfidfVectorizer

start = time.time()

def numeric_form(df):
    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    df = df.replace(GT_conversion)
    df = df.replace(to_replace=r'\d\|\d', value=3, regex=True)
    return df            

df1 = pd.read_csv("exon_pca_df.csv")
print(df1.shape)
print(df1.head())
df2 = pd.read_csv("BUCA_files.csv")
df2 = df2.iloc[:,-1]
print(df2.shape)
print(df1.head())
df2 = numeric_form(df2)
df2 = df2.T
df2 = np.array(df2).reshape(-1, 1)
#cvec2 = df2.index

igsr = pd.read_csv("/mnt/iribhm/people/nifernan/1000_genome/igsr_samples.tsv", delimiter="\t")
igsr_dict = {name:code for name, code in zip(igsr['Sample name'], igsr['Superpopulation code'])}
cvec = [igsr_dict[Id] for Id in df1['ID']]

print("===== Performing PCA =====")
scaler = MinMaxScaler()
X = scaler.fit_transform(df1.drop('ID', axis=1))
D = X.shape[1]
L = int(df2.shape[0]//D)
df3 = np.array(df2)[:int(D*L)].reshape(-1, D)
df = scaler.transform(df3)

# default value works fine for PCA
pca = PCA()
exon_pca = pca.fit_transform(X)
proj = pca.transform(df)

print("===== Plotting =====")
fig1 = px.scatter(exon_pca, x=0, y=1, color=cvec)
fig2 = px.scatter(proj, x=0, y=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], symbol_sequence='x')
fig2.update_traces(marker={'color':'black', 'size': 8})

fig3 = go.Figure(data=fig1.data + fig2.data)
fig3.show()


if not os.path.exists("images"):
    os.mkdir("images")

# html file
plotly.offline.plot(fig3, filename='images/exon_pca_projection.html')    
# png file
fig3.write_image("images/exon_pca_projection.png")
# pdf file
fig3.write_image("images/exon_pca_projection.pdf")
