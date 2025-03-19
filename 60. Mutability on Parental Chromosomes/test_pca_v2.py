import os
import sys
import requests
import io
import gzip
import glob
import time
import ast
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

df1 = pd.read_csv("exon_pca_df.csv")
df1 = df1.iloc[:, 1:]
print(df1.shape)
print(df1.head())

df2 = pd.read_csv("BUCA_files_3.csv")

igsr = pd.read_csv("/mnt/iribhm/people/nifernan/1000_genome/igsr_samples.tsv", delimiter="\t")
igsr_dict = {name:code for name, code in zip(igsr['Sample name'], igsr['Superpopulation code'])}
cvec = [igsr_dict[Id] for Id in df1['ID']]

print("===== Performing PCA =====")
scaler = MinMaxScaler()
X = scaler.fit_transform(df1.drop('ID', axis=1))
df2 = pd.DataFrame([ast.literal_eval(df2.iloc[i,0]) for i in range(df2.shape[0])])
print(df2.shape)
df3 = [np.array(df2.iloc[i])[:X.shape[1]].reshape(-1, X.shape[1]) for i in range(df2.shape[0])]
df4 = pd.DataFrame(np.squeeze(np.array(df3)))
df = scaler.fit_transform(df4)

# default value works fine for PCA
pca = PCA()
exon_pca = pca.fit_transform(X)
proj = pca.fit_transform(df)

arr = [(i, np.mean(proj[i, :])) for i in range(proj.shape[0])]
print(arr)
x1, y1 =  list(zip(*arr))

N = proj.shape[0]
cvec2 = ['SAMP001']*N

print("===== Plotting =====")
fig1 = px.scatter(exon_pca, x=0, y=1, color=cvec)
fig2 = px.scatter(proj, x=x1, y=y1, color=cvec2, symbol_sequence='x')
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

print("Check the images folder now!!!")