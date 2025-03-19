### Pandas
import pandas as pd
from sklearn.datasets import load_iris

# Load Iris dataset
iris = load_iris(as_frame=True)
df = iris.frame

# Display first 5 rows and summary
print(df.head())
print(df.describe())

### Dask
import dask.dataframe as dd
from sklearn.datasets import load_iris

# Load Iris dataset into Dask
iris = load_iris(as_frame=True)
df = dd.from_pandas(iris.frame, npartitions=2)

# Compute mean of features grouped by target
result = df.groupby("target").mean().compute()
print(result)


### Polars
import polars as pl
from sklearn.datasets import load_iris

# Load Iris dataset into Polars
iris = load_iris(as_frame=True)
df = pl.DataFrame(iris.frame)

# Perform a simple filter and selection
print(df.filter(pl.col("sepal length (cm)") > 5.0).select(["sepal length (cm)", "target"]))


### Biopython
from Bio.Seq import Seq

# Short example DNA sequence
seq = Seq("ATGCGTACG")
print("Complement:", seq.complement())
print("Reverse Complement:", seq.reverse_complement())


### PyBigWig
import pyBigWig

# Create and read from a small BigWig file
bw = pyBigWig.open("example.bw", "w")
bw.addHeader([("chr1", 1000)])

# Add a few entries
bw.addEntries(["chr1"], [100], ends=[200], values=[0.5])
bw.close()

# Read values
bw = pyBigWig.open("example.bw")
print(bw.values("chr1", 100, 200))
bw.close()


### Bokeh
from bokeh.plotting import figure, show
from sklearn.datasets import load_iris

# Load Iris dataset
iris = load_iris(as_frame=True)
df = iris.frame

# Scatter plot
plot = figure(title="Iris Dataset: Sepal Length vs. Sepal Width")
plot.scatter(df["sepal length (cm)"], df["sepal width (cm)"], alpha=0.5)
show(plot)


### Holoviews
import holoviews as hv
from sklearn.datasets import load_iris

hv.extension("bokeh")

# Load Iris dataset
iris = load_iris(as_frame=True)
df = iris.frame

# Scatter plot
hv.Scatter(df, "sepal length (cm)", "sepal width (cm)").opts(title="Iris Sepal Dimensions")


### snakemake 
rule all:
    input:
        "processed_iris.csv"

rule process:
    input:
        "iris.csv"
    output:
        "processed_iris.csv"
    shell:
        "tail -n +2 {input} | sort > {output}"


### nextflow 
process filterFASTA {
    input:
        path fasta
    output:
        path "filtered.fasta"

    """
    grep -v "N" ${fasta} > filtered.fasta
    """
}

workflow {
    filterFASTA(fasta: 'example.fasta')
}


### NCBI e-utilities
import requests

# Search NCBI for BRCA1 gene
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
params = {"db": "gene", "term": "BRCA1", "retmax": 5, "retmode": "json"}
response = requests.get(url, params=params)

# Print gene IDs
print(response.json()["esearchresult"]["idlist"])


### BioMart
from bioservices import BioMart

# Use Ensembl BioMart to query human chromosome 1 genes
mart = BioMart(host="www.ensembl.org")

# Construct the XML query
xml_query = """
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="chromosome_name" value="1"/>
    <Attribute name="ensembl_gene_id"/>
  </Dataset>
</Query>
"""

# Query using the XML
result = mart.query(xml_query)

# Display a few genes
print(result.split("\n")[:5])


### scikit-learn
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# Load Iris dataset
X, y = load_iris(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Train and evaluate
model = RandomForestClassifier()
model.fit(X_train, y_train)
print("Test accuracy:", model.score(X_test, y_test))


### XGBoost
import xgboost as xgb
from sklearn.datasets import load_iris

# Load Iris dataset
X, y = load_iris(return_X_y=True)
dtrain = xgb.DMatrix(X, label=y)

# Train XGBoost model
params = {"objective": "multi:softmax", "num_class": 3, "max_depth": 3}
model = xgb.train(params, dtrain, num_boost_round=10)


### SHAP
import shap
import xgboost as xgb
from sklearn.datasets import load_iris

# Train XGBoost model
X, y = load_iris(return_X_y=True)
model = xgb.XGBClassifier().fit(X, y)

# SHAP explanation
explainer = shap.Explainer(model, X)
shap_values = explainer(X)
shap.summary_plot(shap_values, X)


### igraph
import igraph as ig
import matplotlib.pyplot as plt
import random

random.seed(0)
g = ig.Graph.GRG(50, 0.15)

components = g.connected_components(mode='weak')

fig, ax = plt.subplots()
ig.plot(
    components,
    target=ax,
    palette=ig.RainbowPalette(),
    vertex_size=0.07,
    vertex_color=list(map(int, ig.rescale(components.membership, (0, 200), clamp=True))),
    edge_width=0.7
)
plt.show()


### PyVis
from pyvis.network import Network

# Create network
net = Network(notebook=True)
net.add_nodes([1, 2, 3], label=["setosa", "versicolor", "virginica"])
net.add_edges([(1, 2), (2, 3)])
net.show("network.html")

