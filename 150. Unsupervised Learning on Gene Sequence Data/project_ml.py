## Required Libraries
import umap  # for umap method of dimensionality reduction
import numpy as np  # arithmetic operations
import pandas as pd  # handling dataframe objects
from sklearn.decomposition import PCA  # for pca method of dimensionality reduction
from sklearn.manifold import TSNE  # for t-sne method of dimensionality reduction
from sklearn.cluster import KMeans  # KMeans clustering
import matplotlib.pyplot as plt  # for plotting data


## 1. Read and parse the gene sequence data based on the provided genes

dataset = pd.read_csv("geneExpFinal.csv")
#assigned_gene_indices = dataset[dataset['feature'].isin(genes)].index
dataset.set_index('feature', inplace=True)
# list of assigned genes
genes = ['Acads', 'Dusp28', 'Emx1', 'Vip', 'Urah']
# Extract genes from the dataset and transpose it
df = dataset.loc[genes].T

## 2. Clustering on All Genes

# Perform elbow method
values_of_K = range(2, 6)
scores = []
# Loop through K values
for k in values_of_K:
    # Fit the model
    model = KMeans(n_clusters=k, n_init=10, random_state=42)
    labels = model.fit_predict(df)

    # Calculate fitness statistic
    score = model.inertia_
    scores.append(score)

# Plot the results
# Create figure
plt.figure(figsize=(8, 6))
plt.plot(values_of_K, scores, marker='o')
plt.title(f'Elbow Method For KMeans Clustering')
plt.xlabel('Number of Clusters (K)')
plt.ylabel('Inertia')
plt.savefig("Elbow method for all genes.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()

# Fit the KMeans model
n_clusters = 4
model = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels = model.fit_predict(df)
centroids = model.cluster_centers_

# Visualize cluster assignments and centroids
plt.figure(figsize=(10, 6))
# Scatter plot for the original data points
plt.scatter(x=df.iloc[:, 0], y=df.iloc[:, 1], c=labels, cmap='viridis', s=100)
plt.title(f'KMeans Clustering - Cluster Assignments (K={n_clusters})')
plt.xlabel('Gene - Acads')
plt.ylabel('Gene - Dusp28')
# Scatter plot for the centroids
plt.scatter(centroids[:, 0], centroids[:, 1], c='red', marker='*', s=200, label='Centroids')
plt.legend()
plt.savefig("Clustering for all genes.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()


## 3. Clustering for Dimensionality Reduced Data

# Perform PCA
pca = PCA(n_components=2, random_state=42)
pca_result = pca.fit_transform(df)

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(df)

# Perform UMAP
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_result = umap_model.fit_transform(df)

# create subplots
fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharex=True, sharey=True)

# Visualize clustering based on all 5 genes
axes[0].scatter(pca_result[:, 0], pca_result[:, 1], cmap='viridis', alpha=0.5)
axes[0].set_title('PCA')

# Visualize clustering using reduced dimensions
axes[1].scatter(tsne_result[:, 0], tsne_result[:, 1], cmap='viridis', alpha=0.5)
axes[1].set_title('TSNE')

# Visualize clustering using reduced dimensions
axes[2].scatter(umap_result[:, 0], umap_result[:, 1], cmap='viridis', alpha=0.5)
axes[2].set_title('UMAP')

plt.savefig("Dimensionality reduction of 5 genes using pca tsne and umap.png", dpi=300, format="PNG")
plt.show()

# Perform elbow method for reduced data
values_of_K = range(2, 6)
dr_results = [pca_result, tsne_result, umap_result]
dr_methods = ['PCA', 'TSNE', 'UMAP']

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(12, 6))

for i, result in enumerate(dr_results):
    scores = []
    # Loop through K values
    for k in values_of_K:
        # Fit the model
        model = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = model.fit_predict(result)

        # Calculate fitness statistic
        score = model.inertia_
        scores.append(score)

    # Plot the results
    axs[i].plot(values_of_K, scores, marker='o')
    axs[i].set_title(f'Elbow Method For {dr_methods[i]}')
    axs[i].set_xlabel('Number of Clusters (K)')
    axs[i].set_ylabel('Inertia')

plt.savefig("Elbow method for reduced genes from 5.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()

# Fit the KMeans Clustering model
n_clusters = 4
model_pca = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels = model_pca.fit_predict(pca_result)
centroids_pca = model_pca.cluster_centers_

# Perform KMeans clustering on t-SNE result
model_tsne = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_tsne = model_tsne.fit_predict(tsne_result)
centroids_tsne = model_tsne.cluster_centers_

# Perform KMeans clustering on UMAP result
model_umap = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_umap = model_umap.fit_predict(umap_result)
centroids_umap = model_umap.cluster_centers_

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 6), sharex=True, sharey=True)

# Visualize cluster assignments for PCA result
axs[0].scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='viridis', s=100)
axs[0].scatter(centroids_pca[:, 0], centroids_pca[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[0].set_title(f'KMeans Clustering - PCA (K={n_clusters})')
axs[0].set_xlabel('Principal Component 1 (PCA1)')
axs[0].set_ylabel('Principal Component 2 (PCA2)')
axs[0].legend()

# Visualize cluster assignments for t-SNE result
axs[1].scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels_tsne, cmap='viridis', s=100)
axs[1].scatter(centroids_tsne[:, 0], centroids_tsne[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[1].set_title(f'KMeans Clustering - t-SNE (K={n_clusters})')
axs[1].set_xlabel('t-SNE Dimension 1')
axs[1].set_ylabel('t-SNE Dimension 2')
axs[1].legend()

# Visualize cluster assignments for UMAP result
axs[2].scatter(umap_result[:, 0], umap_result[:, 1], c=labels_umap, cmap='viridis', s=100)
axs[2].scatter(centroids_umap[:, 0], centroids_umap[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[2].set_title(f'KMeans Clustering - UMAP (K={n_clusters})')
axs[2].set_xlabel('UMAP Dimension 1')
axs[2].set_ylabel('UMAP Dimension 2')
axs[2].legend()

plt.savefig("Clustering for reduced genes from 5.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()


# ### 3.1 Bonus

# Fit the KMeans Clustering model for K=2
n_clusters = 2
model_pca = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels = model_pca.fit_predict(pca_result)
centroids_pca = model_pca.cluster_centers_

# Perform KMeans clustering on t-SNE result
model_tsne = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_tsne = model_tsne.fit_predict(tsne_result)
centroids_tsne = model_tsne.cluster_centers_

# Perform KMeans clustering on UMAP result
model_umap = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_umap = model_umap.fit_predict(umap_result)
centroids_umap = model_umap.cluster_centers_

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 6), sharex=True, sharey=True)

# Visualize cluster assignments for PCA result
axs[0].scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='viridis', s=100)
axs[0].scatter(centroids_pca[:, 0], centroids_pca[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[0].set_title(f'KMeans Clustering - PCA (K={n_clusters})')
axs[0].set_xlabel('Principal Component 1 (PCA1)')
axs[0].set_ylabel('Principal Component 2 (PCA2)')
axs[0].legend()

# Visualize cluster assignments for t-SNE result
axs[1].scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels_tsne, cmap='viridis', s=100)
axs[1].scatter(centroids_tsne[:, 0], centroids_tsne[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[1].set_title(f'KMeans Clustering - t-SNE (K={n_clusters})')
axs[1].set_xlabel('t-SNE Dimension 1')
axs[1].set_ylabel('t-SNE Dimension 2')
axs[1].legend()

# Visualize cluster assignments for UMAP result
axs[2].scatter(umap_result[:, 0], umap_result[:, 1], c=labels_umap, cmap='viridis', s=100)
axs[2].scatter(centroids_umap[:, 0], centroids_umap[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[2].set_title(f'KMeans Clustering - UMAP (K={n_clusters})')
axs[2].set_xlabel('UMAP Dimension 1')
axs[2].set_ylabel('UMAP Dimension 2')
axs[2].legend()

plt.savefig("Clustering for K=2 for reduced genes.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()

## 4. Clustering for Top 3 Genes

# Calculate the Pearson correlation coefficient for each pair of genes
gene_corr = np.corrcoef(dataset, rowvar=True)
# create dataframe from gene_corr
corr_df = pd.DataFrame(gene_corr, index=dataset.index, columns=dataset.index)
# find out top 3 genes for each assigned genes
top_gene_indices = []
for gene in genes:
    top_3_closest_genes = corr_df[gene].abs().sort_values(ascending=False)[1:4].index.tolist()
    top_gene_indices.extend(top_3_closest_genes)

# create dataframe with all 20 genes
top_genes_df = dataset.loc[top_gene_indices].T
top_genes_df = pd.concat([df, top_genes_df], axis=1)

# Perform PCA
pca = PCA(n_components=2, random_state=42)
pca_result = pca.fit_transform(top_genes_df)

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(top_genes_df)

# Perform UMAP
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_result = umap_model.fit_transform(top_genes_df)

fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharex=True, sharey=True)

# Visualize clustering based on all 5 genes
axes[0].scatter(pca_result[:, 0], pca_result[:, 1], cmap='viridis', alpha=0.5)
axes[0].set_title('PCA')

# Visualize clustering using reduced dimensions
axes[1].scatter(tsne_result[:, 0], tsne_result[:, 1], cmap='viridis', alpha=0.5)
axes[1].set_title('TSNE')

# Visualize clustering using reduced dimensions
axes[2].scatter(umap_result[:, 0], umap_result[:, 1], cmap='viridis', alpha=0.5)
axes[2].set_title('UMAP')

plt.savefig("Dimensionality reduction of 20 genes using pca tsne and umap.png", dpi=300, format="PNG")
plt.show()

# Perform elbow method for reduced data from 20 genes
values_of_K = range(2, 6)
dr_results = [pca_result, tsne_result, umap_result]
dr_methods = ['PCA', 'TSNE', 'UMAP']

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(12, 6))

for i, result in enumerate(dr_results):
    scores = []
    # Loop through K values
    for k in values_of_K:
        # Fit the model
        model = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = model.fit_predict(result)

        # Calculate fitness statistic
        score = model.inertia_
        scores.append(score)

    # Plot the results
    axs[i].plot(values_of_K, scores, marker='o')
    axs[i].set_title(f'Elbow Method For {dr_methods[i]}')
    axs[i].set_xlabel('Number of Clusters (K)')
    axs[i].set_ylabel('Inertia')

plt.savefig("Elbow method for reduced genes from 20.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()

# Fit the KMeans Clustering model
n_clusters = 4
model_pca = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels = model_pca.fit_predict(pca_result)
centroids_pca = model_pca.cluster_centers_

# Perform KMeans clustering on t-SNE result
model_tsne = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_tsne = model_tsne.fit_predict(tsne_result)
centroids_tsne = model_tsne.cluster_centers_

# Perform KMeans clustering on UMAP result
model_umap = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
labels_umap = model_umap.fit_predict(umap_result)
centroids_umap = model_umap.cluster_centers_

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 6), sharex=True, sharey=True)

# Visualize cluster assignments for PCA result
axs[0].scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='viridis', s=100)
axs[0].scatter(centroids_pca[:, 0], centroids_pca[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[0].set_title(f'KMeans Clustering - PCA (K={n_clusters})')
axs[0].set_xlabel('Principal Component 1 (PCA1)')
axs[0].set_ylabel('Principal Component 2 (PCA2)')
axs[0].legend()

# Visualize cluster assignments for t-SNE result
axs[1].scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels_tsne, cmap='viridis', s=100)
axs[1].scatter(centroids_tsne[:, 0], centroids_tsne[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[1].set_title(f'KMeans Clustering - t-SNE (K={n_clusters})')
axs[1].set_xlabel('t-SNE Dimension 1')
axs[1].set_ylabel('t-SNE Dimension 2')
axs[1].legend()

# Visualize cluster assignments for UMAP result
axs[2].scatter(umap_result[:, 0], umap_result[:, 1], c=labels_umap, cmap='viridis', s=100)
axs[2].scatter(centroids_umap[:, 0], centroids_umap[:, 1], c='red', marker='*', s=200, label='Centroids')
axs[2].set_title(f'KMeans Clustering - UMAP (K={n_clusters})')
axs[2].set_xlabel('UMAP Dimension 1')
axs[2].set_ylabel('UMAP Dimension 2')
axs[2].legend()

plt.savefig("Clustering for reduced genes from 20.png", dpi=300, format="PNG")
plt.tight_layout()
plt.show()
