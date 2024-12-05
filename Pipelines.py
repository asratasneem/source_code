# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 17:56:48 2024

@author: manas
"""

import pandas as pd
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score
import umap
import time
import tracemalloc

os.chdir('C:/Computational methods_Bohdan')

# Start tracking memory usage and execution time
start_time = time.time()
tracemalloc.start()

# Load the gene counts data from a file
gene_counts = pd.read_csv('gene_counts.txt', sep='\t', comment='#')

# Clean the gene counts data
gene_counts = gene_counts.set_index('Geneid')  # Set Geneid as index
gene_counts = gene_counts.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])  # Drop metadata columns
gene_counts.columns = [col.split('/')[-1].split('_')[0] for col in gene_counts.columns]  # Simplify column names
gene_counts = gene_counts.groupby(by=gene_counts.columns, axis=1).sum()  # Combine duplicate columns

# Load the updated metadata
metadata = pd.read_csv('SraRunTable.csv')

# Ensure alignment of gene counts columns with metadata samples
gene_counts = gene_counts.groupby(by=gene_counts.columns, axis=1).sum()  # Reorder columns to match metadata samples

# Filter out low-expression genes
filtered_counts = gene_counts[gene_counts.sum(axis=1) > 5]

# Log-transform the data
log_counts = np.log2(filtered_counts + 1)

# Normalize the data
scaled_counts = StandardScaler().fit_transform(log_counts.T)  # Transpose for samples as rows

# Convert the scaled data back to a DataFrame
scaled_df = pd.DataFrame(
    scaled_counts,
    index=log_counts.columns,    # Sample names as index
    columns=log_counts.index     # Gene names as columns
)


#print(scaled_df.head())

# Save the normalized data to a CSV file
scaled_df.to_csv('normalized_data.csv', index=True)

#print("Normalized data has been saved to 'normalized_data.csv'.")

#PIPELINE-1

# Perform PCA (Principal Component Analysis) for dimensionality reduction
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_counts)

# Create a DataFrame for the PCA results, adding metadata columns for visualization
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Condition'] = metadata['Condition']
pca_df['Gender'] = metadata['Gender']
pca_df['Severity'] = metadata['Severity']

# Visualize the PCA results based on condition and severity
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Condition', style='Severity', s=100, palette='Set2')
plt.title('PCA: Gene Expression Data by Condition and Severity')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(title='Condition')
plt.show()

# Plot PCA results with a focus on Severity and Gender
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Severity', style='Gender', s=100, palette='Set3')
plt.title('PCA: Gene Expression Data by Severity and Gender')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(title='Severity')
plt.show()


# Perform Hierarchical Clustering
# Use PCA results for hierarchical clustering
Z = linkage(pca_result, method='ward')  # Ward's method minimizes variance within clusters

# Visualize the Dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z, truncate_mode='lastp', p=30, leaf_rotation=90, leaf_font_size=12)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Sample Index')
plt.ylabel('Distance')
plt.show()

# Extract clusters from the dendrogram
num_clusters = 3 # Specify the number of clusters desired
clusters = fcluster(Z, t=num_clusters, criterion='maxclust')  # Cut the tree into the specified number of clusters

# Add cluster labels to the PCA DataFrame for further analysis
pca_df['Cluster'] = clusters

# Calculate Silhouette Score to evaluate clustering quality
silhouette_avg = silhouette_score(pca_result, clusters)
print(f'Silhouette Score for {num_clusters} clusters: {silhouette_avg:.2f}')

#  Visualize Clusters on PCA
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', style='Condition', s=100, palette='tab10')
plt.title('Hierarchical Clustering on PCA Data')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(title='Cluster')
plt.show()

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(scaled_counts)

# Visualize t-SNE clusters
sns.scatterplot(x=tsne_result[:, 0], y=tsne_result[:, 1], hue=clusters, palette='tab10', s=100)
plt.title('t-SNE Visualization of Clusters')
plt.show()

# Perform UMAP for dimensionality reduction
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
umap_result = umap_model.fit_transform(scaled_counts)

# Create UMAP DataFrame
umap_df = pd.DataFrame(umap_result, columns=['UMAP1', 'UMAP2'])
umap_df['Cluster'] = clusters  # Use clusters from hierarchical clustering
umap_df['Condition'] = metadata['Condition']

# Visualize UMAP Clusters
plt.figure(figsize=(8, 6))
sns.scatterplot(data=umap_df, x='UMAP1', y='UMAP2', hue='Cluster', style='Condition', s=100, palette='tab10')
plt.title('UMAP Visualization of Clusters')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))  # Place legend outside the plot
plt.tight_layout()  # Adjust layout to ensure nothing is cut off
plt.show()

# Add cluster labels to log-transformed data
log_counts_transposed = log_counts.T
log_counts_transposed['Cluster'] = clusters

# Perform clustering stability check
n_runs = 100
cluster_memberships = []

for i in range(n_runs):
    Z = linkage(pca_result, method='ward')
    clusters_run = fcluster(Z, t=num_clusters, criterion='maxclust')
    cluster_memberships.append(clusters_run)

ari_scores = []
for i in range(n_runs):
    for j in range(i + 1, n_runs):
        ari_scores.append(adjusted_rand_score(cluster_memberships[i], cluster_memberships[j]))

print(f"Cluster Stability (ARI): Mean = {np.mean(ari_scores):.2f}, Std = {np.std(ari_scores):.2f}")
# Calculate ARI score between predicted clusters and biological condition labels
ari = adjusted_rand_score(metadata['Condition'], clusters)
print(f"Adjusted Rand Index (ARI): {ari:.2f}")
# Cross-tabulate clusters with biological metadata
print("Clusters vs. Condition:")
print(pd.crosstab(metadata['Condition'], clusters))

print("\nClusters vs. Severity:")
print(pd.crosstab(metadata['Severity'], clusters))

print("\nClusters vs. Gender:")
print(pd.crosstab(metadata['Gender'], clusters))

# Stop the timer and print the runtime
end_time = time.time()
runtime = end_time - start_time
print(f"Runtime: {runtime:.2f} seconds")

# Get the current memory usage
current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage: {current / 10**6:.2f} MB")
print(f"Peak memory usage: {peak / 10**6:.2f} MB")

# Stop tracking memory
tracemalloc.stop()

#PIPELINE-2

import pandas as pd
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
import tracemalloc
import time

# Tracking function for runtime and memory
def track_efficiency(func, *args, **kwargs):
    start_time = time.time()
    tracemalloc.start()
    result = func(*args, **kwargs)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    runtime = end_time - start_time
    memory_usage = peak / 10**6  # Convert to MB
    return result, runtime, memory_usage

# Define the entire workflow in a single function
def complete_workflow():
    # Load and preprocess data
    os.chdir('C:/Computational methods_Bohdan')
    gene_counts = pd.read_csv('gene_counts.txt', sep='\t', comment='#')
    gene_counts = gene_counts.set_index('Geneid')
    gene_counts = gene_counts.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
    gene_counts.columns = [col.split('/')[-1].split('_')[0] for col in gene_counts.columns]
    gene_counts = gene_counts.groupby(by=gene_counts.columns, axis=1).sum()
    filtered_counts = gene_counts[gene_counts.sum(axis=1) > 5]
    log_counts = np.log2(filtered_counts + 1)
    scaled_counts = StandardScaler().fit_transform(log_counts.T)
    
    metadata = pd.read_csv('SraRunTable.csv')
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_counts)
    
    # Perform KMeans
    kmeans = KMeans(n_clusters=3, random_state=42)
    clusters = kmeans.fit_predict(scaled_counts)
    
    # Perform t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(scaled_counts)
    
    # Silhouette score
    silhouette_avg = silhouette_score(scaled_counts, clusters)
    print(f"Silhouette Score: {silhouette_avg:.2f}")
    
    # Visualize PCA
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = clusters
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', s=100)
    plt.title('PCA: Clustering Visualization')
    plt.show()
    
    # Visualize t-SNE
    tsne_df = pd.DataFrame(tsne_result, columns=['tSNE1', 'tSNE2'])
    tsne_df['Cluster'] = clusters
    sns.scatterplot(data=tsne_df, x='tSNE1', y='tSNE2', hue='Cluster', palette='tab10', s=100)
    plt.title('t-SNE: Clustering Visualization')
    plt.show()

# Measure the overall workflow performance
_, total_time, total_memory = track_efficiency(complete_workflow)

print(f"Overall Runtime: {total_time:.2f} seconds")
print(f"Overall Memory Usage: {total_memory:.2f} MB")

#Pipeline 3: t-SNE + K-means
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
tsne_result = tsne.fit_transform(scaled_counts)

# Create a DataFrame for t-SNE results
tsne_df = pd.DataFrame(tsne_result, columns=['tSNE1', 'tSNE2'])

# Perform K-means clustering on t-SNE results
kmeans_tsne = KMeans(n_clusters=3, random_state=42)
tsne_clusters = kmeans_tsne.fit_predict(tsne_result)

# Add clusters and metadata to the t-SNE DataFrame
tsne_df['Cluster'] = tsne_clusters
tsne_df['Condition'] = metadata['Condition']

# Visualize t-SNE clusters
plt.figure(figsize=(8, 6))
sns.scatterplot(data=tsne_df, x='tSNE1', y='tSNE2', hue='Cluster', style='Condition', palette='tab10', s=100)
plt.title('t-SNE + K-means Clustering')
plt.legend(title='Cluster')
plt.show()

# Calculate silhouette score for t-SNE + K-means
silhouette_tsne = silhouette_score(tsne_result, tsne_clusters)
print(f"Silhouette Score for t-SNE + K-means: {silhouette_tsne:.2f}")

#Pipeline 4: t-SNE + Hierarchical Clustering

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
tsne_result = tsne.fit_transform(scaled_counts)

# Create a DataFrame for t-SNE results
tsne_hc_df = pd.DataFrame(tsne_result, columns=['tSNE1', 'tSNE2'])

# Perform hierarchical clustering on t-SNE results
linkage_matrix = linkage(tsne_result, method='ward')
hc_clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')  # 3 clusters

# Add clusters and metadata to the t-SNE DataFrame
tsne_hc_df['Cluster'] = hc_clusters
tsne_hc_df['Condition'] = metadata['Condition']

# Visualize hierarchical clustering on t-SNE results
plt.figure(figsize=(8, 6))
sns.scatterplot(data=tsne_hc_df, x='tSNE1', y='tSNE2', hue='Cluster', style='Condition', palette='tab10', s=100)
plt.title('t-SNE + Hierarchical Clustering')
plt.legend(title='Cluster')
plt.show()

# Visualize the dendrogram
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Samples')
plt.ylabel('Distance')
plt.show()

# Calculate Silhouette Score
silhouette_hc = silhouette_score(tsne_result, hc_clusters)
print(f"Silhouette Score for t-SNE + Hierarchical Clustering: {silhouette_hc:.2f}")


