Below is a README file based on the provided report and Python clustering pipelines code:

---

# Streamlining RNA-seq Data Analysis: A Framework for Comparing Computational Pipelines

## Project Overview

This repository provides the methodology, scripts, and results from our study evaluating clustering pipelines for RNA-seq data analysis. The project focuses on analyzing bulk RNA-seq data from COVID-19 patients and healthy controls using a combination of dimensionality reduction techniques and clustering algorithms.

### Key Objectives:
1. Evaluate clustering performance using various metrics (e.g., silhouette score, ARI).
2. Compare computational efficiency across pipelines.
3. Identify biologically meaningful patterns and perform downstream analysis like differential gene expression.

## Dataset

**Source:** [Gene Expression Omnibus (GSE152418)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)  
**Composition:** 68 samples - 2 Convalescent, 32 COVID-19, 34 Healthy  
**Data Type:** Bulk RNA sequencing (RNA-seq)

## Workflow Overview

1. **Data Preprocessing**
   - Download and quality check raw FASTQ files.
   - Trim reads and confirm quality improvement.
   - Align reads to the human reference genome (`hg38`) using HISAT2.
   - Generate gene-level count matrix using featureCounts.

2. **Dimensionality Reduction and Clustering**
   - PCA: For global variance simplification.
   - t-SNE: To uncover non-linear relationships.
   - UMAP: For fine-grained local/global structure visualization.
   - Clustering Algorithms: K-means and Hierarchical Clustering.

3. **Evaluation Metrics**
   - Silhouette Score
   - Adjusted Rand Index (ARI)
   - Computational metrics: runtime and memory usage.

4. **Biological Insights**
   - Differential gene expression analysis.
   - Gene Ontology (GO) enrichment.
   - Pathway analysis (KEGG).

5. **Best Pipeline**
   - **t-SNE + Hierarchical Clustering** was identified as the optimal approach.

## Repository Structure

- **`Clustering_Pipelines.py`**  
   Contains scripts for dimensionality reduction, clustering, and evaluation. Pipelines include:
   - PCA + K-means
   - PCA + Hierarchical Clustering
   - t-SNE + K-means
   - t-SNE + Hierarchical Clustering

- **`normalized_data.csv`**  
   Preprocessed and log-transformed RNA-seq data used for analysis.

- **`gene_counts.txt`**  
   Raw gene count matrix.

- **Figures and Results**  
   Visualizations of clustering results (PCA, t-SNE, UMAP) and metrics.

## Results

| Pipeline                 | ARI  | Silhouette Score | Runtime (s) | Memory Usage |
|--------------------------|------|------------------|-------------|--------------|
| PCA + K-means            | 0.32 | 0.10             | 5.56        | 103.78 MB    |
| PCA + Hierarchical       | 0.24 | 0.47             | 30.14       | 104.62 MB    |
| t-SNE + K-means          | 0.51 | 0.41             | 13.00       | ~1.5 GB      |
| t-SNE + Hierarchical     | 0.60 | 0.42             | 12.80       | ~1.6 GB      |

### Biological Findings
- **Differential Expression Analysis:** Identified 14,667 significant genes (3,873 upregulated, 10,794 downregulated).  
- **Top GO Terms:** RNA splicing, ribosome biogenesis.  
- **KEGG Pathways:** PD-L1 expression, ubiquitin-mediated proteolysis.

## Usage Instructions

### Prerequisites
- Python 3.8+
- Required Python libraries: `pandas`, `numpy`, `scikit-learn`, `seaborn`, `matplotlib`, `umap-learn`, `scipy`.

### Running the Pipelines
1. Place `gene_counts.txt` and `SraRunTable.csv` in the working directory.
2. Update the working directory in the script (`os.chdir()`).
3. Run the script:
   ```bash
   python Clustering_Pipelines.py
   ```

### Outputs
- **Plots:** Visualizations for PCA, t-SNE, UMAP, dendrograms, and cluster scatterplots.
- **Metrics:** Printed ARI, silhouette scores, runtime, and memory usage.

## License
This project is licensed under the MIT License.

---

Let me know if you'd like to customize any section or include additional details!