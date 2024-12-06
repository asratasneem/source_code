
# **Streamlining RNA-seq Data Analysis Using Clustering Pipelines**

This repository provides code, results, and documentation for a comparative study on clustering pipelines applied to bulk RNA-seq data. By integrating dimensionality reduction techniques with clustering algorithms, the study offers insights into transcriptomic patterns in COVID-19 patients and healthy controls.

---

## **Project Summary**

RNA sequencing (RNA-seq) provides transformative insights into transcriptomic changes, but analyzing high-dimensional RNA-seq data is computationally challenging. In this project, we compare four clustering pipelines to assess their performance in grouping RNA-seq data and their ability to uncover biologically meaningful insights. The best-performing pipeline was applied to identify differentially expressed genes (DEGs) and enriched pathways, shedding light on molecular mechanisms of COVID-19.

---

## **Key Objectives**

1. Evaluate clustering pipelines combining dimensionality reduction techniques (PCA, t-SNE) with clustering algorithms (K-means, Hierarchical Clustering).
2. Identify biologically significant clusters associated with disease states.
3. Perform differential expression analysis (DEA) and pathway enrichment for actionable insights.

---

## **Dataset Overview**

- **Source:** [GEO Accession: GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)  
- **Samples:** 68 total (34 Healthy, 32 COVID-19, 2 Convalescent)
- **Data Type:** Bulk RNA-seq
- **Relevance:** Investigates the transcriptomic landscape of SARS-CoV-2 infection, recovery, and host immune response.

---

## **Pipelines**

### **Pipeline 1: PCA + K-means Clustering**
- **Strengths:** Computationally efficient, captures global variance trends.
- **Weaknesses:** Limited ability to separate biologically meaningful clusters.

### **Pipeline 2: PCA + Hierarchical Clustering**
- **Strengths:** Identifies nested relationships; moderate silhouette score.
- **Weaknesses:** Computationally slower than Pipeline 1.

### **Pipeline 3: t-SNE + K-means Clustering**
- **Strengths:** Excellent visualization of non-linear relationships; moderate ARI.
- **Weaknesses:** Memory-intensive and slower than PCA-based methods.

### **Pipeline 4: t-SNE + Hierarchical Clustering**
- **Best Pipeline:** Combines the strengths of t-SNE for non-linear data relationships and hierarchical clustering for nested insights.
- **Performance:** Highest ARI (0.60) and a good silhouette score (0.42).

---

## **Workflow**

### **1. Preprocessing**
- Quality control and read trimming using FastQC and TrimGalore.
- Alignment to the human genome (`hg38`) with HISAT2.
- Gene count matrix generation using featureCounts.

### **2. Data Normalization**
- Log transformation of gene expression data.
- Filtering low-expression genes (counts > 5).
- Standardization of data for clustering.

### **3. Clustering and Evaluation**
- Dimensionality reduction (PCA, t-SNE).
- Clustering (K-means, Hierarchical Clustering).
- Evaluation metrics:
  - **Adjusted Rand Index (ARI):** Measures agreement with true labels.
  - **Silhouette Score:** Reflects cluster separation.

### **4. Differential Expression Analysis (DEA)**
- Identify significant DEGs between clusters.
- Gene Ontology (GO) and KEGG pathway enrichment.

---

## **Results**

| **Pipeline**            | **ARI** | **Silhouette Score** | **Runtime (s)** | **Memory Usage** |
|--------------------------|---------|----------------------|-----------------|------------------|
| PCA + K-means            | 0.32    | 0.10                 | 5.56            | 103.78 MB        |
| PCA + Hierarchical       | 0.24    | 0.47                 | 30.14           | 104.62 MB        |
| t-SNE + K-means          | 0.51    | 0.41                 | 13.00           | ~1.5 GB          |
| **t-SNE + Hierarchical** | **0.60**| **0.42**             | **12.80**       | **~1.6 GB**      |

### **Biological Insights**
- **DEGs Identified:** 14,667 genes (3,873 upregulated, 10,794 downregulated).
- **Top Pathways:** RNA splicing, ribosome biogenesis, PD-L1 expression.
- **Key Genes:** CYP1B1, HERC3, RNF216.

---

## **Repository Structure**

- **`Clustering_Pipelines.py`**  
   Python script implementing all four pipelines with visualization and metrics.

- **`gene_counts.txt`**  
   Raw gene count matrix used for analysis.

- **`normalized_data.csv`**  
   Preprocessed and normalized data.

- **Figures**  
   Includes PCA, t-SNE, UMAP visualizations, dendrograms, and DEG heatmaps.

---

## **How to Run**

### **Dependencies**
- Python 3.8+
- Libraries: `pandas`, `numpy`, `scikit-learn`, `matplotlib`, `seaborn`, `umap-learn`, `scipy`.

### **Instructions**
1. Clone the repository:
   ```bash
   git clone https://github.com/username/repo_name.git
   ```
2. Place the following files in the working directory:
   - `gene_counts.txt`
   - `SraRunTable.csv`
3. Run the clustering pipeline:
   ```bash
   python Clustering_Pipelines.py
   ```
4. Outputs:
   - Visualizations: PCA, t-SNE, UMAP, dendrograms.
   - Metrics: Silhouette scores, ARI, runtime, memory usage.

---

## **Acknowledgments**

**Team Members:**  
Asra Tasneem Shaik, Muni Manasa Vema, Mahima Mahabaleshwar Siddheshwar, Saranya Guvvala.  

**Course:** Computational Methods for Biomedical Informatics (B536).  

---
