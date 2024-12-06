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

### **Clustering Pipelines**
We implemented and evaluated four clustering pipelines combining dimensionality reduction with clustering algorithms:

#### **Pipeline 1: PCA + K-means Clustering**
- **Dimensionality Reduction:** Principal Component Analysis (PCA) captures the largest variance in the data.
- **Clustering Method:** K-means clustering groups data into clusters by minimizing within-cluster variance.
- **Key Findings:**
  - Computationally efficient with a runtime of **5.56 seconds** and memory usage of **103.78 MB**.
  - Moderate performance with an **ARI of 0.32** and a low silhouette score of **0.10**.
  - Revealed basic trends in the data but struggled with overlapping clusters.

#### **Pipeline 2: PCA + Hierarchical Clustering**
- **Dimensionality Reduction:** PCA simplifies the high-dimensional gene expression data.
- **Clustering Method:** Hierarchical Clustering constructs a dendrogram to identify nested relationships in the data.
- **Key Findings:**
  - Better separation of clusters than K-means alone, with a silhouette score of **0.47**.
  - Computationally more intensive, with a runtime of **30.14 seconds** and memory usage of **104.62 MB**.
  - Highlighted nested relationships, but lacked the resolution needed for biological insights.

#### **Pipeline 3: t-SNE + K-means Clustering**
- **Dimensionality Reduction:** t-SNE captures non-linear patterns in the data, preserving local and global relationships.
- **Clustering Method:** K-means clustering applied to the t-SNE-transformed data.
- **Key Findings:**
  - Better performance in separating clusters, with an **ARI of 0.51** and silhouette score of **0.41**.
  - Runtime of approximately **13 seconds** but memory-intensive, requiring **~1.5 GB**.
  - Enhanced cluster visualization but showed moderate overlap among conditions.

#### **Pipeline 4: t-SNE + Hierarchical Clustering**
- **Dimensionality Reduction:** t-SNE reveals complex, non-linear relationships.
- **Clustering Method:** Hierarchical Clustering constructs a dendrogram to further refine clusters.
- **Key Findings:**
  - The best-performing pipeline with an **ARI of 0.60** and silhouette score of **0.42**.
  - Runtime of **12.80 seconds**, with memory usage of **~1.6 GB**.
  - Provided biologically meaningful clusters that were used for downstream analyses.

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



