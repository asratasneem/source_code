# **Streamlining RNA-seq Data Analysis Using Clustering Pipelines**

This repository provides code, results, and documentation for a comparative study on clustering pipelines applied to bulk RNA-seq data. By integrating dimensionality reduction techniques with clustering algorithms, the study offers insights into transcriptomic patterns in COVID-19 patients and healthy controls.

---

## **Project Summary**

RNA sequencing (RNA-seq) provides transformative insights into transcriptomic changes, but analyzing high-dimensional RNA-seq data is computationally challenging. In this project, we compare four clustering pipelines to assess their performance in grouping RNA-seq data and their ability to uncover biologically meaningful insights. The best-performing pipeline was applied to identify differentially expressed genes (DEGs) and enriched pathways, shedding light on molecular mechanisms of COVID-19.

---

## **Key Objectives**

1. Evaluate clustering pipelines combining dimensionality reduction techniques (PCA, t-SNE) with clustering algorithms (K-means, Hierarchical Clustering).
2. Identify biologically significant clusters associated with disease states.
3. Perform differential expression analysis (DEA) and pathway enrichment.

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
### Input Files:
1. **Gene Expression Matrix** (gene_counts.txt):  
   - Contains raw gene counts with metadata (e.g., gene ID, chromosome, start/end positions).  
2. **Metadata File** (SraRunTable.csv):  
   - Includes biological metadata such as sample conditions, severity, and gender.

### Steps:
1. **Filter genes**:  
   - Genes with total counts > 5 across samples were retained.
2. **Log Transformation**:  
   - Log-transformed data using log2(x + 1) for normalization.
3. **Scaling**:  
   - Standardized samples (rows) using StandardScaler to zero-mean and unit variance.
4. **Output**:  
   - The preprocessed dataset was saved as normalized_data.csv for downstream analysis.
  
### **3. Clustering and Evaluation**
- Dimensionality reduction (PCA, t-SNE).
- Clustering (K-means, Hierarchical Clustering).
- Evaluation metrics:
  - **Adjusted Rand Index (ARI):** Measures agreement with true labels.
  - **Silhouette Score:** Reflects cluster separation.
  - **Run time**
  - **Memory usage**

### **4. Differential Expression Analysis (DEA)**
- Identify significant DEGs between clusters.
- Gene Ontology (GO) and KEGG pathway enrichment.

---

## **Results**

## **Clustering Pipelines**

### **Pipeline 1: PCA + K-means Clustering**
- **Description:** 
  - Principal Component Analysis (PCA) was used to reduce dimensionality by capturing global variance trends in the dataset.
  - K-means clustering grouped samples into 3 clusters by minimizing within-cluster variance.
  - Applied PCA to reduce high-dimensional data to 2 components while retaining maximum variance.
- **Implementation Highlights:**
  - Optimal cluster count determined using the elbow method.
  - PCA explained **21.2% of the variance in PC1** and **7.9% in PC2**.
- **Insights:**
  - Showed basic separation of healthy and COVID-19 states but limited biological interpretability.

---

### **Pipeline 2: PCA + Hierarchical Clustering**
- **Description:**
  - PCA was applied to reduce dimensionality.
  - Hierarchical Clustering used Ward’s method to form clusters, visualized with a dendrogram.
- **Implementation Highlights:**
  - Produced nested relationships among clusters.
  - Optimal cluster count (3) chosen from the dendrogram.
- **Insights:**
  - Highlighted nested structure but suffered from overlapping clusters in biological contexts.

---

### **Pipeline 3: t-SNE + K-means Clustering**
- **Description:**
  - t-distributed Stochastic Neighbor Embedding (t-SNE) preserved local and global structures, revealing non-linear relationships in the data.
  - K-means clustering applied to the t-SNE-transformed data.
  - Performed K-means clustering with 3 clusters (k=3) on t-SNE results
- **Implementation Highlights:**
  - Parameters: Perplexity = 30, Iterations = 1000, Random Seed = 42.
  - Clusters visualized in 2D scatterplots.
- **Insights:**
  - Provided a more accurate depiction of cluster separations but required more computational resources.

---

### **Pipeline 4: t-SNE + Hierarchical Clustering**
- **Description:**
  - t-SNE reduced dimensionality, emphasizing local data relationships.
  - Hierarchical Clustering formed nested clusters from t-SNE results, visualized using a dendrogram.
- **Implementation Highlights:**
  - Parameters: Same t-SNE settings as Pipeline 3.
  - Hierarchical clustering extracted 3 optimal clusters.
- **Insights:**
  - The best-performing pipeline with meaningful cluster separations, used for downstream analysis.

---

| **Pipeline**                | **ARI** | **Silhouette Score** | **Runtime (s)** | **Memory Usage** |
|-----------------------------|---------|----------------------|-----------------|------------------|
| PCA + K-means               | 0.32    | 0.10                 | 5.56            | 103.78 MB        |
| PCA + Hierarchical          | 0.24    | 0.47                 | 30.14           | 104.62 MB        |
| t-SNE + K-means             | 0.51    | 0.41                 | 13.00           | 102.81 GB          |
| t-SNE + Hierarchical        | 0.60    | 0.42                 | 12.80           | 103.62 GB          |

---

## **Key Findings**

- **Best Performing Pipeline:** t-SNE + Hierarchical Clustering
  - Achieved the highest clustering agreement (ARI: 0.60) and revealed biologically meaningful clusters.
  - Enabled downstream analysis for identifying DEGs and enriched pathways.

## **Biological Insights**

### **Differential Expression Analysis (DEA)**
- Identified **14,667 significant DEGs**, including:
  - **3,873 upregulated** genes.
  - **10,794 downregulated** genes.
- Notable genes: CYP1B1, RNF216, HERC3 (upregulated), ZNF573, KIAA1279, NEMP1 (downregulated).

### **Gene Ontology (GO) Enrichment**
- Top biological processes:
  - RNA splicing via the spliceosome (most enriched, -log10 p-value > 12).
  - mRNA processing and ribosome biogenesis.

### **KEGG Pathway Enrichment**
- Significant pathways:
  - PD-L1 expression and PD-1 checkpoint pathway.
  - Ubiquitin-mediated proteolysis.
  - Herpes simplex virus 1 infection.

### **Biological Significance**
- SARS-CoV-2 significantly impacts RNA processing, immune modulation, and protein regulation.
- Insights suggest therapeutic targets like RNA splicing and immune checkpoint pathways.

---

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
