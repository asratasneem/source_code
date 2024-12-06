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

## **Clustering Pipelines**

### **Pipeline 1: PCA + K-means Clustering**
- **Description:** 
  - Principal Component Analysis (PCA) was used to reduce dimensionality by capturing global variance trends in the dataset.
  - K-means clustering grouped samples into 3 clusters by minimizing within-cluster variance.
  - Applied PCA to reduce high-dimensional data to 2 components while retaining maximum variance.
- **Implementation Highlights:**
  - Optimal cluster count determined using the elbow method.
  - PCA explained **21.2% of the variance in PC1** and **7.9% in PC2**.
- **Evaluation Results:**
  - **Runtime:** 103.07 seconds
  - **ARI:** 0.32
  - **Silhouette Score:** 0.10
  - **Memory usage:** 218.43 MB
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
- **Evaluation Results:**
  - **Runtime:** 30.14 seconds
  - **ARI:** 0.24
  - **Silhouette Score:** 0.47
  - **Memory usage:** 104.62 MB
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
- **Evaluation Results:**
  - **Runtime:** 13 seconds
  - **ARI:** 0.51
  - **Silhouette Score:** 0.41
  - **Memory usage:** 102.81 MB
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
- **Evaluation Results:**
  - **Runtime:** 12.8 seconds
  - **ARI:** 0.60
  - **Silhouette Score:** 0.42
  - **Memory usage:** 103.62 MB
- **Insights:**
  - The best-performing pipeline with meaningful cluster separations, used for downstream analysis.

---

## **Key Findings**

- **Best Performing Pipeline:** t-SNE + Hierarchical Clustering
  - Achieved the highest clustering agreement (ARI: 0.60) and revealed biologically meaningful clusters.
  - Enabled downstream analysis for identifying DEGs and enriched pathways.

- **Biological Insights:**
  - Top GO Terms: RNA splicing, ribosome biogenesis.
  - Significant Pathways: PD-L1 expression, ubiquitin-mediated proteolysis.
  - Highlighted key genes: CYP1B1, HERC3, RNF216.

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












