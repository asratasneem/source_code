# Computational-Biology----Project-Report
Group Information
Group Number: 09
Members: Asra Tasneem Shaik, Muni Manasa Vema, Mahima Mahabaleshwar Siddheshwar, Saranya Guvvala
Course: Computational Methods for Biomedical Informatics (B536)
Project Title: Streamlining RNA-seq Data Analysis: A Framework for Comparing Computational Pipelines

# Scope
This project evaluates and compares RNA-seq pipelines for analyzing gene expression variability in COVID-19 patients.
Objective: Compare four pipelines using PCA + K-means, PCA + Hierarchical Clustering, t-SNE + K-means, and t-SNE + Hierarchical Clustering. Perform Differential Expression Analysis (DEA) on clusters from the best pipeline.

# Datasets & Dependencies
Disease States: 2 Convalescent, 32 COVID-19, 34 Healthy samples (GEO Accession: GSE152418).
Reference Genome: UCSC hg38.
Tools/Dependencies:
Preprocessing: SRA Toolkit, FastQC, Trim Galore, HISAT2, SAMtools, featureCounts.
Python Packages: pandas, numpy, scikit-learn, matplotlib, seaborn, DESeq2, ggplot2.

## Workflow
# Data Preprocessing
Download Data: Used prefetch to fetch RNA-seq data and fasterq-dump for FASTQ conversion.
Quality Control: Conducted with FastQC on raw and trimmed reads.
Trimming: Trim Galore was used for adapter trimming.
# Alignment & Quantification
Reference Preparation: Built HISAT2 index for UCSC hg38 genome.
Read Alignment: Aligned reads using HISAT2 and generated sorted BAM files with SAMtools.
Gene Quantification: Used featureCounts to produce gene count matrices.
# 1. Data Preparation
Gene Counts Cleaning:
Drops unnecessary metadata columns (e.g., Chr, Start, End).
Groups duplicate sample columns and filters out low-expression genes.
Log Transformation and Scaling:
Log2 transformation is applied to normalize the data.
StandardScaler ensures all features have a mean of 0 and a standard deviation of 1.
# Pipeline 1: PCA + K-means
PCA reduced data to 2 dimensions for clustering.
K-means clustering (k=3) identified clusters, evaluated with Silhouette Score and ARI.
# Pipeline 2: PCA + Hierarchical Clustering
PCA reduced dimensionality.
Hierarchical clustering (Ward's method) formed 3 clusters visualized in a dendrogram.
Clusters were evaluated.
# Pipeline 3: t-SNE + K-means
t-SNE reduced data to 2 dimensions with parameters (perplexity=30, iterations=1000).
K-means clustering (k=3) was applied, and clusters evaluated with Silhouette Score and ARI.
# Pipeline 4: t-SNE + Hierarchical Clustering
t-SNE for dimensionality reduction (same parameters as Pipeline 3).
Hierarchical clustering (Ward's method) generated clusters and a dendrogram.
Evaluated using Silhouette Score and ARI, followed by DEG identification.


# Key Outputs
Preprocessed Data: Quality-trimmed FASTQ files and gene count matrices.
Performance Metrics:
Runtime, memory usage, Silhouette Score, ARI.
Visualizations:
Elbow method plots, PCA/t-SNE scatterplots, dendrogram.
Differential Expression Analysis:
Identifies DEGs across clusters using ANOVA.
Outputs significant genes for further biological interpretation.

# Execution Notes
Environment Setup: Modules (sra-toolkit, HISAT2, FastQC, subread) and conda environments are required.
Scripts: Provided for data processing, alignment, quantification, and pipeline evaluation. Ensure paths are updated before execution.
Data Locations:
Gene Count File: gene_counts.txt
Scripts: .py or .ipynb files

#File Paths
Raw Data Directory: /N/slate/msiddhe/Computational_Project/Raw_Data
FASTQ Directory: /N/slate/msiddhe/Computational_Project/FASTQ_Data
Trimmed Data Directory: /N/slate/msiddhe/Computational_Project/Trimmed_Data
Output Directories: Alignment, FastQC, and featureCounts directories.

# Next Steps: Differential Gene Expression Analysis
Pipeline Selection
The t-SNE + Hierarchical Clustering pipeline was identified as the most effective method for clustering and visualizing gene expression patterns in COVID-19 samples. This pipeline forms the basis for downstream analyses.
Differential Gene Expression (DEG) Analysis

# Key Findings:
Identified 14,667 significant genes, including:
3,873 upregulated genes (e.g., CYP1B1, RNF216, HERC3, POLR2B).
10,794 downregulated genes (e.g., ZNF573, KIAA1279, NEMP1).
Heatmap of the top 50 DEGs showed clear expression clusters with distinct upregulation and downregulation patterns.
Volcano plot illustrated strong statistical significance, with -log10 p-values up to 25.
Figures:

Figure 14: Volcano Plot of Differentially Expressed Genes.
Figure 15: Heatmap of Top 50 Differentially Expressed Genes.
Gene Ontology (GO) Enrichment Analysis

# Top Enriched Biological Processes:
mRNA splicing via the spliceosome (-log10 adjusted p-value > 12).
mRNA processing, RNA splicing via transesterification reactions, and ribosome biogenesis.
Emphasizes the role of RNA processing in COVID-19 pathogenesis.
Figure 16: Top 10 Enriched Biological Processes from GO Analysis.

# KEGG Pathway Analysis
Significant Pathways Identified:
Herpes simplex virus 1 infection.
Acute myeloid leukemia.
Colorectal cancer.
Ubiquitin-mediated proteolysis.
PD-L1 expression and PD-1 checkpoint pathway.
Highlights molecular mechanisms related to immune response modulation and cellular regulatory processes.
Figure 17: Top 10 Enriched Pathways from KEGG Analysis.

# Insights and Importance
Molecular Mechanisms: RNA processing plays a crucial role in the disease mechanism of COVID-19.
Immune Modulation: Pathways such as PD-1 checkpoint emphasize immune regulation in the disease.
Disease Relevance: Enrichment in pathways related to cancer and viral infections suggests overlapping mechanisms.


