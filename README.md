# Olink-based-proteomics-analysis
This repository contains a case study of an Olink-based proteomics analysis pipeline. The analysis was performed on an online dataset of serum from allergic patients. The objective of the analysis is to identify the differentially expressed proteins across different levels of allergy. Starting from Olink NPX data, differential expression analysis and functional enrichment analysis were performed in R.
 
The code is divided into 3 main sections:
1. Explorative analysis: PCA, t-SNE, Heatmap with comments
2. Differential expression analysis: statistical analysis using limma
3. Functional enrichment analysis: ORA (Over-Representation Analysis)
 
The results are constantly commented throughout the code, and visualized through a series of graphs, that can be found in the pdf "Output Graph.pdf", within the repository.
 
## Experimental Design and Dataset
The dataset is made of serum samples from patients with mild and severe allergy and healthy controls (three conditions: "Control", "Mild", "Severe"). These samples were processed for Olink Target 96 analysis on three protein panels:
- Cardiometabolic (Cm) panel: 92 target proteins on 29 samples (10 "Control", 9 "Mild", 10 "Severe")
- Cardiovascular II (Cv) panel: 92 target proteins on 29 samples (same as Cm)
- Immuno-Oncology (Io) panel: 92 target proteins on 30 samples (same as Cm, with one additional "Severe" sample. This sample was excluded for consistent comparison across the panels)
 
Original dataset on Mendeley: https://doi.org/10.17632/9b64psp627.1
 
## Analysis Results
Through limma statistical analysis, 28 protein-contrast pairs were found significant at raw p-value < 0.05 (across 22 unique proteins). Functional enrichment ORA on these did not show any significant GO term. In conclusion, results show that only a mild difference in protein expression is present among the different conditions.
 
## R Package Versions
| Package | Version |
|---|---|
| readxl | 1.4.5 |
| reshape2 | 1.4.5 |
| irlba | 2.3.7 |
| factoextra | 1.0.7 |
| ggplot2 | 4.0.2 |
| cowplot | 1.2.0 |
| Rtsne | 0.17 |
| pheatmap | 1.0.13 |
| tidyverse | 2.0.0 |
| car | 3.1.5 |
| ggpubr | 0.6.2 |
| limma | 3.64.3 |
| clusterProfiler | 4.16.0 |
| org.Hs.eg.db | 3.21.0 |
| enrichplot | 1.28.4 |
