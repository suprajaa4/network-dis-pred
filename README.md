# Network-Based Disease Gene Prediction

Predict disease-associated genes using protein–protein interaction (PPI) networks, random-walk proximity, and graph-topological features.  
This pipeline ranks candidate genes for any disease given a seed set of known associated genes.

---

## Overview

This project implements a complete pipeline for network-based disease gene discovery:

1. **Fetch PPI network** from STRING for seed genes  
2. **Build induced graph** (largest connected component)  
3. **Compute node features**  
   - degree  
   - betweenness  
   - closeness  
   - PageRank  
   - random-walk with restart (RWR) using seeds  
4. **Rank all genes** via z-score aggregation  
5. **Validate** via held-out seeds  
   - AUROC  
   - Recall@50 / Recall@100  

Output includes a ranked list of top 100 candidate genes and validation metrics.

---

##  Data

This project uses publicly available biological resources.

STRING Database
Szklarczyk et al. (2023). STRING v12: protein–protein association networks.

DisGeNET
Piñero et al. (2020). The DisGeNET knowledge platform for disease genomics.

##  Example: Alzheimer’s Disease

Default seeds (APP, PSEN1, APOE, MAPT…) are bundled.  
You can supply any disease-gene list (e.g., DisGeNET).

---

## Installation

Requires R ≥ 4.2.

Install required packages:

```r
install.packages(c(
  "httr","jsonlite","dplyr","readr","purrr",
  "stringr","igraph","pROC","optparse"
))
Optional (recommended):

install.packages("renv")
renv::restore()

## How to run
```bash
conda env create -f env/environment.yml
conda activate ppi-rank
jupyter lab

```

Example Results (Alzheimer’s disease)

Using a set of well-established Alzheimer’s genes as seeds, the STRING-based PPI network and RWR-derived features achieve:

AUROC = 0.9785

Recall@50 = 1.0, Recall@100 = 1.0

All held-out Alzheimer’s genes are recovered within the top 50 ranked nodes, indicating that the integrated score (graph centralities + random-walk proximity) is highly effective at prioritizing known disease genes within this PPI subnetwork.

The top 13 ranked genes are exactly the canonical AD genes used as seeds (e.g. APOE, APP, MAPT, TREM2, BIN1, CLU, PICALM, PSEN1/2, ABCA7, SORL1, CD33, PLCG2), followed by biologically plausible neighbors such as SRC, GRB2, C3, SNCA, BACE1, ERBB4 and LRP8. This pattern suggests that the method captures both:

Core AD biology (re-ranking the seed genes to the very top), and network-proximal candidates involved in signaling, immune response and neurodegeneration.

Note that these metrics are likely optimistic, because the evaluation is performed on a STRING network expanded around the seed genes (non-seed nodes are not guaranteed to be true negatives). A more stringent assessment would require an external test set or cross-disease benchmarking.


Please reach out for any collaborations or feedback:
Suprajaa V
suprajaav4@gmail.com
