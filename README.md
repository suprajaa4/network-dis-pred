# network-dis-pred
Network based Disease Gene Prediction


 Rank candidate genes for a chosen disease by leveraging PPI topology, random-walk proximity, and node embeddings, with validation against known disease genes.


## Data
- STRING PPI (filtered by confidence)
- DisGeNET disease–gene associations (seed/labels)
- Optional: GEO/GTEx features


## Pipeline
1. Download & harmonize IDs
2. Build induced disease subgraph
3. Compute features (centralities, RWR, embeddings)
4. Aggregate ranks and score candidates
5. Validate with held-out positives (AUROC/AUPRC)


## How to run
```bash
conda env create -f env/environment.yml
conda activate ppi-rank
jupyter lab
