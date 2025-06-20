# CES_Discovery
# Reproducibility Repository for “Cohort-level analysis of human de novo mutations points to drivers of clonal expansion in spermatogonia”

This repository contains Jupyter notebooks to reproduce key analyses, statistical tests, and visualizations from our manuscript. The relevant input data files required to run these notebooks are available via Zenodo:  
**https://zenodo.org/records/15660433)**

---

## Jupyter Notebooks

### 1️⃣ **Enrichment_plots.ipynb**  
Implements key visualizations and analyses reported in the manuscript:
- Statistical tests for homogeneity of de novo mutation sampling across cohorts  
- Enrichment and recurrence analyses of loss-of-function (LoF) and gain-of-function (GoF) variants  
- Investigation of the "winner’s curse" effect  
- Plots of observed vs. expected LoF counts  
- Relationship between gene constraint (LOEUF) and mutation burden  
- CES effects of LoF-2 genes across cohorts  

**Dependencies:**  
- `pandas`  
- `numpy`  
- `seaborn`  
- `gzip`  
- `scipy.optimize` (for `minimize`)  
- `scipy.special` (for `gammaln`)  
- `scipy.stats` (for `multinomial`, `nbinom`, `poisson`, `gamma`)  
- `matplotlib.pyplot`  

---

### 2️⃣ **paternal_bias.ipynb**  
Contains code for calculating the power of the paternal overtransmission test for detecting clonal expansions (CES) in spermatogonia:
- Models power in relation to number of tested mutations *N_vars* and mutation rate elevation *κ*  
- Assumes binomial test with *p₀ = 3/4* under the null  
- Models CES effect as: *p_alt = κ / (1/3 + κ)*  

**Dependencies:**  
- `numpy`  
- `matplotlib`  
- `scipy`  
- `seaborn`  

---

### 3️⃣ **Expression_sc.ipynb**  
Processes Human Protein Atlas single-cell RNAseq data and generates auxiliary plots and differential expression tests:
- Barplots of CES candidate gene expression in spermatogonia (average nTPM and relative nTPM)  
- Barplots for oogonia  
- Heatmap of CES candidate gene expression in fetal DAG clusters  

**Dependencies:**  
- `pandas`  
- `numpy`  
- `matplotlib`  
- `seaborn`  
- `scipy`  

---

### 4️⃣ **Roulette_QC.ipynb**  
Performs quality control of the mutation rate model (Roulette):
- Per-site observed vs. expected synonymous mutation counts  
- Poisson vs. Negative Binomial variance inflation  
- Per-gene variance inflation analyses  
- Model fit diagnostics and overdispersion tests  

**Dependencies:**  
- `pandas`  
- `numpy`  
- `matplotlib`  
- `scipy`  
- `gzip`  

---

### 5️⃣ **CES_SFS.ipynb**  
Implements a population-genetic test for validating CES candidates:
- Based on mean and variance estimates of the site frequency spectrum (SFS)

**Dependencies:**  
- `pandas`  
- `numpy`  
- `matplotlib`  
- `scipy`  
- `gzip`  

---


- Please see the **Methods** section of the paper for details on models, tests, and assumptions.  


