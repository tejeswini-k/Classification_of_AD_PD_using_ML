# Classification of Alzheimer’s and Parkinson’s Disease Using Machine Learning

## Overview
This project focuses on the classification of **Alzheimer’s Disease (AD)** and **Parkinson’s Disease (PD)** using **gene expression data** obtained from public GEO datasets.  
A complete end-to-end pipeline is implemented, covering **data preprocessing, differential expression analysis, machine learning modeling, and model interpretability using SHAP**.

The goal is to identify disease-specific molecular biomarkers and evaluate multiple machine learning algorithms in distinguishing AD and PD samples from healthy controls and from each other.

---

## Dataset
- Source: NCBI Gene Expression Omnibus (GEO)
- Data type: Microarray gene expression
- Classes:
  - Alzheimer’s Disease (AD)
  - Parkinson’s Disease (PD)
  - Healthy Controls
- Samples from multiple GEO studies were merged to increase statistical power.

---

## Data Preprocessing
The following preprocessing steps were performed:

1. **Raw Data Cleaning**
   - Removal of invalid probes and missing values
   - Gene symbol mapping and duplicate handling

2. **Log Transformation**
   - Log2 transformation applied to stabilize variance

3. **Normalization**
   - Quantile normalization to ensure comparable distributions

4. **Batch Effect Correction**
   - Applied **ComBat** to remove inter-dataset batch effects

5. **Filtering Low-Expression Genes**
   - Very low expressed genes removed
   - Threshold-based filtering to reduce noise

---

## 📊 Differential Expression Analysis
- Tool used: **limma**
- Comparisons performed:
  - AD vs Healthy
  - PD vs Healthy
  - AD vs PD
- Differentially expressed genes (DEGs) were used for downstream feature selection and model interpretation.

---

## Machine Learning Models
The following models were implemented and evaluated:

- **Support Vector Machine (SVM)**
- **Random Forest (RF)**
- **LightGBM**
- **XGBoost**

### Classification Tasks
- Binary classification:
  - AD vs Rest
  - PD vs Rest
- Multiclass classification:
  - AD vs PD vs Healthy

---

## Model Training & Evaluation
- Train–test split
- Evaluation metrics:
  - Accuracy
  - Precision
  - Recall
  - F1-score
  - ROC-AUC (binary and multiclass)

---

## Model Interpretability (SHAP)
To interpret model predictions and identify key biomarkers:

- **SHAP values** computed for:
  - Binary classifiers
  - Multiclass models (class-wise SHAP)
- Top contributing genes identified for:
  - AD
  - PD
  - Healthy controls
- SHAP-based feature importance used for feature selection refinement.
