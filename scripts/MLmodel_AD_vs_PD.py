# Step 1: Segregate the samples into healthy and diseased and remove the metadata 


import pandas as pd

# File path
file_path = r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\PD\csv_files\GSE57475_series_matrix.csv"

# Load gene expression data
df = pd.read_csv(file_path, sep=",", comment="!", header=1, index_col=0)  # Set 'IDREF' as index

# Load metadata (First 105 rows contain metadata)
metadata_df = pd.read_csv(file_path, sep=",", header=None, nrows=105)

# Extract relevant metadata (Ensure the correct row indices)
patient_status_row = metadata_df.iloc[[60, 38]]  # Adjust if necessary

# Transposing & cleaning metadata
patient_status_df = patient_status_row.T.iloc[1:].reset_index(drop=True)

# Ensure correct column naming
patient_status_df.columns = ["PatientID", "DiseaseStatus"]

# Set PatientID as index
patient_status_df = patient_status_df.set_index("PatientID")

# Transpose gene expression data so IDREFs remain in rows
gene_df = df.T  # Now, IDREFs are rows, and PatientIDs are columns

# Merge expression data with metadata
merge_df = gene_df.merge(patient_status_df, left_index=True, right_index=True)

# Function to classify samples
def classify_sample(description):
    description = str(description).lower()  # Normalize text
    if "healthy" in description or "control" in description:
        return "HEALTHY"
    elif "mci" in description:
        return "MCI"
    elif "pd" in description:
        return "PD"
    else:
        return "Unknown"


# Apply classification
merge_df["Category"] = merge_df["DiseaseStatus"].apply(classify_sample)

# Segregate based on category
healthy_samples = merge_df[merge_df["Category"] == "HEALTHY"]
mci_samples = merge_df[merge_df["Category"] == "MCI"]
ad_samples = merge_df[merge_df["Category"] == "PD"]

combined_samples = merge_df[merge_df["Category"].isin(["PD", "HEALTHY"])]

labels = combined_samples["Category"]

new_column_names = [
    f"{sample_id}_{labels[sample_id]}"
    for sample_id in combined_samples.index]


combined_samples_T = combined_samples.T

#print(new_column_names)
#combined_samples = new_column_names.T

# Assign new column names
combined_samples_T.columns = new_column_names

# Save the labeled data
combined_samples_T.to_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\PD\csv_files\GSE57475_labelled.csv")

#print(combined_samples_T.head())





# Step 2: Find the maximum minimum values, group duplicated genes


import pandas as pd

# Load the GEO matrix table
def load_geo_data(file_path):
    """Load the GEO dataset, remove metadata, and convert to numeric."""
    df = pd.read_csv(file_path, low_memory=False)
    #df = df.iloc[:, 1:]
    df.set_index(df.columns[0], inplace=True)  # Set first column as index (Gene/Probe IDs)
    df = df.apply(pd.to_numeric, errors="coerce")# Convert all values to numeric 
    return df

# Find maximum and minimum values
def find_min_max(df):
    """Find the max and min expression values in the dataset."""
    max_value = df.max().max()  # Maximum expression value
    min_value = df.min().min()  # Minimum expression value
    return max_value, min_value

# Group duplicated gene IDs
def group_duplicate_genes(df):
    """If multiple probes map to the same gene, group by averaging their values."""
    df_grouped = df.groupby(df.index).mean()  # Average duplicate rows (genes)
    return df_grouped


# Main function
def process_geo_data(file_path, index = True):
    """Complete pipeline for processing GEO data."""
    print("Loading data...")
    df = load_geo_data(file_path)
    #df_selected = df.iloc[:, 2:] 

    print("Finding min and max values...")
    max_value, min_value = find_min_max(df)
    print(f"Maximum Value: {max_value}, Minimum Value: {min_value}")

    print("Grouping duplicate genes...")
    df = group_duplicate_genes(df)

    # Save the full dataset (including first two columns) to a new CSV file
    #df.to_csv("output_file", index=False)

    print("Processing completed successfully!")
    return df


# Example
df = process_geo_data(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\24.04.2025\AD_PD_merged.csv")



#Log transformations

import pandas as pd
import numpy as np 


# Apply log2 transformation (log2(x + 1) to avoid log(0))
df1_log = np.log2(df1 + 1)

# Save the transformed dataset
df1_log.to_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\PD\csv_files\log_transformed\GSE6613.csv")

print("Log2 transformation applied successfully!")

print("Max value after log transformation")
max_value = df1_log.max().max()
print(max_value)

print("Min value after log transformation")
min_value = df1_log.min().min()
print(min_value)



#Merging two GEO datasets together into one file -  common genes only
import pandas as pd
import glob
import os

# Define file paths (update these with actual paths)
diseased_files = glob.glob(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_1\Merged_AD_files\*.csv")
healthy_files = glob.glob(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\PD\Processed_PD_data\Merged_PD_files\*.csv")

def load_geo_data(file_paths):
    """
    Reads and processes multiple GEO files, setting 'Gene Symbol' as index.
    Returns a merged DataFrame for all files in file_paths.
    """
    dataframes = []

    for file in file_paths:
        df = pd.read_csv(file, dtype=str, low_memory=False) 
        df.columns = df.columns.str.strip().str.lower()  

        # Use 'Gene Symbol' or 'Probe ID' as index
        if "gene_name" in df.columns:
            df.set_index("gene_name", inplace=True)
        else:
            raise ValueError(f"Neither 'Gene Symbol' nor 'Probe ID' found in {file}")

        df = df.apply(pd.to_numeric, errors='coerce')  # Convert values to numeric
        df = df.groupby(df.index).mean()  # Handle duplicate gene symbols by averaging

        dataset_name = os.path.splitext(os.path.basename(file))[0]
        #df = df.add_suffix(f"_{dataset_name}")  # Add dataset name as suffix

        dataframes.append(df)

    if dataframes:
        merged_df = pd.concat(dataframes, axis=1, join="inner")  # Outer join to keep all genes
        merged_df.fillna(0, inplace=True)  # Fill missing values with 0
        return merged_df
    else:
        return pd.DataFrame()  # Return empty DataFrame if no files found

# Load diseased and healthy datasets
diseased_df = load_geo_data(diseased_files)
healthy_df = load_geo_data(healthy_files)

# Ensure all genes are included
all_genes = sorted(set(diseased_df.index).union(set(healthy_df.index)))

# Reindex both dataframes to have the same gene index
diseased_df = diseased_df.reindex(all_genes).fillna(0)
healthy_df = healthy_df.reindex(all_genes).fillna(0)

# Merge horizontally (side-by-side) so diseased and healthy samples stay separate
final_merged_df = pd.concat([diseased_df, healthy_df], axis=1)

# Save the final merged CSV file
output_path = r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_merged.csv"
final_merged_df.to_csv(output_path)

print(f" Merging complete! Final dataset saved as: {output_path}")



# Combined: Total of NaNs + Zeros
na_total = df.isna().sum().sum()
zero_total = (df == 0).sum().sum()

print(f"NaNs: {na_total}, Zeros: {zero_total}")

df_imputed = df.copy()

# Replace 0 with row-wise mean (excluding zeros and NaNs)
for idx in df_imputed.index:
    row = df_imputed.loc[idx]
    non_zero_mean = row[(row != 0) & (~row.isna())].mean()
    df_imputed.loc[idx] = row.replace(0, non_zero_mean)




#Quantile normalization

import pandas as pd
import scanpy as sc
from sklearn.preprocessing import quantile_transform

# Load your merged GEO data
df = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_imputed.csv", index_col=0)

# -------- Step 1: Quantile Normalization --------
qn_data = quantile_transform(df.values, axis=0, n_quantiles=100, output_distribution='normal', copy=True)
df_qn = pd.DataFrame(qn_data, index=df.index, columns=df.columns)


#df_qn.to_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_quantile.csv")
print(" Quantile normalization complete and saved.")


# Histogram for quantile normalization

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Histogram before normalization
plt.subplot(1, 2, 1)
sns.histplot(df_imputed.values.flatten(), bins=100, kde=True, color='salmon')
plt.title("Before Quantile Normalization")
plt.xlabel("Expression Value")
plt.ylabel("Frequency")

# Histogram after normalization
plt.subplot(1, 2, 2)
sns.histplot(df_qn.values.flatten(), bins=100, kde=True, color='skyblue')
plt.title("After Quantile Normalization")
plt.xlabel("Expression Value")
plt.ylabel("Frequency")

plt.tight_layout()
plt.show()


#Batch effect correction

import pandas as pd
import numpy as np
import scanpy as sc
import re

# Load quantile-normalized dataframe
df_qn = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_quantile.csv", index_col=0)

# Extract metadata: batch = gse ID, condition = disease type
samples = df_qn.columns.to_series()
batch = samples.str.extract(r'_?(gse\d+)_', flags=re.IGNORECASE)[0]
condition = samples.str.extract(r'_gse\d+_(\w+)_', flags=re.IGNORECASE)[0]

# Transpose for scanpy (samples = rows)
df_T = df_qn.T.copy()
df_T['batch'] = batch.values
df_T['condition'] = condition.values

# Drop samples with missing metadata
df_T = df_T.dropna(subset=['batch', 'condition'])

# Separate numeric expression data only
expression_data = df_T.drop(columns=['batch', 'condition']).copy()

# expression_data = expression_data.fillna(expression_data.mean())
# Convert all values to numeric, coercing errors to NaN

expression_data = expression_data.apply(pd.to_numeric, errors='coerce')

# Fill missing values with column-wise mean (gene-wise)
#expression_data = expression_data.fillna(expression_data.mean())


# Align metadata with filtered expression_data
filtered_df_T = df_T.loc[expression_data.index]
filtered_batch = filtered_df_T['batch'].astype(str)

# Confirm dimensions before AnnData
print("Expression data shape:", expression_data.shape)
print("Batch labels:", filtered_batch.shape)

# Create AnnData object
adata = sc.AnnData(expression_data)
adata.obs['batch'] = filtered_batch.values

# Run ComBat for batch correction
sc.pp.combat(adata, key='batch')

# Convert result back to gene x sample format
corrected_df = pd.DataFrame(adata.X.T, index=df_qn.index, columns=expression_data.index)

# Save corrected data
#corrected_df.to_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_combat_corrected.csv")

print(" Batch correction complete and saved.")

# Histogram

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming df_corrected is your batch-corrected DataFrame (genes x samples)

# Flatten all values into a 1D array
allvalues = corrected_df.values.flatten()

# Plot histogram or KDE
plt.figure(figsize=(10, 5))
sns.histplot(allvalues, bins=100, kde=True, color='skyblue')
plt.title("Histogram of All Expression Values (After Batch Correction)")
plt.xlabel("Expression Value")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()

#PCA plot - Before and after batch correction

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def plot_pca(df, title):
    from matplotlib.cm import get_cmap

    # Extract metadata again
    samples = df.columns.to_series()
    batch = samples.str.extract(r'_?(gse\d+)_', flags=re.IGNORECASE)[0]
    condition = samples.str.extract(r'_gse\d+_(\w+)_', flags=re.IGNORECASE)[0]

    # Run PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df.T)

    # Create PCA dataframe
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df['batch'] = batch.values
    pca_df['condition'] = condition.values

    # Plot
    plt.figure(figsize=(10, 7))
    unique_conditions = pca_df['condition'].unique()
    
    for cond in unique_conditions:
        subset = pca_df[pca_df['condition'] == cond]
        plt.scatter(subset["PC1"], subset["PC2"], label=cond, alpha=0.7)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    plt.legend(title="Condition")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    
plot_pca(df_qn, "PCA Before Batch Correction")
plot_pca(corrected_df, "PCA After Batch Correction")




# Label healthy and diseased using labelencoder

import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score

# Load your gene expression data
df = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\AD_PD_combat_corrected.csv", index_col=0)

# Extract labels from the column names
original_columns = df.columns
sample_ids = [col.split()[0] for col in original_columns]  # clean sample names
labels = []
for col in original_columns:
    col_lower = col.lower()
    if "ad" in col_lower:
        labels.append("ad")
    elif "pd" in col_lower:
        labels.append("pd")
    else:
        labels.append("healthy")

# Update column names to just sample IDs
df.columns = sample_ids

# Transpose data so samples are rows
df_transposed = df.T

#  Add the label column
df_transposed['label'] = labels

#  Prepare data for ML
X = df_transposed.drop(['label'], axis=1)
y = df_transposed['label']

#  Encode labels
le = LabelEncoder()
y_encoded = le.fit_transform(y)
print("Label encoding map:", dict(zip(le.classes_, le.transform(le.classes_))))

print(df_transposed['label'].value_counts())


#Machine learning models

# Train-test split with stratification
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
    X, y_encoded, test_size=0.2, stratify=y, random_state=42)

# Scale the features
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train classifiers
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

svm = SVC(kernel="rbf", probability=True, random_state=42)
rf = RandomForestClassifier(n_estimators=100, random_state=42)
xgb = XGBClassifier(use_label_encoder=False, eval_metric="mlogloss", random_state=42)

svm.fit(X_train_scaled, y_train)
rf.fit(X_train_scaled, y_train)
xgb.fit(X_train_scaled, y_train)

y_pred_svm = svm.predict(X_test_scaled)
y_pred_rf = rf.predict(X_test_scaled)
y_pred_xgb = xgb.predict(X_test_scaled)

print(f"SVM Accuracy: {accuracy_score(y_test, y_pred_svm):.4f}")
print(f"Random Forest Accuracy: {accuracy_score(y_test, y_pred_rf):.4f}")
print(f"XGBoost Accuracy: {accuracy_score(y_test, y_pred_xgb):.4f}")



# Confusion Matrices
import seaborn as sns
import matplotlib.pyplot as plt

def plot_confusion_matrix(y_true, y_pred, title):
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(5, 4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=le.classes_, yticklabels=le.classes_)
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title(title)
    plt.show()

plot_confusion_matrix(y_test, y_pred_svm, "SVM")
plot_confusion_matrix(y_test, y_pred_rf, "Random Forest")
plot_confusion_matrix(y_test, y_pred_xgb, "XGBoost")

# Classification Reports
models = {"SVM": svm, "Random Forest": rf, "XGBoost": xgb}
for name, model in models.items():
    print(f"\n{name} Classification Report:")
    print(classification_report(y_test, model.predict(X_test_scaled), target_names=le.classes_))

from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np

# Assuming y_test has 3 classes: 0 = AD, 1 = Healthy, 2 = PD
# Binarize the labels for ROC
classes = [0, 1, 2]
y_test_bin = label_binarize(y_test, classes=classes)

# Get predicted probabilities from your model
y_score = xgb.predict_proba(X_test_scaled)  # Replace lgbm with your model name



# ROC curve for each class
plt.figure(figsize=(8, 6))
for i in range(len(classes)):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_score[:, i])
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=2, label=f"Class {i} (AUC = {roc_auc:.2f})")

plt.plot([0, 1], [0, 1], "k--", lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Multiclass ROC Curve")
plt.legend(loc="lower right")
plt.grid(True)
plt.show()


#SHAP analysis

import pandas as pd
import shap

# Assuming X_train_scaled is a DataFrame
explainer = shap.TreeExplainer(lgb_model)
shap_values = explainer.shap_values(X_train_scaled)

# Save average SHAP values per feature, per class

for i, class_vals in enumerate(shap_values):
    #print(f"Class {i} SHAP shape: {class_vals.shape}")
    
    # Ensure correct shape: (n_samples, n_features)
    if class_vals.shape[1] != len(X.columns):
        #print("Transposing SHAP values...")
        class_vals = class_vals.T

    mean_abs_shap_vals = np.abs(class_vals).mean(axis=0)  # mean SHAP per feature
    mean_abs_shap = pd.DataFrame({
        "feature": X.columns,
        "mean_abs_shap": mean_abs_shap_vals
    }).sort_values(by="mean_abs_shap", ascending=False)

    mean_abs_shap.to_csv(f"mean_abs_shap_class_{i}.csv", index=False)
    
    
# Convert X_test_scaled (numpy) to DataFrame for SHAP plotting
X_test_df = pd.DataFrame(X_test_scaled, columns=X.columns)

# Loop over classes and show SHAP summary per class
for i in range(shap_values.shape[2]):  # loop over 3 classes
    print(f"🔍 SHAP Summary for Class {i} ({le.classes_[i]})")
    shap.summary_plot(shap_values[:, :, i], X_test_df, feature_names=X.columns)
                  
				  
				  

"""

import pandas as pd

# Load DEG comparison files
deg_ad_healthy = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\DEG_AD_vs_Healthy_significant_strict.csv")
deg_pd_healthy = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\DEG_PD_vs_Healthy_significant_strict.csv")
deg_ad_pd = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\DEG_AD_vs_PD_significant_strict.csv")


# Filter significant DEGs
sig_deg_ad_healthy = set(deg_ad_healthy[deg_ad_healthy['adj.P.Val'] < 0.05]['gene'])
sig_deg_pd_healthy = set(deg_pd_healthy[deg_pd_healthy['adj.P.Val'] < 0.05]['gene'])
sig_deg_ad_pd = set(deg_ad_pd[deg_ad_pd['adj.P.Val'] < 0.05]['gene'])


# Load SHAP top genes per class
shap_ad = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\top_genes_ad.csv")['gene']
shap_pd = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\top_genes_finalpd.csv")['gene']
shap_healthy = pd.read_csv(r"C:\Users\tejes\OneDrive\Documents\BIOINFORMATICS\Msc_Bioinformatics_project\Data\top_genes_healthy.csv")['gene']


# Intersections
overlap_ad = list(set(shap_ad).intersection(sig_deg_ad_healthy))
overlap_pd = list(set(shap_pd).intersection(sig_deg_pd_healthy))
overlap_ad_vs_pd = list(set(shap_ad).intersection(sig_deg_ad_pd)) + \
                   list(set(shap_pd).intersection(sig_deg_ad_pd))
				   
"""

#Reactome_2022


enr = gp.enrichr(gene_list=overlap_ad,
           gene_sets='Reactome_2022',
           organism='Human',
           outdir='enrich_AD',
           cutoff=0.1)  # <-- changed from 0.05
				  