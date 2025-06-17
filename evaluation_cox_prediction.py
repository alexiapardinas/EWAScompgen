import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, roc_auc_score, confusion_matrix, precision_score, recall_score, precision_recall_curve, matthews_corrcoef

# Get command line arguments
if len(sys.argv) < 3:
    raise ValueError("Please provide sex to be tested, level and disease code.")

sex = sys.argv[1]
level = sys.argv[2]
ICD10_code = sys.argv[3]

# Load data
results_cox = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/cox_results/{sex}/{level}_{ICD10_code}_cox_results_{sex}.txt", sep="\t")
test_data = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/survival_dataframes/{sex}/{level}_ICD10_{sex}_survival_test_{ICD10_code}.tsv", sep="\t")
train_data = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/survival_dataframes/{sex}/{level}_ICD10_{sex}_survival_train_{ICD10_code}.tsv", sep="\t")

# Define variables
covariates = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "TOWNSEND_22189", "age_recruitment"]
event_col = f"event_{ICD10_code}".replace("-", "_")
age_diag_col = f"age_{ICD10_code}".replace("-", "_")
exclude_cols = {event_col, age_diag_col}

# Data Cleaning
test_data.replace(["", ".-."], np.nan, inplace=True)
test_data.set_index("eid", inplace=True)
train_data.replace(["", ".-."], np.nan, inplace=True)
train_data.set_index("eid", inplace=True)

# Remove PC1-PC7 and TOWNSEND_22189
test_data.drop(columns=[col for col in test_data if "PC" in col or "TOWNSEND_22189" in col], inplace=True, errors='ignore')
train_data.drop(columns=[col for col in train_data if "PC" in col or "TOWNSEND_22189" in col], inplace=True, errors='ignore')

# Open a file to write the results
with open(f"/gpfs/scratch/bsc05/bsc695622/EWAS/evaluation_results/{sex}/{level}_{ICD10_code}_test_prediction_{sex}.txt", "a") as f:
    # Define p-value thresholds
    p_values = [10 ** -i for i in range(2, 31, 2)]
    roc_curves = {}
    metrics_list = {}
    risk_scores_by_p = {}

    for p in p_values:
        print(f"Processing for p-value threshold: {p}")
        results_filtered = results_cox[(results_cox["Bonferroni_pValue"] < p) & (results_cox["SE"] <= results_cox["Coefficient"].abs()) & (results_cox["HazardRatio"].astype(float) >= 0.01) & (results_cox["Bonferroni_pValue"] != 0)]
        beta_values = results_filtered.set_index("Variable")["Coefficient"].to_dict()
    
        train_data_copy = train_data.copy()
        test_data_copy = test_data.copy()

        train_data_copy = train_data_copy.rename(columns=lambda col: f'X{col}' if col not in exclude_cols else col)
        test_data_copy = test_data_copy.rename(columns=lambda col: f'X{col}' if col not in exclude_cols else col)
        
        common_vars_train = list(set(train_data_copy.columns) & set(beta_values.keys()))
        common_vars_test = list(set(test_data_copy.columns) & set(beta_values.keys()))
        
        for var in common_vars_train:
            train_data_copy[var].fillna(train_data_copy[var].median(), inplace=True)
        for var in common_vars_test:
            test_data_copy[var].fillna(test_data_copy[var].median(), inplace=True)

        beta_values_numeric = np.array([beta_values[var] for var in common_vars_train])
        train_data_copy[common_vars_train] = train_data_copy[common_vars_train].apply(pd.to_numeric, errors='coerce')
        test_data_copy[common_vars_test] = test_data_copy[common_vars_test].apply(pd.to_numeric, errors='coerce')
    
        risk_col = f"Risk_Score_p{p}"
        train_data_copy[risk_col] = train_data_copy[common_vars_train].dot(beta_values_numeric)
        test_data_copy[risk_col] = test_data_copy[common_vars_test].dot(beta_values_numeric)
    
        risk_scores_by_p[p] = train_data_copy[[risk_col]].copy()

        test_data_clean = test_data_copy.dropna(subset=[risk_col, event_col])
        train_data_clean = train_data_copy.dropna(subset=[risk_col, event_col])

        X_train = train_data_clean[[risk_col]]
        y_train = train_data_clean[event_col].astype(int)
    
        X_test = test_data_clean[[risk_col]]
        y_test = test_data_clean[event_col].astype(int)
    
        # Threshold tuning
        threshold = X_train.median()
        test_data_clean["predicted_class"] = (X_test > threshold).astype(int)
	
	# Evaluation
        cm = confusion_matrix(y_test, test_data_clean["predicted_class"])
    	tn, fp, fn, tp = cm.ravel()

    	fpr, tpr, _ = roc_curve(y_test, X_test[risk_col])
    	roc_auc = auc(fpr, tpr)
    	roc_curves[p] = (fpr, tpr, roc_auc)
    	accuracy_w = balanced_accuracy_score(y_test, test_data_clean["predicted_class"])
    	sensitivity = recall_score(y_test, test_data_clean["predicted_class"])
    	precision = precision_score(y_test, test_data_clean["predicted_class"])
    	specificity = tn / (tn + fp)
    	npv = tn / (tn + fn)
    	f1 = f1_score(y_test, test_data_clean["predicted_class"])
    	bookmakers_informedness = sensitivity + specificity - 1
    	bias_mk = precision + npv - 1
    	mcc = matthews_corrcoef(y_test, test_data_clean["predicted_class"])
        
        metrics_list[p] = {"Accuracy_w": accuracy_w, "Sensitivity": sensitivity, "Precision": precision, "Specificity": specificity, "BM": bookmakers_informedness, "AUC": roc_auc, "Bias": bias_mk, "MCC": mcc}
    
    
        print(f"RESULTS FOR P-VALUE {p}", file=f)
        print(f"Number of predictors: {len(beta_values_numeric)}", file=f)
        print(f"Confusion matrix: {cm}", file=f)
        print(f"AUC: {roc_auc:.4f}", file=f)
        print(f"Accuracy_w: {accuracy_w:.4f}", file=f)
        print(f"Sensitivity: {sensitivity:.4f}", file=f)
        print(f"Precision: {precision:.4f}", file=f)
        print(f"BM: {bookmakers_informedness:.4f}", file=f)
        print(f"Bias: {bias_mk:.4f}", file=f)
        print(f"MCC: {mcc:.4f}", file=f)
        print(f"F1: {f1:.4f}\n", file=f)

if sex == "M":
    sex_text = "Males"
elif sex == "F":
    sex_text = "Females"

# Plot ROC Curves
plt.figure(figsize=(8, 6))
for p, (fpr, tpr, roc_auc) in roc_curves.items():
    plt.plot(fpr, tpr, label=f"p={p}, AUC={roc_auc:.4f}")
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"ROC Curves for different p-value thresholds on {ICD10_code} diagnosis prediction in {sex_text}")
plt.legend()
plt.show()

