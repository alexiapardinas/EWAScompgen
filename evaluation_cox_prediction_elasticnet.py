import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, precision_recall_curve, auc, roc_auc_score, confusion_matrix, precision_score, recall_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold

# Get command line arguments
if len(sys.argv) < 4:
    raise ValueError("Please provide sex, level, and ICD10_code as arguments.")

sex = sys.argv[1]
level = sys.argv[2]
ICD10_code = sys.argv[3]

# Load data
results_cox = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/cox_results/{sex}/{level}_{ICD10_code}_cox_results_{sex}.txt", sep="\t")
test_data = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/survival_dataframes/{sex}/{level}_ICD10_{sex}_survival_test_{ICD10_code}.tsv", sep="\t")
train_data = pd.read_csv(f"/gpfs/scratch/bsc05/bsc695622/EWAS/survival_dataframes/{sex}/{level}_ICD10_{sex}_survival_train_{ICD10_code}.tsv", sep="\t")

# Define variables
event_col = f"event_{ICD10_code}".replace("-", "_")
age_diag_col = f"age_{ICD10_code}".replace("-", "_")
exclude_cols = {event_col, age_diag_col}

# Filter significant features from results_cox
results_filtered = results_cox[
    (results_cox["Bonferroni_pValue"] < 0.05) &
    (results_cox["SE"] <= results_cox["Coefficient"].abs()) &
    (results_cox["HazardRatio"].astype(float) >= 0.01) &
    (results_cox["Bonferroni_pValue"] != 0)
]

# Extract feature names from results_cox, removing leading "X"
selected_feature_names = results_filtered["Variable"].str.replace("^X", "", regex=True).tolist()

# Clean test/train data
test_data.replace(["", ".-."], np.nan, inplace=True)
test_data.set_index("eid", inplace=True)
train_data.replace(["", ".-."], np.nan, inplace=True)
train_data.set_index("eid", inplace=True)

# Select features from train data
X_train = train_data[selected_feature_names].apply(pd.to_numeric, errors='coerce')
X_train.dropna(axis=1, thresh=int(0.8 * len(X_train)), inplace=True)  # keep features with enough data
X_train.fillna(X_train.median(), inplace=True)

# Get outcome
y_train = train_data[event_col].astype(int).loc[X_train.index]

# Standardize features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)

# Elastic net logistic regression with cross validation
clf = LogisticRegressionCV(
    cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
    penalty='elasticnet',
    solver='saga',
    l1_ratios=[0.1, 0.5, 0.9],
    scoring='roc_auc',
    max_iter=5000,
    n_jobs=-1,
    verbose=1,
    random_state=42
)
clf.fit(X_train_scaled, y_train)

# Get non-zero coefficients
coef = clf.coef_.flatten()
nonzero_mask = coef != 0
selected_features_final = X_train.columns[nonzero_mask]
selected_coefs = coef[nonzero_mask]

print(f"\nSelected {len(selected_features_final)} features:")
sorted_feats = sorted(zip(selected_features_final, selected_coefs), key=lambda x: abs(x[1]), reverse=True)
for feat, weight in sorted_feats:
    print(f"{feat}: {weight:.4f}")

# Prepare train set
train_data_clean = train_data[selected_features_final].apply(pd.to_numeric, errors='coerce')
train_data_clean.fillna(train_data_clean.median(), inplace=True)
X_train = scaler.transform(train_data_clean)
y_train = train_data[event_col].astype(int).loc[train_data_clean.index]

# Prepare test set
test_data_clean = test_data[selected_features_final].apply(pd.to_numeric, errors='coerce')
test_data_clean.fillna(test_data_clean.median(), inplace=True)
X_test = scaler.transform(test_data_clean)
y_test = test_data[event_col].astype(int).loc[test_data_clean.index]

# Compute risk scores
test_risk_score = X_test @ selected_coefs
test_risk_score = pd.Series(test_risk_score, index=test_data_clean.index, name="Risk_Score")
train_risk_score = X_train @ selected_coefs
train_risk_score = pd.Series(train_risk_score, index=train_data_clean.index, name="Risk Score")

# Threshold tuning
threshold = train_risk_score.median()

# Get predictions
test_data_clean["predicted_class"] = (test_risk_score > threshold).astype(int)

# Evaluation
cm = confusion_matrix(y_test, test_data_clean["predicted_class"])
tn, fp, fn, tp = cm.ravel()

fpr, tpr, _ = roc_curve(y_test, test_risk_score)
roc_auc = auc(fpr, tpr)
accuracy = np.mean(test_data_clean["predicted_class"] == y_test)
accuracy_w = balanced_accuracy_score(y_test, test_data_clean["predicted_class"])
sensitivity = recall_score(y_test, test_data_clean["predicted_class"])
precision = precision_score(y_test, test_data_clean["predicted_class"])
specificity = tn / (tn + fp)
npv = tn / (tn + fn)
bookmakers_informedness = sensitivity + specificity - 1
bias_mk = precision + npv - 1
mcc = matthews_corrcoef(y_test, test_data_clean["predicted_class"])

# Print metrics
print("RESULTS")
print(f"Number of predictors: {len(selected_features_final)}")
print(f"Confusion matrix: {cm}")
print(f"AUC: {roc_auc:.4f}")
print(f"Accuracy_w: {accuracy_w:.4f}")
print(f"Sensitivity: {sensitivity:.4f}")
print(f"Precision: {precision:.4f}")
print(f"BM: {bookmakers_informedness:.4f}")
print(f"Bias: {bias_mk:.4f}")
print(f"MCC: {mcc:.4f}\n")

if sex == 'M':
    sex_text = 'Males'
elif sex == 'F':
    sex_text = 'Females'

# Plot ROC curve
plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, label=f"AUC={roc_auc:.4f}")
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"ROC Curve for Elastic Net on {ICD10_code} diagnosis prediction in {sex_text}")
plt.legend()
plt.show()

# Plot Precision-Recall curve
precision_vals, recall_vals, _ = precision_recall_curve(y_test, test_risk_score)
pr_auc = auc(recall_vals, precision_vals)
plt.figure(figsize=(8, 6))
plt.plot(recall_vals, precision_vals, label=f"AUC={pr_auc:.4f}")
plt.plot([0, 1], [1, 0], linestyle='--', color='gray')
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title(f"Precision-Recall Curve for Elastic Net prediction on {ICD10_code} diagnosis in {sex_text}")
plt.legend()
plt.show()
