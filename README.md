# EWAScompgen
Environment-Wide Association Study on > 3500 ICD-10 diseases using UK Biobank data

In this repository, on the one hand, you can find the scripts that were used to run the EWAS analysis. However, the data they use comes from UK Biobank and is not open access (only available upon request).
* cox_for_all_lifestyle_factors.R: Runs Cox regression analysis for a certain ICD-10 disease (outcome) and all actionable lifestyle factors included in the analysis. 
* evaluation_cox_prediction.py: For different sets of lifestyle factors (predictors), obtained using different significance p-value thresholds, it calculates a linear combination of their effects and values to build Environmental Risk Scores (ERS) and evaluates the prediction of these scores.
* evaluation_cox_prediction_elasticnet.py: It applies an Elastic Net Regression on the significant lifestyle factors (p_adjusted < 0.05), and uses the non-zero coefficients to build Environmental Risk Scores (ERS) and evaluates the prediction of these scores.

On the other hand, you can find Supplementary Material of the corresponding master's thesis.
