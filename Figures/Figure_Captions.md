# Main Manuscript Figures

This folder contains figures from the main manuscript.

### Figure_1.png
* Title: Statistical analysis workflow.
* Description: The initial UKB cohort was splitted into a training and a test sets. The Cox regression analysis was done in the training set. For the first model, coefficients of different sets of factors (for different p-value thresholds) were used to compute Environmental Risk Scores (ERS). For the second model, an Elastic Net Regression was applied to significant factors (adjusted p-value < 0.05) for feature selection. The median ERS of the training set was used as a threshold for individual stratification on the test set. 

### Figure_2.png
* Title: Associated lifestyle and environmental factors across ICD-10 disease chapters (males).
* Description: This bubble diagram shows the number of times a factor was found associated with each ICD-10 disease Chapter. The color of the bubble indicates the direction of the effect (red = risk; blue = protective). Factors are grouped by categories. Only those factors associated with at least 10% of the phenotypes of the Chapter are shown.

### Figure_3.png
* Title: Exposome map of significant associations between environmental and lifestyle factors and ICD-10 disease Chapters.
* Description: This chord diagram shows the significant associations that were found in each ICD-10 disease Chapter grouped by lifestyle categories. (A) Exposome map for males. (B) Exposome map for females.

### Figure_4.png
* Title: Manhattan plots of the association between environmental and lifestyle factors and type 2 diabetes in males (A) and in females (B).
* Description: The x axis shows the category domains and the y axis represents statistical significance. The horizontal dashed line indicates p-value threshold after Bonferroni correction.

### Figure_5.png
* Title: ROC-AUC for two models of ERS prediction on T2D in males and females.
* Description: The x axis shows the false positive rate (FPR) and the y axis represents the true positive rate (TPR). (A) Prediction with ERS derived from lifestyle factors coefficients in males (p < 1 x 10-24). (B) Prediction with ERS derived from lifestyle factors coefficients in females (p < 1 x 10-30). (C) Prediction with ERS derived from Elastic Net feature selection in males. (D) Prediction with ERS derived from Elastic Net feature selection in females.
