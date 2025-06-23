# Supplementary Information Overview

This folder contains supplementary figures and tables referenced in the main manuscript.

==============================
## Supplementary Figures (PNG)
==============================

### Supplementary_Figure_S1.png
* Title: Number of phenotypes per ICD-10 Chapter.
* Description: Bar plot of counts of phenotypes belonging to each ICD-10 Chapter. Phenotypes include both those common between sexes (males and females) and sex-unique.

### Supplementary_Figure_S2.png
* Title: Number of lifestyle and environmental factors grouped by category.
* Description: Bar plot of counts of factors belonging to each lifestyle category. Factors included in each category can be found in Supplementary Table S1.

### Supplementary_Figure_S3.png
* Title: Associated lifestyle and environmental factors across ICD-10 disease chapters (females).
* Description: This bubble diagram shows the number of times a factor was found associated with each ICD-10 disease Chapter. The color of the bubble indicates the direction of the effect (red = risk; blue = protective). Factors are grouped by categories. Only those factors associated with at least 10% of the phenotypes of the Chapter are shown.

### Supplementary_Figure_S4.png
* Title: Lifestyle category enrichment per ICD-10 Chapter.
* Description: This heatmap shows the enrichment score (ES) in log2 scale of the lifestyle categories (diet, nutrients, physical activity, etc.) per ICD-10 Chapter. ES is calculated as the ratio of observed count of associations belonging to each category in diseases of the Chapter and expected ones (if they were distributed randomly and independently).

### Supplementary_Figure_S5.png
* Title: Risk and protective lifestyle factors associated with T2D grouped by lifestyle categories.
* Description: Representation of risk and protective lifestyle categories found associated with T2D, calculated as a proportion within categories (number of associated factors in the category / total number of factors in the category). (A) Males. (B) Females.

### Supplementary_Figure_S6.png
* Title: Sex-differential analysis on type 2 diabetes.
* Description: Representation of significantly different effect sizes associations between males (blue) and females (coral). (A) Pure (sex-specific) significant associations, found significantly associated with T2D in one of the sexes. (B) Quantitative (different effect size magnitude) significant associations, found significantly lassociated with T2D in both sexes but exhibiting different effect sizes towards the same direction (protective or adverse in both sexes).

### Supplementary_Figure_S7.png
* Title: Predictive models' performance comparison.
* Description: Comparison of performance metrics between Model 1 (ERS derived from Cox p-value threshold subsetting) and Model 2 (ERS derived from elastic net feature selection). (A) Males. (B) Females.

==============================
## Supplementary Tables (TSV)
==============================

### Supplementary_Table_S1.tsv
* Title: Actionable Lifestyle Factors.
* Description: List of actionable lifestyle factors tested for association with ICD-10 diseases, including information about UK Biobank ID, type (binary or quantitative), sex in which is present, lifestyle category, etc..

### Supplementary_Table_S2.tsv
* Title: Top 10 risk and protective lifestyle factors associated across ICD-10 disease Chapters.
* Description: List of the top 10 risk and protective lifestyle factors that were found associated more times across all the ICD-10 disease Chapters, including information about the category they belong to, and the number and name of chapters where association was found.

### Supplementary_Table_S3.tsv
* Title: Set of lifestyle factors used as predictors.
* Description: List of the selected lifestyle factors used as predictors in each of the predictive models (Cox p-value threshold subsetting and Elastic Net Regression).

### Supplementary_Table_S4.tsv
* Title: Models' prediction performance metrics and confusion matrices.
* Description: Summary of Cox's ERS and Elastic net's ERS performance descriptive metrics.

### Supplementary_Table_S5.tsv
* Title: Sex-differential analysis significant results.
* Description: List of significantly different factors (effect sizes) between males and females. The table includes the results for both the full evaluation and the first-filtering approaches.
