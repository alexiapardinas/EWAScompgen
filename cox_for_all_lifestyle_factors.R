#!/usr/bin/env Rscript
library(survival, lib.loc = "/gpfs/home/bsc/bsc695622/R/x86_64-pc-linux-gnu-library/4.3/")
library(dplyr, lib.loc = "/gpfs/home/bsc/bsc695622/R/x86_64-pc-linux-gnu-library/4.3/")
library(nnet, lib.loc = "/gpfs/home/bsc/bsc695622/R/x86_64-pc-linux-gnu-library/4.3/")
library(WeightIt, lib.loc = "/gpfs/home/bsc/bsc695622/R/x86_64-pc-linux-gnu-library/4.3/")
library(ggplot2, lib.loc = "/gpfs/home/bsc/bsc695622/R/x86_64-pc-linux-gnu-library/4.3/")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure argument is provided
if(length(args) == 0) {
  stop("Please provide sex, level and code.")
}

sex = args[1]
level = args[2]
ICD10_code = args[3]

# Sanitize the ICD10_code variable function
sanitize_ICD10_code <- function(ICD10_code) {
  gsub("-", "_", ICD10_code)  # Replace hyphens with underscores
}

sanitized_ICD10_code <- sanitize_ICD10_code(ICD10_code)

# Load data
data_name <- paste0("/gpfs/scratch/bsc05/bsc695622/EWAS/survival_dataframes/", sex, "/", level, "_ICD10_", sex, "_survival_train_", ICD10_code, ".tsv")
data <- read.csv(data_name, sep="\t", header=TRUE, fileEncoding="UTF-8")
var_info_file_name <- sprintf("/gpfs/scratch/bsc05/bsc695622/EWAS/actionable_lifestyle_factors/lifestyle_factors_list_%s.txt", sex)
var_info <- read.csv(var_info_file_name, sep="\t", header=TRUE)

# Clean file of missing values and assign rownames
data[data == ""] <- NA
data[data == ".-."] <- NA
rownames(data) <- data$eid
data <- data[, -1]

# Add "X" in front of variable ids of var_type
var_info$id <- paste0("X", var_info$id)

# Define variables
covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "TOWNSEND_22189", "age_recruitment")
event_col <- sprintf("event_%s", sanitized_ICD10_code)
age_diag_col <- sprintf("age_%s", sanitized_ICD10_code)

# Transform negative values from all variables except covariates, event and age diagnosis columns to NA (since UKB coding codes our coding's NAs in negative values)
data[, -which(names(data) %in% c(covariates, event_col, age_diag_col))] <- lapply(data[, -which(names(data) %in% c(covariates, event_col, age_diag_col))], function(x) ifelse(x < 0, NA, x))

# Prepare the survival object
data$SurvObj <- with(data, Surv(time = get(age_diag_col), event = get(event_col)))

# Initialize a results dataframe to store output for each kind of variable
results_bin <- data.frame(
  Variable = character(),
  Type = character(),
  Size = numeric(),
  Cases = numeric(),
  Coefficient = numeric(),
  HazardRatio = numeric(),
  pValue = numeric(),
  SE = numeric(),
  CI95 = character(),
  stringsAsFactors = FALSE)

results_cat <- results_bin
results_quant <- results_bin

# Function to validate variable existence and sufficient categories to analyse
check_variable <- function(var_name, data) {
  if (!(var_name %in% colnames(data))) {
    return(FALSE)
  }
  if (any(is.na(data[[var_name]]))) {
    data <- data %>% filter(!is.na(.data[[var_name]]))
  }
  return(nrow(data) > 0)
}

# Create a list for cox_summaries
cox_summaries <- list()

# Run a cox analysis for all lifestyle factors
for (i in 1:nrow(var_info)) {
  var_name <- var_info$id[i]
  var_type <- var_info$type[i]
  if (!check_variable(var_name, data)) {
    warning(paste("Variable", var_name, "not found in dataset or has insufficient data. Skipping."))
    next
  }
  
  # Create a filtered dataset with individuals with info
  data_filt <- data[complete.cases(data[, c(var_name, covariates)]), ]
  
  # Ensure enough cases exist
  N_cases <- sum(data_filt[[event_col]] == 1)
  if (N_cases < 100) {
    warning(paste("Not enough cases for", var_name, ". Skipping."))
    next
  }

  # Calculate IPW weights based on variable type
  print(paste("Calculating IPW weights for", var_name))
  tryCatch({
    if (var_type == "quant") {
      # Quantitative variable
      # Normalization
      data_filt[[var_name]] <- scale(data_filt[[var_name]], center = TRUE, scale = TRUE)

      weight_model <-  weightit(as.formula(paste(var_name, "~", paste(covariates, collapse = "+"))),
                                data = data_filt,
                                method = "glm",
                                estimand = "ATE",
                                density = "dt_2",
                                link = "identity")
      
      data_filt$ipw_weights <- weight_model$weights
      
    } else if (var_type == "cat") {
      # Categorical variable
      data_filt[[var_name]] <- as.numeric(as.character(data_filt[[var_name]]))
      
      weight_model <-  weightit(as.formula(paste(var_name, "~", paste(covariates, collapse = "+"))),
                                data = data_filt,
                                method = "glm",
                                estimand = "ATE",
                                density = "dt_2",
                                link = "identity")
      
      data_filt$ipw_weights <- weight_model$weights
      
    } else if (var_type == "bin") {
      # Binary variable
      # Be sure that the variable has two possible values (0 and 1)
      if (length(unique(data_filt[[var_name]])) != 2) {
        warning(paste("Binary variable", var_name, "does not have two possible values. Skipping."))
        next
      }
      weight_model <- weightit(as.formula(paste(var_name, "~", paste(covariates, collapse = "+"))),
                               data = data_filt,
                               estimand = "ATE",
                               method = "glm",
                               link = "logit")
      
      data_filt$ipw_weights <- weight_model$weights
      
    } else {
      warning(paste("Unknown variable type for", var_type))
      next
    }
    
    # Create the formula for Cox regression
    formula <- as.formula(paste(
      "SurvObj ~",
      var_name)
    )
    
    # Perform Cox regression
    cox_model <- coxph(formula, data = data_filt, weights = data_filt$ipw_weights)

    # Check PH assumption
    ph_test <- cox.zph(cox_model)
    ph_p <- tryCatch(ph_test$table[var_name, "p"], error = function(e) NA)

    if (!is.na(ph_p) && ph_p < 0.05) {
      print(paste("PH assumption violated for", var_name, "- adding time interaction."))
      cox_model <- coxph(as.formula(paste("SurvObj ~", var_name, "+ tt(", var_name, ")")), 
			 data = data_filt, 
			 tt = function(x, t, ...) x * log(t),
			 weights = data_filt$ipw_weights)
    }
	    
    # Extract results
    summary_model <- summary(cox_model)
    cox_summaries[[var_name]] <- summary_model
    
    # Identify all coefficient rows corresponding to the variable (important for categorical variables)
    coef_rows <- grep(paste0("^", var_name), rownames(summary_model$coefficients), value = TRUE)
    
    if (length(coef_rows) > 0) {
      for (coef_row in coef_rows) {
        coef <- summary_model$coefficients[coef_row, "coef"]
        hr <- summary_model$coefficients[coef_row, "exp(coef)"]
        pval <- summary_model$coefficients[coef_row, "Pr(>|z|)"]
        error <- summary_model$coefficients[coef_row, "se(coef)"]
        lower95 <- summary_model$conf.int[coef_row, "lower .95"]
        upper95 <- summary_model$conf.int[coef_row, "upper .95"]

        # Append results to the corresponding dataframe
        if (var_type == "quant") {
          results_quant <- results_quant %>% add_row(
            Variable = var_name,
            Type = var_type,
            Size = cox_model$n,
            Cases = N_cases,
            Coefficient = coef,
            HazardRatio = hr,
            pValue = pval,
            SE = error,
            CI95 = paste0(round(lower95, 2), " - ", round(upper95, 2)) 
          )
        } else if (var_type == "cat") {
          results_cat <- results_cat %>% add_row(
            Variable = var_name,
            Type = var_type,
            Size = cox_model$n,
            Cases = N_cases,
            Coefficient = coef,
            HazardRatio = hr,
            pValue = pval,
            SE = error,
            CI95 = paste0(round(lower95, 2), " - ", round(upper95, 2))
          )
        } else if (var_type == "bin") {
          results_bin <- results_bin %>% add_row(
            Variable = var_name,
            Type = var_type,
            Size = cox_model$n,
            Cases = N_cases,
            Coefficient = coef,
            HazardRatio = hr,
            pValue = pval,
            SE = error,
            CI95 = paste0(round(lower95, 2), " - ", round(upper95, 2))
          )
        }
      }
    } else {
      warning(paste("Variable", var_name, "not found in Cox regression coefficients. Skipping."))
    }
  }, error = function (e) {
	  warning(paste("Error in processing variable", var_name, ":", e$message))
  })
}

# Create a dataframe classifying the positive and negative significant associations (pvalue < 0.05)
results_bin$Significance <- ifelse(results_bin$pValue < 0.05, ifelse(results_bin$Coefficient > 0, "Adverse", "Protective"), "Non-significant")
results_cat$Significance <- ifelse(results_cat$pValue < 0.05, ifelse(results_cat$Coefficient > 0, "Adverse", "Protective"), "Non-significant")
results_quant$Significance <- ifelse(results_quant$pValue < 0.05, ifelse(results_quant$Coefficient > 0, "Adverse", "Protective"), "Non-significant")

# Save all results dataframes into a same dataframe
results <- rbind(results_bin, results_cat, results_quant)

# Function to apply Bonferroni correction
results$Bonferroni_pValue <- p.adjust(results$pValue, method = "bonferroni")
results$Significance_Bonferroni <- results$Bonferroni_pValue < 0.05

output_file <- paste0("/gpfs/scratch/bsc05/bsc695622/EWAS/cox_results/", sex, "/", level, "_", ICD10_code, "_cox_results_", sex, ".txt")
write.table(results, file = output_file, sep = "\t", row.names = FALSE)

# Save all summaries to a text file
summary_list <- lapply(names(cox_summaries), function(factor_name) {
  c(paste("\n=============================\nSummary for:", factor_name, "\n=============================\n"),
    capture.output(cox_summaries[[factor_name]]))
})
summaries_file <- paste0("/gpfs/scratch/bsc05/bsc695622/EWAS/cox_results/", sex, "/", level, "_", ICD10_code, "_cox_summaries_", sex, ".txt"
)			 

# Save as an RDS file for later analysis
saveRDS(cox_summaries, "cox_summaries.rds")
